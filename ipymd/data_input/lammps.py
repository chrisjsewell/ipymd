# -*- coding: utf-8 -*-
"""
Created on Mon May 16 01:15:56 2016

@author: cjs14
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""
import numpy as np
import pandas as pd
import os
import glob
import re

from .base import DataInput
        
class LAMMPS_Input(DataInput):
    """ file format according to http://lammps.sandia.gov/doc/read_data.html """
    def setup_data(self, atom_path='',atom_style='atomic'):
        """ get data from file
        
        Parameters
        ----------
        atom_style : 'atomic', 'charge'
            defines how atomic data is listed:
            atomic; atom-ID atom-type x y z
            charge; atom-ID atom-type q x y z
        
        """
        #TODO add more atom styles
        assert atom_style in ['atomic', 'charge']
        assert os.path.exists(atom_path) or not atom_path, 'atom_path does not exist'
        self._atom_path = atom_path
        self._atom_style = atom_style
        self._data_set = True
        
    def _get_atom_data(self,step):        
        """ """
        num_atoms = None
        atom_data = []
        
        with open(self._atom_path, 'r') as f:
            for line in f:
                if len(line.split()) == 0: continue

                if line[0] == '#':
                    continue
                #TODO read massess, etc
                if line.split()[0] == 'Atoms':
                    if num_atoms is None: 
                        raise IOError('the file has not specified how many atoms it has')
                    line = self._skiplines(f,2) # skip blank line
                    for i in range(num_atoms):
                        if self._atom_style == 'atomic':
                            aid, atype, x, y, z = (int(line.split()[0]), 
                                                   int(line.split()[1]), 
                                                   float(line.split()[2]),
                                                   float(line.split()[3]),
                                                   float(line.split()[4]))
                            atom_data.append([aid, atype, x, y, z])
                        elif self._atom_style == 'charge':
                            aid, atype, q, x, y, z = (int(line.split()[0]), 
                                                   int(line.split()[1]), 
                                                   float(line.split()[2]),
                                                   float(line.split()[3]),
                                                   float(line.split()[4]),
                                                   float(line.split()[5]))
                            atom_data.append([aid, atype, q, x, y, z])
                        if i != num_atoms-1:
                            line = self._skiplines(f,1)
                        continue
        
                if len(line.split()) < 2: continue

                if line.split()[1] == 'atoms':
                    num_atoms = int(line.split()[0])
                    continue

        if self._atom_style == 'atomic':
            atom_df = pd.DataFrame(atom_data,columns=['id', 'type', 'x', 'y', 'z'])
        elif self._atom_style == 'charge':
            atom_df = pd.DataFrame(atom_data,columns=['id', 'type', 'q', 'x', 'y', 'z'])

        self._add_colors(atom_df)
        self._add_radii(atom_df)
            
        return atom_df
    
    def _get_meta_data(self, step):
        """ return pandas.Series of origin, a, b & c coordinates """
        xy, xz, yz = 0., 0., 0.
        with open(self._atom_path, 'r') as f:
            for line in f:
                if len(line.split()) == 0: continue

                if line[0] == '#':
                    continue

                if len(line.split()) < 3: continue

                if line.split()[2] == 'xlo':
                    xlo, xhi = float(line.split()[0]), float(line.split()[1])
                    continue
                if line.split()[2] == 'ylo':
                    ylo, yhi = float(line.split()[0]), float(line.split()[1])
                    continue
                if line.split()[2] == 'zlo':
                    zlo, zhi = float(line.split()[0]), float(line.split()[1])
                    continue

                if len(line.split()) < 4: continue

                if line.split()[3] == 'xy':
                    xy, xz, yz = float(line.split()[0]), float(line.split()[1]), float(line.split()[2])
                    continue
                    
        return pd.Series([(xlo,ylo,zlo),(xhi-xlo,0.,0.),(xy,yhi-ylo,0.),(xz,yz,zhi-zlo)],
                          index=['origin','a','b','c'])

    def _count_configs(self):
        return 1

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order, 
    e.g. ['1','2','100'] instead of ['1','100','2']
    http://nedbatchelder.com/blog/200712/human_sorting.html
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
def atoi(text):
    return int(text) if text.isdigit() else text

class LAMMPS_Output(DataInput):
    """
    Data divided into two levels; sytem and atom
    
    System level data created with `fix print`, e.g.;

        fix sys_info all print 100 "${t} ${natoms} ${temp}" &    
        title "time natoms temp" file system.dump screen no
        
    Atom level data created with `dump`, e.g.;
    
        dump atom_info all custom 100 atom.dump id type x y z mass q
        OR (file per configuration)
        dump atom_info all custom 100 atom_*.dump id type xs ys zs mass q

    """
    #TODO option to ensure timesteps of atom and sys are the same
    def setup_data(self, atom_path='', sys_path='', 
                   unscale_coords=True, sys_sep=' ',
                   incl_atom_step=False,incl_sys_data=True):
        """
        Data divided into two levels; meta and atom
        
        Properties
        ----------
        unscale_coords : bool
            By default, atom coords are written in a scaled format (from 0 to 1), 
            i.e. an x value of 0.25 means the atom is at a location 1/4 of the 
            box boundaries 'a' vector. 
            http://lammps.sandia.gov/doc/dump.html?highlight=dump
        sys_sep : str
            the separator between variables in the system data file
        incl_atom_time : bool
            include time according to atom file in column 'atom_step' of meta
        incl_sys_data : bool
            include system data in the single step meta data
        
        Notes
        -----
        
        Meta level data created with `fix print`, e.g.;
    
            fix sys_info all print 100 "${t} ${natoms} ${temp}" &    
            title "time natoms temp" file system.dump screen no
            
        Atom level data created with `dump`, e.g.;
        
            dump atom_info all custom 100 atom.dump id type x y z mass q
            OR (file per configuration)
            dump atom_info all custom 100 atom_*.dump id type xs ys zs mass q

        """
        if sys_path:
            assert os.path.exists(sys_path), 'sys_path does not exist'
        self._sys_path = sys_path
        self._sys_sep = sys_sep
        
        if '*' in atom_path:            
            self._single_atom_file = False  
            self._atom_path = glob.glob(atom_path)
            assert len(self._atom_path)>0, 'atom_path does not exist'
            self._atom_path.sort(key=natural_keys)
            self._configs = len(self._atom_path)
        else:   
            self._configs = 0
            if atom_path:
                assert os.path.exists(atom_path), 'atom_path does not exist'
                with open(atom_path, 'r') as f:
                    for line in f:
                        if 'ITEM: TIMESTEP' in line: 
                            self._configs += 1
            self._single_atom_file = True
            self._atom_path = atom_path
        
        self._unscale = unscale_coords
        self._data_set = True
        self._incl_atom_step = incl_atom_step
        self._incl_sys_data = incl_sys_data
        
    def _count_configs(self):
            return self._configs

    #TODO include_bb
    def _get_meta_data_all(self, incl_bb=False):
        """ return pandas.DataFrame 

        incl_bb : bool
            include bounding box parameters 

        """
        if incl_bb:
            raise NotImplemented('must run with incl_bb=False')
        
        if self._sys_path:
            sys_df = pd.read_csv(self._sys_path, sep=self._sys_sep)
            # no data output for initial configuration
            sys_df.index += 2
        else:
            sys_df = pd.DataFrame(index=[i+1 for i in range(self.count_configs())])  

        if self._incl_atom_step and self._atom_path:
            sys_df.loc[1] = [np.nan for _ in sys_df.columns]
            sys_df.sort_index(inplace=True)
            atimes = [self._get_atom_timestep(i+1) for i in range(self.count_configs())]        
            sys_df['atom_time'] = atimes 
        
        sys_df.index.name = 'config'
            
        return sys_df

    def _get_meta_data(self, step):
        """ pandas.Series  """
        if self._sys_path and self._incl_sys_data:
            sys_df = pd.read_csv(self._sys_path, sep=self._sys_sep)
            # no data output for initial configuration
            sys_df.index += 2
            sys_df.loc[1] = [np.nan for _ in sys_df.columns]
            sys_df.sort_index(inplace=True)
            
            if sys_df.shape[0] < step:
                raise RuntimeError('the system data does not contain data for each step, \
                                    perhaps use the incl_sys_data=False in setup_data method')
                
            s1 = sys_df.loc[step]            
        else:
            s1 = pd.Series()
        
        if self._atom_path:            
            # ensure systems data doesn't already contain bounding box variable names
            for var in ['origin','a','b','c']:
                if var in s1.index:
                    old_var = s1.pop(var)
                    s1['sys_{0}'.format(var)] = old_var  

            origin,a,b,c = self._get_simulation_box(step)
            s2 = pd.Series([origin,a,b,c],index=['origin','a','b','c'])
        else:
            s2 = pd.Series()
        
        return pd.concat([s1,s2])
                        
    def _get_atom_data(self, step):
        """ return pandas.DataFrame         
        """
        if self._single_atom_file:
            current_step = 1
            with open(self._atom_path, 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        
                        if step==current_step:
                            atoms_df =  self._extract_atom_data(f, self._unscale)
                            self._add_colors(atoms_df)
                            self._add_radii(atoms_df)
                            return atoms_df
                        else:
                            current_step+=1
                            # find the number of atoms and skip that many lines
                            line = self._skiplines(f, 3)
                            self._skiplines(f, int(line.split()[0])+5) 
            
            if current_step>1:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
        elif not self._single_atom_file:
            if len(self._atom_path) < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step-1], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        atoms_df =  self._extract_atom_data(f, self._unscale)
                        self._add_colors(atoms_df)
                        self._add_radii(atoms_df)
                        return atoms_df
                        
            raise IOError("atom file of wrong format")
        
    
    def _extract_atom_data(self, f, unscale_coords=True):
        """ """
        xy, xz, yz = 0., 0., 0.
        
        line = self._skiplines(f, 1)
        #time = int(line.split()[0])
        line = self._skiplines(f, 2)
        num_atoms = int(line.split()[0])
        line = self._skiplines(f, 2)
        xlo_bound, xhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xy = float(line.split()[2])
        line = self._skiplines(f, 1)
        ylo_bound, yhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xz = float(line.split()[2])
        line = self._skiplines(f, 1)
        zlo_bound, zhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            yz = float(line.split()[2])
        line = self._skiplines(f, 1)
        headers = line.split()[2:]
        atoms = []
        for atom in range(num_atoms):
            line = self._skiplines(f, 1)
            atoms.append(np.array(line.split(),dtype=float))
        atoms_df = pd.DataFrame(atoms, columns=headers)
        
        #fix legacy issue
        atoms_df.rename(columns={'xs': 'x', 'ys': 'y', 'zs':'z'}, inplace=True)

        if unscale_coords:
            self._unscale_coords(atoms_df, 
                                 xlo_bound, xhi_bound, 
                                 ylo_bound, yhi_bound, 
                                 zlo_bound, zhi_bound, 
                                 xy, xz, yz)        

        return atoms_df
        
    def _unscale_coords(self, atoms_df, 
                        xlo_bound, xhi_bound, ylo_bound, yhi_bound, 
                        zlo_bound, zhi_bound, xy, xz, yz):
        """
        By default, atom coords are written in a scaled format (from 0 to 1), 
        i.e. an x value of 0.25 means the atom is at a location 1/4 of the 
        box boundaries 'a' vector. 
        http://lammps.sandia.gov/doc/dump.html?highlight=dump
        """
        xlo, xhi = xlo_bound - min(0.0,xy,xz,xy+xz), xhi_bound - max(0.0,xy,xz,xy+xz)
        ylo, yhi = ylo_bound - min(0.0,yz), yhi_bound - max(0.0,yz)
        zlo, zhi = zlo_bound, zhi_bound

        a,b,c = np.array([[xhi-xlo,0.,0.],[xy,yhi-ylo,0.],[xz,yz,zhi-zlo]])
        origin = np.array([xlo,ylo,zlo])
        
        new_coords = (np.array([atoms_df['x'].values]).transpose() * a + 
               np.array([atoms_df['y'].values]).transpose() * b + 
               np.array([atoms_df['z'].values]).transpose() * c + 
               origin) 
        atoms_df[['x','y','z']] = new_coords

    def _get_atom_timestep(self, step):
        """ return simulation step, according to atom data """
        if self._single_atom_file:
            current_step = 1
            with open(self._atom_path, 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        
                        if step==current_step:
                            return self._extract_time_data(f)
                        else:
                            current_step+=1
                            # find the number of atoms and skip that many lines
                            line = self._skiplines(f, 3)
                            self._skiplines(f, int(line.split()[0])+5) 
            
            if current_step>1:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
        else:
            if len(self._atom_path) < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step-1], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        return self._extract_time_data(f)
            raise IOError("atom file of wrong format")

    def _extract_time_data(self, f):
        """ """
        line = self._skiplines(f, 1)
        time = int(line.split()[0])
        
        return time
        
    def _get_simulation_box(self, step):
       """ return list of coordinates origin,a,b,c """
       if self._single_atom_file:
            current_step = 1
            with open(self._atom_path, 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        
                        if step==current_step:
                            return self._extract_simulation_box(f)
                        else:
                            current_step+=1
                            # find the number of atoms and skip that many lines
                            line = self._skiplines(f, 3)
                            self._skiplines(f, int(line.split()[0])+5) 
            
            if current_step>1:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
       else:
            if len(self._atom_path) < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step-1], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        return self._extract_simulation_box(f)
            raise IOError("atom file of wrong format")

    def _extract_simulation_box(self, f):
        """ """
        xy, xz, yz = 0., 0., 0.
        
        line = self._skiplines(f, 1) # to time
        line = self._skiplines(f, 2) # to number of atoms
        line = self._skiplines(f, 2) # to simulation box
        xlo_bound, xhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xy = float(line.split()[2])
        line = self._skiplines(f, 1)
        ylo_bound, yhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xz = float(line.split()[2])
        line = self._skiplines(f, 1)
        zlo_bound, zhi_bound = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            yz = float(line.split()[2])
       
        xlo, xhi = xlo_bound - min(0.0,xy,xz,xy+xz), xhi_bound - max(0.0,xy,xz,xy+xz)
        ylo, yhi = ylo_bound - min(0.0,yz), yhi_bound - max(0.0,yz)
        zlo, zhi = zlo_bound, zhi_bound

        return (xlo,ylo,zlo), (xhi-xlo,0.,0.),(xy,yhi-ylo,0.),(xz,yz,zhi-zlo)