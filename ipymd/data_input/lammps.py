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
    def __init__(self, atom_path=''):
        assert os.path.exists(atom_path) or not atom_path, 'atom_path does not exist'
        self._atom_path = atom_path
    def get_atom_data(self, atom_style='atomic'):        
        """ get data from file
        
        atom_style : 'atomic', 'charge'
            defines how atomic data is listed:
            atomic; atom-ID atom-type x y z
            charge; atom-ID atom-type q x y z
        
        """
        #TODO add more atom styles
        assert atom_style in ['atomic', 'charge']
        
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
                        if atom_style == 'atomic':
                            aid, atype, x, y, z = (int(line.split()[0]), 
                                                   int(line.split()[1]), 
                                                   float(line.split()[2]),
                                                   float(line.split()[3]),
                                                   float(line.split()[4]))
                            atom_data.append([aid, atype, x, y, z])
                        elif atom_style == 'charge':
                            aid, atype, q, x, y, z = (int(line.split()[0]), 
                                                   int(line.split()[1]), 
                                                   float(line.split()[2]),
                                                   float(line.split()[3]),
                                                   float(line.split()[4]),
                                                   float(line.split()[5]))
                            atom_data.append([aid, atype, q, x, y, z])
                        line = self._skiplines(f,1)
                        continue
        
                if len(line.split()) < 2: continue

                if line.split()[1] == 'atoms':
                    num_atoms = int(line.split()[0])
                    continue

        if atom_style == 'atomic':
            atom_df = pd.DataFrame(atom_data,columns=['id', 'type', 'xs', 'ys', 'zs'])
        elif atom_style == 'charge':
            atom_df = pd.DataFrame(atom_data,columns=['id', 'type', 'q', 'xs', 'ys', 'zs'])
            
        return atom_df
    
    def get_simulation_box(self):
        """ return list of coordinates (np.array(3)) for origin and [a,b,c] """
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
                    
        return np.array([[xhi-xlo,0.,0.],[xy,yhi-ylo,0.],[xz,yz,zhi-zlo]]), np.array([xlo,ylo,zlo])

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
    
        dump atom_info all custom 100 atom.dump id type xs ys zs mass q
        OR (file per configuration)
        dump atom_info all custom 100 atom_*.dump id type xs ys zs mass q

    """
    def __init__(self, sys_path='', atom_path='', unscale_coords=True):
        """
        unscale_coords : bool
            By default, atom coords are written in a scaled format (from 0 to 1), 
            i.e. an x value of 0.25 means the atom is at a location 1/4 of the 
            distance from xlo to xhi of the box boundaries. 
            http://lammps.sandia.gov/doc/dump.html?highlight=dump
        """
        assert os.path.exists(sys_path) or not sys_path, 'sys_path does not exist'
        self._sys_path = sys_path
        
        if '*' in atom_path:            
            self._single_atom_file = False  
            self._atom_path = glob.glob(atom_path)
            assert len(self._atom_path)>0, 'atom_path does not exist'
            self._atom_path.sort(key=natural_keys)
        else:        
            assert os.path.exists(atom_path) or not atom_path, 'atom_path does not exist'
            self._single_atom_file = True
            self._atom_path = atom_path
        
    def get_system_data(self, step=None):
        """ return pandas.DataFrame or pandas.Series (for single step) """
        sys_df = pd.read_csv(self._sys_path, sep=' ',)
        sys_df.index += 1
        if step is None:
            return sys_df
        else:
            return sys_df.loc[step]
                        
    def get_atom_data(self, step):
        """ return pandas.DataFrame """
        if self._single_atom_file:
            current_step = 0
            with open(self._atom_path, 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        
                        if step==current_step:
                            return self._extract_atom_data(f)
                        else:
                            current_step+=1
                            # find the number of atoms and skip that many lines
                            line = self._skiplines(f, 3)
                            self._skiplines(f, int(line.split()[0])+5) 
            
            if current_step>0:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
        else:
            if len(self._atom_path)-1 < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        return self._extract_atom_data(f)
            raise IOError("atom file of wrong format")
    
    def _extract_atom_data(self, f):
        """ """
        line = self._skiplines(f, 1)
        #time = int(line.split()[0])
        line = self._skiplines(f, 2)
        num_atoms = int(line.split()[0])
        line = self._skiplines(f, 2)
        xlo, xhi = [float(line.split()[0]), float(line.split()[1])]
        line = self._skiplines(f, 1)
        ylo, yhi = [float(line.split()[0]), float(line.split()[1])]
        line = self._skiplines(f, 1)
        zlo, zhi = [float(line.split()[0]), float(line.split()[1])]
        line = self._skiplines(f, 1)
        headers = line.split()[2:]
        atoms = []
        for atom in range(num_atoms):
            line = self._skiplines(f, 1)
            atoms.append(np.array(line.split(),dtype=float))
        atoms_df = pd.DataFrame(atoms, columns=headers)
        
        self._unscale_coords(atoms_df, xlo, xhi, ylo, yhi, zlo, zhi)        
        
        return atoms_df
        
    def _unscale_coords(self, atoms_df, x0, x1, y0, y1, z0, z1):
        
        atoms_df['xs'] = atoms_df['xs'] * abs(x1-x0) + x0
        atoms_df['ys'] = atoms_df['ys'] * abs(y1-y0) + y0
        atoms_df['zs'] = atoms_df['zs'] * abs(z1-z0) + z0

    def get_atom_timestep(self, step):
        """ return simulation step, according to atom data """
        if self._single_atom_file:
            current_step = 0
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
            
            if current_step>0:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
        else:
            if len(self._atom_path)-1 < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        return self._extract_time_data(f)
            raise IOError("atom file of wrong format")

    def _extract_time_data(self, f):
        """ """
        line = self._skiplines(f, 1)
        time = int(line.split()[0])
        
        return time
        
    def get_simulation_box(self, step):
       """ return list of coordinates (np.array(3)) for origin and [a,b,c] """
       if self._single_atom_file:
            current_step = 0
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
            
            if current_step>0:
                raise IOError("timestep {0} exceeds maximum ({1})".format(
                                                        step, current_step-1))
            else:
                raise IOError("atom file of wrong format")
       else:
            if len(self._atom_path)-1 < step:
                raise IOError("timestep {0} exceeds maximum ({1})".format
                                                (step, len(self._atom_path)-1))
            with open(self._atom_path[step], 'r') as f:
                for line in f:
                    if 'ITEM: TIMESTEP' in line: 
                        return self._extract_simulation_box(f)
            raise IOError("atom file of wrong format")

    def _extract_simulation_box(self, f):
        """ """
        xy, xz, yz = 0., 0., 0.
        
        line = self._skiplines(f, 1) # to time
        line = self._skiplines(f, 2) # to nummper of atoms
        line = self._skiplines(f, 2) # to simulation box
        xlo, xhi = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xy = float(line.split()[2])
        line = self._skiplines(f, 1)
        ylo, yhi = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            xz = float(line.split()[2])
        line = self._skiplines(f, 1)
        zlo, zhi = [float(line.split()[0]), float(line.split()[1])]
        if len(line.split()) == 3:
            yz = float(line.split()[2])
        
        return np.array([[xhi-xlo,0.,0.],[xy,yhi-ylo,0.],[xz,yz,zhi-zlo]]), np.array([xlo,ylo,zlo])

    def count_timesteps(self):
        if not self._single_atom_file:
            return len(self._atom_path) - 1       
        elif self._sys_path:
            return sum(1 for line in open(self._sys_path)) - 1
        else:
            raise Exception('cannot compute from current data')         