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

import re

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]

def skiplines(f, num):
    """ skip line(s) in an open file """
    for n in range(num):
        line = next(f)
    return line

class MD_Data(object):
    """
    Data can be divided into two levels; sytem and atom

    """
    def __init__(self):
        raise NotImplemented
    def get_system_data(self, step=None):
        """ return pandas.DataFrame or pandas.Series (for single step) """
        raise NotImplemented
    def get_atom_data(self, step):
        """ return pandas.DataFrame """
        raise NotImplemented
        
class LAMMPS_Data(MD_Data):
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
                            line = skiplines(f, 3)
                            skiplines(f, int(line.split()[0])+5) 
            
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
        line = skiplines(f, 1)
        time = int(line.split()[0])
        line = skiplines(f, 2)
        num_atoms = int(line.split()[0])
        line = skiplines(f, 2)
        x0, x1 = [float(line.split()[0]), float(line.split()[1])]
        line = skiplines(f, 1)
        y0, y1 = [float(line.split()[0]), float(line.split()[1])]
        line = skiplines(f, 1)
        z0, z1 = [float(line.split()[0]), float(line.split()[1])]
        line = skiplines(f, 1)
        headers = line.split()[2:]
        atoms = []
        for atom in range(num_atoms):
            line = skiplines(f, 1)
            atoms.append(np.array(line.split(),dtype=float))
        atoms_df = pd.DataFrame(atoms, columns=headers)
        
        self._unscale_coords(atoms_df, x0, x1, y0, y1, z0, z1)        
        
        return atoms_df, time, np.array([[x0,x1], [y0,y1], [z0,z1]])
        
    def _unscale_coords(self, atoms_df, x0, x1, y0, y1, z0, z1):
        
        atoms_df['xs'] = atoms_df['xs'] * abs(x1-x0) + x0
        atoms_df['ys'] = atoms_df['ys'] * abs(y1-y0) + y0
        atoms_df['zs'] = atoms_df['zs'] * abs(z1-z0) + z0
        
    def count_timesteps(self):
        if not self._single_atom_file:
            return len(self._atom_path) - 1       
        elif self._sys_path:
            return sum(1 for line in open(self._sys_path)) - 1
        else:
            raise Exception('cannot compute from current data') 
