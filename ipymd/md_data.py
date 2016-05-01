# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""
import numpy as np
import pandas as pd

class MD_Data(object):
    """
    Data can be divided into two levels; sytem and atom

    """
    def __init__(self):
        raise NotImplemented
    def get_system_data_all(self):
        """ return pandas.DataFrame """
        raise NotImplemented
    def get_system_data(self, step):
        """ return pandas.Series """
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

    For style *real*, these are the units:
    
        mass = grams/mole
        distance = Angstroms
        time = femtoseconds
        energy = Kcal/mole
        velocity = Angstroms/femtosecond
        force = Kcal/mole-Angstrom
        torque = Kcal/mole
        temperature = Kelvin
        pressure = atmospheres
        dynamic viscosity = Poise
        charge = multiple of electron charge (1.0 is a proton)
        dipole = charge*Angstroms
        electric field = volts/Angstrom
        density = gram/cm^dim

    """
    def __init__(self, sys_path='', atom_path=''):
        self.sys_path = sys_path
        self.atom_path = atom_path
        
    #TODO include simulation bounds
    def get_system_data_all(self):
        """ return pandas.DataFrame """
        sys_df = pd.read_csv(self.sys_path, sep=' ',)
        sys_df.index += 1
        return sys_df
        
    def get_system_data(self, step):
        """ return pandas.Series """
        sys_df = pd.read_csv(self.sys_path, sep=' ',)
        sys_df.index += 1
        return sys_df.loc[step]
        
    def _skiplines(self, f, num):
        """ skip line(s) in a file """
        for n in range(num):
            line = next(f)
        return line
        
    def _convert_units(self, atom_df):
        """ convert x,y,z from Angstrom to nm """
        for head in ['xs','ys','zs']:
            atom_df[head] = atom_df[head]*10
        return atom_df
        
    def get_atom_data(self, step):
        """ return pandas.DataFrame """
        found_timestep = False
        with open(self.atom_path, 'r') as f:
            for line in f:
                if 'ITEM: TIMESTEP' in line: 
                    line = self._skiplines(f, 1)
                    if step==int(line.split()[0]):
                        found_timestep = True
                        line = self._skiplines(f, 2)
                        num_atoms = int(line.split()[0])
                        line = self._skiplines(f, 5)
                        headers = line.split()[2:]
                        atoms = []
                        for atom in range(num_atoms):
                            line = self._skiplines(f, 1)
                            atoms.append(np.array(line.split(),dtype=float))
                        atoms_df = pd.DataFrame(atoms, columns=headers)
    
                        break
                    else:
                        # find the number of atoms and skip that many lines
                        line = self._skiplines(f, 2)
                        self._skiplines(f, int(line.split()[0])+5) 
        
        if found_timestep:
            return self._convert_units(atoms_df)
        else:
            raise IOError("couldn't find required timestep")

    def count_steps(self):
        return sum(1 for line in open(self.sys_path)) - 1