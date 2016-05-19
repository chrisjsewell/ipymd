# -*- coding: utf-8 -*-
"""
Created on Mon May 16 01:23:11 2016

@author: cjs14

Adapted from chemlab which, in turn, was adapted from
ASE https://wiki.fysik.dtu.dk/ase/
Copyright (C) 2010, Jesper Friis

"""
import os
import numpy as np
import pandas as pd

from chemlab.core.spacegroup import Spacegroup
from chemlab.core.spacegroup.cell import cellpar_to_cell

from .base import DataInput

def get_spacegroup_df():
    """ dataframe of spacegroup mappings """
    sg_path =  os.path.join(os.path.dirname(os.path.abspath(__file__)), 
                            'spacegroups.csv')
    return pd.read_csv(sg_path, index_col=0)

class Crystal(DataInput):
    """Build a crystal from atomic positions, space group and cell
    parameters.

    """
    def __init__(self, positions, atom_type, group,
            cellpar=[1.0, 1.0, 1.0, 90, 90, 90], repetitions=[1, 1, 1],
            mass_map={}, charge_map={}):
        """Build a crystal from atomic positions, space group and cell
        parameters (in Angstroms)
        
        Parameters
        -----------
    
        positions : list of coordinates
            A list of the fractional atomic positions 
        atom_type : list of atom type
            The atom types corresponding to the positions, the atoms will be
            translated in all the equivalent positions.
        group : int | str
            Space group given as its number in International Tables
            NB: to see mappings from Hermann–Mauguin notation, etc, 
            use the get_spacegroup_df function in this module
        repetitions :
            Repetition of the unit cell in each direction
        cellpar :
            Unit cell parameters (in nm and degrees)
        mass_map : dict of atom masses
            mapping of atom masses to atom types
        charge_map : dict of atom charges
            mapping of atom charges to atom types
    
        This function was taken and adapted from the *spacegroup* module 
        found in `ASE <https://wiki.fysik.dtu.dk/ase/>`_.
    
        The module *spacegroup* module was originally developed by Jesper
        Frills.
        
        Example
        -------
        
        from ipymd.data_input import crystal
        c = crystal.Crystal([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
            ['Na', 'Cl'], 225,
            cellpar = [.54, .54, .54, 90, 90, 90],
            repetitions = [5, 5, 5])
        c.get_atom_data()
        c.get_simulation_box()
    
        """
        sp = Spacegroup(group)
        sites, kind = sp.equivalent_sites(positions)
        
        nx, ny, nz = repetitions
        
        atoms = []
        aid = 1
        
        # Unit cell parameters
        a,b,c = cellpar_to_cell(cellpar) * 10.
        
        self._sim_box = (np.array([a*nx, b*ny, c*nz]),np.array([0.,0.,0.]))
        
        for rx in range(nx):
            for ry in range(ny):
                for rz in range(nz):
                    for s, ki in zip(sites, kind):
                        atype = atom_type[ki]
                        x,y,z = s[0]*a +s[1]*b + s[2]*c + a*rx + b*ry + c*rz
                        atoms.append([aid,atype,x,y,z])
                        aid+=1
                
        self._atoms = pd.DataFrame(atoms,columns=['id','type','xs','ys','zs'])
        
        if mass_map:
            self._atoms['mass'] = np.nan
            for typ in atom_type:
                self._atoms.loc[self._atoms['type']== typ,'mass'] = mass_map[typ]
        if charge_map:
            self._atoms['q'] = np.nan
            for typ in atom_type:
                self._atoms.loc[self._atoms['type']== typ,'q'] = charge_map[typ]

        self._add_colors(self._atoms)
        self._add_radii(self._atoms)
  
    def get_atom_data(self):
        """ return atom data """
        return self._atoms.copy()

    def get_simulation_box(self):
        """ return list of coordinates origin & [a,b,c] """
        return self._sim_box

