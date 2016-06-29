# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""
import itertools
from .. import _colors

class DataInput(object):
    """
    Data is divided into two levels; sytem and atom
    
    system is simply a table containing variables (columns) for each timestep (rows)
    atom is a series of tables, one for each timestep, 
        containing variables (columns) for each atom (rows)
    systems bounding box is also supplied for each time step
    """
    def __init__(self):
        raise NotImplemented
    def get_system_data(self, step=None):
        """ return pandas.DataFrame or pandas.Series (for single step) """
        raise NotImplemented
    def get_atom_data(self, step):
        """ return pandas.DataFrame """
        raise NotImplemented
    def get_simulation_box(self, step):
        """ return list of coordinates (np.array(3)) for [a,b,c] & origin"""
        raise NotImplemented
        
    def count_timesteps(self):
        """ return int of total number of timesteps """
        raise NotImplemented

    def _add_radii(self, atom_df):
        atom_df['radius'] = 1.
        
    #TODO cycle through specific colors in each list
    def _add_colors(self, atom_df):
        """ add colors to atom_df, with different color for each atom type """

        atom_df['transparency'] = 1.
        
        col_keys = _colors.col_dict.keys()
        col_cycle = itertools.cycle(col_keys)
        for typ in atom_df['type'].unique():            
            atom_df.loc[atom_df['type']==typ,'color'] = _colors.col_dict[col_cycle.next()][0]
        
    def _skiplines(self, f, num=1):
        """ skip line(s) in an open file """
        for n in range(num):
            line = next(f)
        return line