# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""

class DataInput(object):
    """
    Data is divided into two levels; sytem and atom
    
    system is simply a table containing variables (columns) for each timestep (rows)
    atom is a series of tables, one for each timestep, 
        containing variables (columns) for each atom (rows)
    the atom data also often contains a description of the systems bounding box
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

    def _skiplines(self, f, num):
        """ skip line(s) in an open file """
        for n in range(num):
            line = next(f)
        return line