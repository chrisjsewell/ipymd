# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""
import itertools
from ..shared import colors

class DataInput(object):
    """data input base class
    
    subclasses should override methods; 
    read_data, _get_atom_data, _get_meta_data and _count_configs
    and optionally _get_meta_data_all
    
    """
    def __init__(self):
        """
        Data is divided into two levels; atomic and meta
        
            - atom is a series of tables, one for each timestep, containing variables (columns) for each atom (rows)
            - meta is a table containing variables (columns) for each configuration (rows)
        """
        self._data_set = False

    def setup_data(self):
        """a method to setup the data and variables """
        self._data_set = True

    def get_atom_data(self, config=1):
        """ return pandas.DataFrame of atomic data """
        if not self._data_set:
            raise RuntimeError('must call setup_data method first')
        if config>self.count_configs():
            raise ValueError('only {} configurations available'.format(
                                                            self.count_configs()))
        return self._get_atom_data(config)

    def _get_atom_data(self, config):
        raise NotImplemented        

    def get_meta_data(self, config=1):
        """ return pandas.Series of meta data for the atomic configuration """
        if not self._data_set:
            raise RuntimeError('must call setup_data method first')
        if config>self.count_configs():
            raise ValueError('only {} configurations available'.format(
                                                            self.count_configs()))
        return self._get_meta_data(config)
            
    def _get_meta_data(self, config):
        raise NotImplemented

    def get_meta_data_all(self, incl_bb=False, **kwargs):
        """ return pandas.DataFrame of meta data for the atomic configuration

        Properties
        ----------
        incl_bb : bool
            whether to include bounding box coordinates in DataFrame
        kwargs : dict
            kew word arguments relevant to specific input data
        
        """
        if not self._data_set:
            raise RuntimeError('must call setup_data method first')
        return self._get_meta_data_all(incl_bb, **kwargs)
    
    def _get_meta_data_all(self, incl_bb):
        raise NotImplemented
    
    def count_configs(self):
        """ return int of total number of atomic configurations """
        if not self._data_set:
            raise RuntimeError('must call setup_data method first')
        return self._count_configs()
    
    def _count_configs(self):
        raise NotImplemented

    def _add_radii(self, atom_df):
        atom_df['radius'] = 1.
        
    #TODO cycle through specific colors in each list
    def _add_colors(self, atom_df):
        """ add colors to atom_df, with different color for each atom type """

        atom_df['transparency'] = 1.
        
        col_keys = colors.col_dict.keys()
        col_cycle = itertools.cycle(col_keys)
        for typ in atom_df['type'].unique():            
            atom_df.loc[atom_df['type']==typ,'color'] = colors.col_dict[col_cycle.next()][0]
        
    def _skiplines(self, f, num=1):
        """ skip line(s) in an open file """
        for n in range(num):
            line = next(f)
        return line