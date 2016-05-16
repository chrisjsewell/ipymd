# -*- coding: utf-8 -*-
"""
Created on Mon May 16 08:15:13 2016

@author: cjs14
"""
import pandas as pd

class Atom_Manipulation(object):
    """ a class to manipulate atom data
    """
    def __init__(self):
        pass

    def change_variable(self, atom_df, old, new, vtype='type'):
        atom_copy = atom_df.copy()
        atom_copy.loc[atom_copy[vtype]==old, vtype] = new
        return atom_copy

    def filter_variable(self, atom_df, values, vtype='type'):
        if isinstance(values, basestring):
            values = [values]
        return atom_df[atom_df[vtype].isin(values)].copy()

    def repeat_cell(self, atom_df, vectors, repetitions=(1,1,1)):
        """ repeat atoms along vectors a, b, c  """
        dfs = []        
        for i in range(repetitions[0]+1):
            for j in range(repetitions[1]+1):
                for k in range(repetitions[2]+1):
                    atom_copy = atom_df.copy()
                    atom_copy[['xs','ys','zs']] = (atom_copy[['xs','ys','zs']]
                                + i*vectors[0]  + j*vectors[1] + k*vectors[2])
                    dfs.append(atom_copy)
        return pd.concat(dfs)
        
    def slice_x(self, atom_df, minval=None, maxval=None):
        if minval is not None:
            atom_copy = atom_df[atom_df['xs']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_df[atom_df['xs']<=maxval].copy()
        if minval is None and maxval is None:
            atom_copy = atom_df.copy()
        return atom_copy

    def slice_y(self, atom_df, minval=None, maxval=None):
        if minval is not None:
            atom_copy = atom_df[atom_df['ys']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_df[atom_df['ys']<=maxval].copy()
        if minval is None and maxval is None:
            atom_copy = atom_df.copy()
        return atom_copy

    def slice_z(self, atom_df, minval=None, maxval=None):
        if minval is not None:
            atom_copy = atom_df[atom_df['zs']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_df[atom_df['zs']<=maxval].copy()
        if minval is None and maxval is None:
            atom_copy = atom_df.copy()
        return atom_copy
                        
            
        