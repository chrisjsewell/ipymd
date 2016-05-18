# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:19:07 2016

@author: cjs14
"""

#import math
#import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

class Atom_Analysis(object):
    """ a class to analyse atom data
    
    atom_df : pandas.DataFrame
        containing columns; xs, ys, zs, type, mass
    """
    def __init__(self, atom_df, bounding_vectors=np.array([[1,0,0],[0,1,0],[0,0,1]])):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns based on analysis, e.g.; xs, ys, zs, type, mass
        """
        assert set(atom_df.columns).issuperset(['xs','ys','zs','type'])
         
        self._atom_df = atom_df.copy()
        self._vectors = bounding_vectors

    def calc_volume(self):
        """ calculate volume of the bounding box        
        """
        a,b,c = self._vectors
        return a.dot(np.cross(b,c))

    def calc_density(self):
        """ calculate density of the bounding box (assuming all atoms are inside)
        """
        mass = self._atom_df['mass'].sum()
        vol = self.calc_volume()
        
        return mass/vol
        

    def calc_volume_points(self):
        """ calculate volume of the shape encompasing all atom coordinates """
        points = self._atom_df[['xs','ys','zs']].values
        hull = ConvexHull(points)
        return hull.volume
    
        