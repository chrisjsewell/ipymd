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
        containing columns; xs, ys, zs, type
    """
    def __init__(self, atom_df):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns; xs, ys, zs, type
        """
        self._atom_df

    def calc_volume_pp(self, a, b, c):
        """ calculate volume of a parallelepiped determined by 3 vectors

        a,b,c : np.array(3)        
        
        """
        return a.dot(np.cross(b,c))

    def calc_volume(self):
        """ calculate volume of all points in atom_df """
        points = self._atom_df[['xs','ys','zs']].values
        hull = ConvexHull(points)
        return hull.volume