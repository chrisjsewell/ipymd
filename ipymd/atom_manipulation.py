# -*- coding: utf-8 -*-
"""
Created on Mon May 16 08:15:13 2016

@author: cjs14
"""
import math
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull

class Atom_Manipulation(object):
    """ a class to manipulate atom data
    
    atom_df : pandas.DataFrame
        containing columns; xs, ys, zs, type
    """
    def __init__(self):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns; xs, ys, zs, type
        """
        pass

    def change_variable(self, atom_df, old, new, vtype='type'):
        atom_copy = atom_df.copy()
        atom_copy.loc[atom_copy[vtype]==old, vtype] = new
        return atom_copy

    def filter_variable(self, atom_df, values, vtype='type'):
        if isinstance(values, basestring):
            values = [values]
        return atom_df[atom_df[vtype].isin(values)].copy()

    def _pnts_in_pointcloud(self, points, new_pts):
        """2D or 3D
        
        returns np.array(dtype=bool)
        """
        hull = ConvexHull(points)
        vol = hull.volume
        
        inside = []
        for pt in new_pts:
            new_hull = ConvexHull(np.append(points, [pt],axis=0))
            inside.append(vol == new_hull.volume)
        return np.array(inside)
    
    def filter_inside_pts(self, atom_df, points):
        """return only atoms inside points 

        points : np.array((N,3))        
        """ 
        inside = self._pnts_in_pointcloud(points, atom_df[['xs','ys','zs']].values)
        return atom_df[inside].copy()

    def filter_inside_box(self, atom_df, vectors, origin=np.zeros(3)):
        """return only atoms inside box
        
        vectors : np.array((3,3))
            a, b, c vectors
        origin : np.array((1,3))
        
        """
        a,b,c = vectors + origin
        points = [origin, a, b, a+b, c, a+c, b+c, a+b+c]
        return self.filter_inside_pts(atom_df, points)
    
    def _rotate(self, v, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta degrees.
        """
        axis = np.asarray(axis)
        theta = np.asarray(theta)*np.pi/180.
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2.0)
        b, c, d = -axis*math.sin(theta/2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        rotation_matrix = np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])    
        return np.dot(rotation_matrix, v)

    def filter_inside_hexagon(self, atom_df, vectors, origin=np.zeros(3)):
        """return only atoms inside hexagonal prism
        
        vectors : np.array((2,3))
            a, c vectors
        origin : np.array((1,3))
        
        """
        a, c = vectors
        points = [self._rotate(a, c, angle) for angle in [0,60,120,180,240,300]]
        points += [p + c for p in points]
        points = np.array(points) + origin
        return self.filter_inside_pts(atom_df, points)

    def repeat_cell(self, atom_df, vectors, repetitions=((0,1),(0,1),(0,1))):
        """ repeat atoms along vectors a, b, c  """
        dfs = []        
        for i in range(repetitions[0][0], repetitions[0][1]+1):
            for j in range(repetitions[1][0], repetitions[1][1]+1):
                for k in range(repetitions[2][0], repetitions[2][1]+1):
                    atom_copy = atom_df.copy()
                    atom_copy[['xs','ys','zs']] = (atom_copy[['xs','ys','zs']]
                                + i*vectors[0]  + j*vectors[1] + k*vectors[2])
                    dfs.append(atom_copy)
        return pd.concat(dfs)
        
    def slice_x(self, atom_df, minval=None, maxval=None):
        atom_copy = atom_df.copy()
        if minval is not None:
            atom_copy = atom_copy[atom_copy['xs']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_copy[atom_copy['xs']<=maxval].copy()
            
        return atom_copy

    def slice_y(self, atom_df, minval=None, maxval=None):
        atom_copy = atom_df.copy()
        if minval is not None:
            atom_copy = atom_copy[atom_copy['ys']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_copy[atom_copy['ys']<=maxval].copy()
            
        return atom_copy

    def slice_z(self, atom_df, minval=None, maxval=None):
        atom_copy = atom_df.copy()
        if minval is not None:
            atom_copy = atom_copy[atom_copy['zs']>=minval].copy()
        if maxval is not None:
            atom_copy = atom_copy[atom_copy['zs']<=maxval].copy()
            
        return atom_copy
                        
    #TODO slice along arbitrary direction
        