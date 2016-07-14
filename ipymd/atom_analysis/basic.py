# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:06:40 2016

@author: cjs14

functions to calculate basic properties of the atoms

"""
import numpy as np
from scipy.spatial import ConvexHull

def _unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def _angle_between(v1, v2, rounded=None):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = _unit_vector(v1)
    v2_u = _unit_vector(v2)
    angle =  np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    if rounded is not None:
        angle = round(angle,rounded)
    return angle

def volume_bb(vectors=np.array([[1,0,0],[0,1,0],[0,0,1]]), rounded=None):
    """ calculate volume of the bounding box        
    """
    a,b,c = vectors
    vol = a.dot(np.cross(b,c))
    if rounded is not None:
        vol = round(vol, rounded)
    return vol

def lattparams_bb(vectors=np.array([[1,0,0],[0,1,0],[0,0,1]]),
                     rounded=None, cells=(1,1,1)):
    """ calculate unit cell parameters of the bounding box 
    
    Parmeters
    ---------
    rounded : int
        the number of decimal places to return
    cells : (int,int,int)
        how many unit cells the vectors represent in each direction
    
    Returns
    -------
    a, b, c, alpha, beta, gamma
    
    """
    
    a,b,c = vectors
    a0,b0,c0 = cells
    a = a/a0
    b = b/b0
    c = c/c0
    
    a_length = np.linalg.norm(a)
    b_length = np.linalg.norm(b)
    c_length = np.linalg.norm(c)
    
    if rounded is not None:
        a_length = round(a_length, rounded)
        b_length = round(b_length, rounded)
        c_length = round(c_length, rounded)

    alpha = _angle_between(b,c, rounded)
    beta = _angle_between(a,c, rounded)
    gamma = _angle_between(a,b, rounded)
    
    return a_length, b_length, c_length, alpha, beta, gamma
    
    
def density_bb(atoms_df, vectors=np.array([[1,0,0],[0,1,0],[0,0,1]])):
    """ calculate density of the bounding box (assuming all atoms are inside)
    """
    assert set(atoms_df.columns).issuperset(['mass'])
    mass = atoms_df['mass'].sum()
    vol = volume_bb(vectors)
    
    return mass/vol
    
def volume_points(atoms_df):
    """ calculate volume of the shape encompasing all atom coordinates """
    assert set(atoms_df.columns).issuperset(['x','y','z'])
    points = atoms_df[['x','y','z']].values
    hull = ConvexHull(points)
    return hull.volume
