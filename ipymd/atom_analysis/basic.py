# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 14:06:40 2016

@author: cjs14

functions to calculate basic properties of the atoms

"""
import numpy as np
from scipy.spatial import ConvexHull

from ..shared.transformations import angle_between_vectors

def volume_bb(vectors=[[1,0,0],[0,1,0],[0,0,1]], rounded=None,
              cells=(1,1,1)):
    """ calculate volume of the bounding box        

    Parmeters
    ---------
    rounded : int
        the number of decimal places to return
    cells : (int,int,int)
        how many unit cells the vectors represent in each direction
    
    Returns
    -------
    volume : float
    
    """
    a,b,c = np.asarray(vectors)
    a,b,c = a/cells[0],b/cells[1],c/cells[2]
    vol = a.dot(np.cross(b,c))
    if rounded is not None:
        vol = round(vol, rounded)
    return vol

def lattparams_bb(vectors=[[1,0,0],[0,1,0],[0,0,1]],
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
    a, b, c, alpha, beta, gamma (in degrees)
    
    """
    
    a,b,c = np.asarray(vectors)
    a0,b0,c0 = cells
    a = a/a0
    b = b/b0
    c = c/c0
    
    a_length = np.linalg.norm(a)
    b_length = np.linalg.norm(b)
    c_length = np.linalg.norm(c)
    
    alpha = np.degrees(angle_between_vectors(b,c))
    beta = np.degrees(angle_between_vectors(a,c))
    gamma = np.degrees(angle_between_vectors(a,b))

    if rounded is not None:
        a_length = round(a_length, rounded)
        b_length = round(b_length, rounded)
        c_length = round(c_length, rounded)
        alpha = round(alpha, rounded)
        beta = round(beta, rounded)
        gamma = round(gamma, rounded)
    
    return a_length, b_length, c_length, alpha, beta, gamma
    
    
def density_bb(atoms_df, vectors=[[1,0,0],[0,1,0],[0,0,1]]):
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
