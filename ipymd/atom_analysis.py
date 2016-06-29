# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:19:07 2016

@author: cjs14
"""

#import math
#import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull, cKDTree

from .atom_manipulation import Atom_Manipulation

class Atom_Analysis(object):
    """ a class to analyse atom data
    
    atom_df : pandas.DataFrame
        containing columns; x, y, z, type, mass
    """
    def __init__(self):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns based on analysis, e.g.; x, y, z, type, mass,...
        """
        pass
    
    def calc_volume_bb(self, vectors=np.array([[1,0,0],[0,1,0],[0,0,1]])):
        """ calculate volume of the bounding box        
        """
        a,b,c = vectors
        return a.dot(np.cross(b,c))

    def calc_density_bb(self,atoms_df, vectors=np.array([[1,0,0],[0,1,0],[0,0,1]])):
        """ calculate density of the bounding box (assuming all atoms are inside)
        """
        assert set(atoms_df.columns).issuperset(['mass'])
        mass = atoms_df['mass'].sum()
        vol = self.calc_volume_bb(vectors)
        
        return mass/vol
        
    def calc_volume_points(self,atoms_df):
        """ calculate volume of the shape encompasing all atom coordinates """
        assert set(atoms_df.columns).issuperset(['x','y','z'])
        points = atoms_df[['x','y','z']].values
        hull = ConvexHull(points)
        return hull.volume
    
    def calc_coordination(self,coord_atoms_df, lattice_atoms_df, max_dist=3, max_coord=10,
                          min_dist=0.01, leafsize=100):
        """ calculate the coordination number of each atom in coords_atoms, w.r.t lattice_atoms
        
        coords_atoms_df : pandas.Dataframe
            atoms to calcualte coordination of
        lattice_atoms_df : pandas.Dataframe
            atoms to act as lattice for coordination
        max_dist : float
            maximum distance for coordination consideration
        max_coord : float
            maximum possible coordination number
        min_dist : float
            lattice points within this distance of the atom will be ignored (assumed self-interaction)
        leafsize : int
            points at which the algorithm switches to brute-force (kdtree specific)
        
        Returns
        -------
        coords : list
            list of coordination numbers
        
        """
        lattice_tree = cKDTree(lattice_atoms_df[['x','y','z']].values, leafsize=leafsize)
        coords = []
        for atom in coord_atoms_df[['x','y','z']].values:
            dists,ids = lattice_tree.query(atom, k=max_coord, distance_upper_bound=max_dist)
            coords.append(np.count_nonzero(np.logical_and(dists>min_dist, dists<np.inf)))
        return coords
    
    def calc_type_coordination(self, atoms_df, coord_type, lattice_type, max_dist=3, max_coord=10,
                          repeat_vectors=None, min_dist=0.01, leafsize=100):
        """ returns dataframe with additional column for the coordination number of 
        each atom of coord type, w.r.t lattice_type atoms
        
        effectively an extension of calc_df_coordination
        
        atoms_df : pandas.Dataframe
            all atoms
        coord_type : string
            atoms to calcualte coordination of
        lattice_atoms_df : string
            atoms to act as lattice for coordination
        max_dist : float
            maximum distance for coordination consideration
        max_coord : float
            maximum possible coordination number
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors
        min_dist : float
            lattice points within this distance of the atom will be ignored (assumed self-interaction)
        leafsize : int
            points at which the algorithm switches to brute-force (kdtree specific)
        
        Returns
        -------
        df : pandas.Dataframe
            copy of atoms_df with new column named coord_{coord_type}_{lattice_type}
        
        """
        df = atoms_df.copy()
        df['coord_{0}_{1}'.format(coord_type, lattice_type)] = np.nan      
        
        coord_df = Atom_Manipulation(df)
        coord_df.filter_variables(coord_type)

        
        lattice_df = Atom_Manipulation(df)
        lattice_df.filter_variables(lattice_type)
        
        if repeat_vectors is not None:
            lattice_df.repeat_cell(repeat_vectors,((-1,1),(-1,1),(-1,1)))
        
        coords = self.calc_coordination(coord_df.df,lattice_df.df,max_dist, max_coord,
                                        min_dist, leafsize)
                                        
    
        df.loc[df['type']==coord_type,'coord_{0}_{1}'.format(coord_type, lattice_type)] = coords
        
        return df
        
    #TODO vacancy identification
    def compare_to_lattice(self, atoms_df, lattice_atoms_df, max_dist=10,leafsize=100):
        """ calculate the minimum distance of each atom in atoms_df from a lattice point in lattice_atoms_df
        
        atoms_df : pandas.Dataframe
            atoms to calculate for
        lattice_atoms_df : pandas.Dataframe
            atoms to act as lattice points
        max_dist : float
            maximum distance for consideration in computation
        leafsize : int
            points at which the algorithm switches to brute-force (kdtree specific)
        
        Returns
        -------
        distances : list
            list of distances to nearest atom in lattice
        
        """
        lattice_tree = cKDTree(lattice_atoms_df[['x','y','z']].values, leafsize=leafsize)
        min_dists = []
        for atom in atoms_df[['x','y','z']].values:
            dist,idnum = lattice_tree.query(atom, k=1, distance_upper_bound=max_dist)
            min_dists.append(dist)
        return min_dists
            

        
                              