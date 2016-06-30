# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:19:07 2016

@author: cjs14
"""

#import math
#import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull, cKDTree
from collections import Counter
import pandas as pd

from .atom_manipulation import Atom_Manipulation

def _chain(start, nodes,nondirectional=True):
    """a recursive function to compute the maximum unbroken chain given a set of nodes
    
    e.g. start=[0,1], nodes=[[1,2],[0,1],[2,3],[8,9]]
        -> [[0,1],[1,2],[2,3]]
        
    nondirectional : bool
        whether nodes are directional, e,g, if [0,2],[1,2] can form a chain
    
    """
    remaining = list(nodes)
    if nondirectional and start not in remaining:
        del remaining[remaining.index(sorted(start))]
    else:
        del remaining[remaining.index(start)]
    possibles = [x for x in remaining if start[1] == x[0]] 
    if nondirectional:
        possibles += [[x[1],x[0]] for x in remaining if start[1] == x[1]]
    maxchain = []
    for c in possibles:
        l = _chain(c, remaining,nondirectional)
        if len(l) > len(maxchain):
            maxchain = l
    return [start] + maxchain
    

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
    
    def calc_coordination(self,coord_atoms_df, lattice_atoms_df, max_dist=4, max_coord=16,
                          repeat_vectors=None, min_dist=0.01, leafsize=100):
        """ calculate the coordination number of each atom in coords_atoms, w.r.t lattice_atoms
        
        coords_atoms_df : pandas.Dataframe
            atoms to calcualte coordination of
        lattice_atoms_df : pandas.Dataframe
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
        coords : list
            list of coordination numbers
        
        """
        lattice_df = Atom_Manipulation(lattice_atoms_df)
        
        if repeat_vectors is not None:
            lattice_df.repeat_cell(repeat_vectors,((-1,1),(-1,1),(-1,1)))

        lattice_tree = cKDTree(lattice_df.df[['x','y','z']].values, leafsize=leafsize)
        coords = []
        for atom in coord_atoms_df[['x','y','z']].values:
            dists,ids = lattice_tree.query(atom, k=max_coord, distance_upper_bound=max_dist)
            coords.append(np.count_nonzero(np.logical_and(dists>min_dist, dists<np.inf)))
        return coords
    
    def calc_type_coordination(self, atoms_df, coord_type, lattice_type, max_dist=4, max_coord=16,
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
                
        coords = self.calc_coordination(coord_df.df,lattice_df.df,max_dist, max_coord,
                                        repeat_vectors, min_dist, leafsize)
                                        
    
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
            
    #http://www.ovito.org/manual/particles.modifiers.common_neighbor_analysis.html
    #https://www.quora.com/Given-a-set-of-atomic-types-and-coordinates-from-an-MD-simulation-is-there-a-good-algorithm-for-determining-its-likely-crystal-structure?__filter__=all&__nsrc__=2&__snid3__=179254150
    #TODO add adaptive cna http://iopscience.iop.org/article/10.1088/0965-0393/20/4/045021/pdf
    def common_neighbour_analysis(self, atoms_df, upper_bound=4, max_neighbours=24,
                                  repeat_vectors=None, leafsize=100):
        """ compute atomic environment of each atom in atoms_df
        
        Based on Faken, Daniel and Jónsson, Hannes,
        'Systematic analysis of local atomic structure combined with 3D computer graphics',
        March 1994, DOI: 10.1016/0927-0256(94)90109-0
        
        ideally:
        - FCC = 12 x 4,2,1
        - HCP = 6 x 4,2,1 & 6 x 4,2,2
        - BCC = 6 x 6,6,6 & 8 x 4,4,4
        - icosahedral = 12 x 5,5,5
        
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors

        Returns
        -------
        df : pandas.Dataframe
            copy of atoms_df with new column named cna
        """
        df = atoms_df.copy()
        max_id = df.shape[0] - 1 # starts at 0
        
        if repeat_vectors is not None:
            repeat = Atom_Manipulation(df)
            repeat.repeat_cell(repeat_vectors,((-1,1),(-1,1),(-1,1)),original_first=True)
            lattice_df = repeat.df
        else:
            lattice_df = df
        
        # create nearest neighbours dictionary
        lattice_tree = cKDTree(lattice_df[['x','y','z']].values, leafsize=leafsize)
        nn_ids = {}
        #nn_dists = {}
        for atom_xyz in lattice_df[['x','y','z']].values:
            dists,ids = lattice_tree.query(atom_xyz, k=max_neighbours+1, distance_upper_bound=upper_bound)
            mask = np.logical_and(dists>0.01, dists<np.inf)
            # assume first id is of that atom, i.e. dists[0]==0
            assert dists[0]==0, dists
            nn_ids[ids[0]] = ids[mask]
            #nn_dists[ids[0]] = dists[mask]
    
        jkls = {}
        for lid, nns in nn_ids.iteritems():
            if lid > max_id:
                continue
            jkls[lid] = []
            for nn in nns:
                # j is number of shared nearest neighbours
                common_nns = set(nn_ids[nn]).intersection(nns)
                j = len(common_nns)
                # k is number of bonds between nearest neighbours
                nn_bonds = []
                for common_nn in common_nns:
                    for nn_bond in set(nn_ids[common_nn]).intersection(common_nns):
                        if sorted((common_nn, nn_bond)) not in nn_bonds:
                            nn_bonds.append(sorted((common_nn, nn_bond)))
                k = len(nn_bonds)
                # l is longest chain of nearest neighbour bonds
                chain_lengths = [0]
                for nn_bond in nn_bonds:
                    chain_lengths.append(len(_chain(nn_bond, nn_bonds)))
                l = max(chain_lengths)
    
                jkls[lid].append('{0},{1},{2}'.format(j,k,l))
            
            jkls[lid] = Counter(jkls[lid])
    
        df['cna'] = [jkls[key] for key in sorted(jkls)]
        
        return df
        
    def cna_sum(self, atoms_df, upper_bound=4, max_neighbours=24,
                    repeat_vectors=None, leafsize=100):
        """ compute summed atomic environments of each atom in atoms_df
        
        Based on Faken, Daniel and Jónsson, Hannes,
        'Systematic analysis of local atomic structure combined with 3D computer graphics',
        March 1994, DOI: 10.1016/0927-0256(94)90109-0
        
        common signatures:
        - FCC = 12 x 4,2,1
        - HCP = 6 x 4,2,1 & 6 x 4,2,2
        - BCC = 6 x 6,6,6 & 8 x 4,4,4
        - Diamond = 12 x 5,4,3 & 4 x 6,6,3
        - Icosahedral = 12 x 5,5,5
        """
        df = self.common_neighbour_analysis(atoms_df, upper_bound, max_neighbours, 
                                            repeat_vectors, leafsize=leafsize)
        
        cnas = df.cna.values
        return sum(cnas,Counter())
