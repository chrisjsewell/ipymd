# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:19:07 2016

@author: cjs14
"""

#import math
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull, cKDTree
from collections import Counter
from IPython.core.display import clear_output
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from .atom_manipulation import Atom_Manipulation

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2, rounded=None):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle =  np.degrees(np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)))
    if rounded is not None:
        angle = round(angle,rounded)
    return angle
    
def _createTreeFromEdges(edges):
    """    
    e.g. _createTreeFromEdges([[1,2],[0,1],[2,3],[8,9],[0,3]])
     -> {0: [1], 1: [2, 0], 2: [1, 3], 3: [2,0], 8: [9], 9: [8]}
    """
    tree = {}
    for v1, v2 in edges:
        tree.setdefault(v1, []).append(v2)
        tree.setdefault(v2, []).append(v1)
    return tree
    
def _longest_path(start,tree,lastnode=None):
    """a recursive function to compute the maximum unbroken chain given a tree
    
    e.g. start=0, tree={0: [1], 1: [2, 0], 2: [1, 3], 3: [2,0], 8: [9], 9: [8]}
     -> [0, 1, 2, 3, 0]
    
    """
    if not tree.has_key(start):
        return []
    new_tree = tree.copy()
    #nodes = new_tree.pop(start) # can use if don't want to complete loops
    nodes = new_tree[start]
    new_tree[start] = []
        
    path = []
    for node in nodes:        
        if node==lastnode:
            continue # can't go back to lastnode, e.g. 1->2->1
        new_path = _longest_path(node,new_tree,start)
        if len(new_path) > len(path):
            path = new_path
    path.append(start)
    return path    

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
    
    def calc_volume_bb(self, vectors=np.array([[1,0,0],[0,1,0],[0,0,1]]), rounded=None):
        """ calculate volume of the bounding box        
        """
        a,b,c = vectors
        vol = a.dot(np.cross(b,c))
        if rounded is not None:
            vol = round(vol, rounded)
        return vol

    def calc_lattparams_bb(self, vectors=np.array([[1,0,0],[0,1,0],[0,0,1]]),
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

        alpha = angle_between(b,c, rounded)
        beta = angle_between(a,c, rounded)
        gamma = angle_between(a,b, rounded)
        
        return a_length, b_length, c_length, alpha, beta, gamma
        
        
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
    
    def calc_bond_lengths(self,atoms_df, coord_type, lattice_type, max_dist=4, max_coord=16,
                          repeat_vectors=None, rounded=2, min_dist=0.01, leafsize=100):
        """ calculate the unique bond lengths atoms in coords_atoms, w.r.t lattice_atoms
        
        atoms_df : pandas.Dataframe
            all atoms
        coord_type : string
            atoms to calcualte coordination of
        lattice_type : string
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
        distances : set
            list of unique distances
        
        """
        coord_df = Atom_Manipulation(atoms_df)
        coord_df.filter_variables(coord_type)
       
        lattice_df = Atom_Manipulation(atoms_df)
        lattice_df.filter_variables(lattice_type)
        
        if repeat_vectors is not None:
            lattice_df.repeat_cell(repeat_vectors,((-1,1),(-1,1),(-1,1)))

        lattice_tree = cKDTree(lattice_df.df[['x','y','z']].values, leafsize=leafsize)
        all_dists,all_ids = lattice_tree.query(coord_df.df[['x','y','z']].values, k=max_coord, distance_upper_bound=max_dist)
        
        distances = []
        for dists in all_dists:
            for d in dists:
                if d > min_dist and not np.isinf(d):
                    distances.append(round(d,rounded))
            
        return sorted(set(distances))

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
        all_dists,all_ids = lattice_tree.query(coord_atoms_df[['x','y','z']].values, k=max_coord, distance_upper_bound=max_dist)
        
        coords = []
        for dists in all_dists:
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
        lattice_type : string
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
        dists,idnums = lattice_tree.query(atoms_df[['x','y','z']].values, k=1, distance_upper_bound=max_dist)
        return dists

    def vacancy_identification(self, atoms_df, res=0.2, nn_dist=2., repeat_vectors=None, remove_dups=True,
                 color='red',transparency=1.,radius=1, type_name='Vac', leafsize=100, 
                 n_jobs=1, ipython_progress=False, ):
            """ identify vacancies
            
            atoms_df : pandas.Dataframe
                atoms to calculate for
            res : float
                resolution of vacancy identification, i.e. spacing of reference lattice
            nn_dist : float
                maximum nearest-neighbour distance considered as a vacancy 
            repeat_vectors : np.array((3,3))
                include consideration of repeating boundary idenfined by a,b,c vectors
            remove_dups : bool
                only keep one vacancy site within the nearest-neighbour distance
            leafsize : int
                points at which the algorithm switches to brute-force (kdtree specific)
            n_jobs : int, optional
                Number of jobs to schedule for parallel processing. If -1 is given all processors are used. 
            ipython_progress : bool
                print progress to IPython Notebook
            
            Returns
            -------
            vac_df : pandas.DataFrame
                new atom dataframe of vacancy sites as atoms
            
            """
            xmin, xmax = atoms_df.x.min(),atoms_df.x.max()
            ymin, ymax = atoms_df.y.min(),atoms_df.y.max()
            zmin, zmax = atoms_df.z.min(),atoms_df.z.max()
            xyz = np.mgrid[xmin:xmax:res, ymin:ymax:res, zmin:zmax:res].reshape(3,-1).T
    
            if repeat_vectors is not None:
                repeat = Atom_Manipulation(atoms_df)
                repeat.repeat_cell(repeat_vectors,((-1,1),(-1,1),(-1,1)),original_first=True)
                lattice_df = repeat.df
            else:
                lattice_df = atoms_df

            if ipython_progress:
                clear_output()
                print('creating nearest neighbour tree')
            
            lattice_tree = cKDTree(lattice_df[['x','y','z']].values, leafsize=leafsize)

            if ipython_progress:
                clear_output()
                print('assessing nearest neighbours')

            dists,idnums = lattice_tree.query(xyz, k=1, distance_upper_bound=nn_dist,n_jobs=n_jobs)
    
            vac_list = []
            for atom,dist in zip(xyz,dists):
                if np.isinf(dist):
                    x,y,z = atom
                    vac_list.append([type_name,x,y,z,radius,color,transparency])
                    
            df = pd.DataFrame(vac_list,columns=['type','x','y','z','radius','color','transparency'])
            
            if remove_dups and df.shape[0]>0:
                vac_tree = cKDTree(df[['x','y','z']].values)
                pairs = np.asarray(list(vac_tree.query_pairs(nn_dist)))
                #drop first atom of each pair
                if pairs.shape[0] > 0:
                    df.drop(pairs[:,0],inplace=True)

            if ipython_progress:
                clear_output()
            
            return df
        
    #TODO group atoms into specified molecules e.g. S2 or CaCO3
    # http://chemwiki.ucdavis.edu/Textbook_Maps/Inorganic_Chemistry_Textbook_Maps/Map%3A_Inorganic_Chemistry_(Wikibook)/Chapter_08%3A_Ionic_and_Covalent_Solids_-_Structures/8.2%3A_Structures_related_to_NaCl_and_NiAs
    # maybe supply central atom type(s) and 'other' atoms type(s), filter df by required atom types, 
    # then find nearest neighbours of central (removing molecule each time)
    # create molecule x,y,z from average of central atoms
            
    #http://www.ovito.org/manual/particles.modifiers.common_neighbor_analysis.html
    #https://www.quora.com/Given-a-set-of-atomic-types-and-coordinates-from-an-MD-simulation-is-there-a-good-algorithm-for-determining-its-likely-crystal-structure?__filter__=all&__nsrc__=2&__snid3__=179254150
    # http://iopscience.iop.org/article/10.1088/0965-0393/20/4/045021/pdf            
    def common_neighbour_analysis(self, atoms_df, upper_bound=4, max_neighbours=24,
                                  repeat_vectors=None, leafsize=100, ipython_progress=False):
        """ compute atomic environment of each atom in atoms_df
        
        Based on Faken, Daniel and Jónsson, Hannes,
        'Systematic analysis of local atomic structure combined with 3D computer graphics',
        March 1994, DOI: 10.1016/0927-0256(94)90109-0
        
        ideally:
        - FCC = 12 x 4,2,1
        - HCP = 6 x 4,2,1 & 6 x 4,2,2
        - BCC = 6 x 6,6,6 & 8 x 4,4,4
        - icosahedral = 12 x 5,5,5
        
        Paramaters
        ----------
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors
        ipython_progress : bool
            print progress to IPython Notebook

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

        if ipython_progress:
            print('creating nearest neighbours dictionary')
        
        # create nearest neighbours dictionary
        lattice_tree = cKDTree(lattice_df[['x','y','z']].values, leafsize=leafsize)
        all_dists,all_ids = lattice_tree.query(lattice_df[['x','y','z']].values, k=max_neighbours+1, distance_upper_bound=upper_bound)
        
        nn_ids = {}
        #nn_dists = {}
        for dists,ids in zip(all_dists,all_ids):
            
            mask = np.logical_and(dists>0.01, dists<np.inf)
            # assume first id is of that atom, i.e. dists[0]==0
            assert dists[0]==0, dists
            nn_ids[ids[0]] = ids[mask]
            #nn_dists[ids[0]] = dists[mask]
            
        jkls = {}
        for lid, nns in nn_ids.iteritems():
            if lid > max_id:
                continue
            if ipython_progress:
                clear_output()
                print('assessing nearest neighbours: {0} of {1}'.format(lid,max_id))
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
                tree = _createTreeFromEdges(nn_bonds)
                chain_lengths = [0]
                for node in tree.iterkeys():
                    chain_lengths.append(len(_longest_path(node, tree))-1)
                l = max(chain_lengths)
    
                jkls[lid].append('{0},{1},{2}'.format(j,k,l))
            
            jkls[lid] = Counter(jkls[lid])
    
        df['cna'] = [jkls[key] for key in sorted(jkls)]
        
        if ipython_progress:
            clear_output()
        
        return df
        
    
    def _equala(self,i, j, accuracy):
        return j*accuracy <= i <= j+j*(1-accuracy)
        
    def cna_categories(self, atoms_df, accuracy=1., upper_bound=4, max_neighbours=24,
                    repeat_vectors=None, leafsize=100, ipython_progress=False):
        """ compute summed atomic environments of each atom in atoms_df
        
        Based on Faken, Daniel and Jónsson, Hannes,
        'Systematic analysis of local atomic structure combined with 3D computer graphics',
        March 1994, DOI: 10.1016/0927-0256(94)90109-0
        
        signatures:
        - FCC = 12 x 4,2,1
        - HCP = 6 x 4,2,1 & 6 x 4,2,2
        - BCC = 6 x 6,6,6 & 8 x 4,4,4
        - Diamond = 12 x 5,4,3 & 4 x 6,6,3
        - Icosahedral = 12 x 5,5,5
        
        Parameters
        ----------
        accuracy : float
            0 to 1 how accurate to fit to signature
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors
        ipython_progress : bool
            print progress to IPython Notebook

        Returns
        -------
        df : pandas.Dataframe
            copy of atoms_df with new column named cna

        """
        df = self.common_neighbour_analysis(atoms_df, upper_bound, max_neighbours, 
                                            repeat_vectors, leafsize=leafsize, 
                                            ipython_progress=ipython_progress)
        
        cnas = df.cna.values
        
        atype = []
        for counter in cnas:
            if self._equala(counter['4,2,1'],6,accuracy) and self._equala(counter['4,2,2'],6,accuracy):
                atype.append('HCP')
            elif self._equala(counter['4,2,1'],12,accuracy):
                atype.append('FCC')
            elif self._equala(counter['6,6,6'],8,accuracy) and self._equala(counter['4,4,4'],6,accuracy):
                atype.append('BCC')
            elif self._equala(counter['5,4,3'],12,accuracy) and self._equala(counter['6,6,3'],4,accuracy):
                atype.append('Diamond')
            elif self._equala(counter['5,5,5'],12,accuracy):
                atype.append('Icosahedral')
            else:
                atype.append('Other')
        df.cna = atype
        return df

    def cna_sum(self, atoms_df, upper_bound=4, max_neighbours=24,
                    repeat_vectors=None, leafsize=100, ipython_progress=False):
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

        Parameters
        ----------
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors
        ipython_progress : bool
            print progress to IPython Notebook

        Returns
        -------
        counter : Counter
            a counter of cna signatures

        """
        df = self.common_neighbour_analysis(atoms_df, upper_bound, max_neighbours, 
                                            repeat_vectors, leafsize=leafsize, 
                                            ipython_progress=ipython_progress)
        
        cnas = df.cna.values
        return sum(cnas,Counter())
    
    def cna_plot(self, atoms_df, upper_bound=4, max_neighbours=24,
                    repeat_vectors=None, leafsize=100, 
                    barwidth=1, ipython_progress=False):
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

        Parameters
        ----------
        repeat_vectors : np.array((3,3))
            include consideration of repeating boundary idenfined by a,b,c vectors
        ipython_progress : bool
            print progress to IPython Notebook

        Returns
        -------
        plot : matplotlib.pyplot
            a matplotlib plot

        """
        df = self.common_neighbour_analysis(atoms_df, upper_bound, max_neighbours, 
                                            repeat_vectors, leafsize=leafsize, 
                                            ipython_progress=ipython_progress)
        
        cnas = df.cna.values
        counter = sum(cnas,Counter())
        
        labels, values = zip(*counter.items())
        indexes = np.arange(len(labels))

        colors = []
        patches = []
        d = {'4,2,1':['orange','FCC or HCP (1 of 2)'],
             '4,2,2':['red','HCP (1 of 2)'],
             '6,6,6':['green','BCC (1 of 2)'],
             '4,4,4':['green','BCC (2 of 2)'],
            '5,5,5':['purple','Icosahedral'],
             '5,4,3':['grey','Diamond (1 of 2)'],
             '6,6,3':['grey','Diamond (1 of 2)']}
        for label in labels:
            if d.has_key(label):
                colors.append(d[label][0])
                patches.append(mpatches.Patch(color=d[label][0], label=d[label][1]))
            else:
                colors.append('blue')
               
        plt.barh(indexes, values, barwidth, color=colors)
        plt.yticks(indexes + barwidth * 0.5, labels)
        plt.grid(True)
        if patches:
            plt.legend(handles=patches)

        plt.ylabel('i,j,k')
        
        return plt

    #TODO _group_molecules needs work
    def _group_molecules(self,atom_df,moltypes,maxdist=3,repeat_vectors=None,
                        mean_xyz=True,remove_atoms=True,
                        color='red',transparency=1.,radius=1.,
                        leafsize=100):
        molname = ''.join(['{}_{}'.format(k,v) 
                        for k,v in Counter(moltypes).iteritems()])
    
        search_df = atom_df[atom_df.type.isin([moltypes[0]])].copy()
        old_index = search_df.index
        search_df.reset_index(inplace=True)
    
        if repeat_vectors is not None:
            manip = Atom_Manipulation(search_df)
            manip.repeat_cell(repeat_vectors,repetitions=((-1,1), (-1,1), (-1, 1)),original_first=True)
            lattice_df = manip.df
            lattice_df.reset_index(inplace=True,drop=True)
            rep_map = dict(zip(range(lattice_df.shape[0]),list(search_df.index)*27))
        else:
            lattice_df = search_df.copy()    
            rep_map = dict(zip(range(search_df.shape[0]),list(search_df.index)))
    
        lattice_tree = cKDTree(lattice_df[['x','y','z']].values, leafsize=leafsize)
        dists,idnums = lattice_tree.query(search_df[['x','y','z']].values, k=len(moltypes), distance_upper_bound=maxdist)
    
        mol_data = []
        used_repeat = []
        for i,dist,idnum in zip(search_df.index, dists,idnums):
            #print i, old_index[i], dist, [rep_map[m] for m in idnum]
            if i in used_repeat:
                continue    
            mol = [i]
            for j, d in zip(idnum[1:],dist[1:]):
                if not j in mol and not np.isinf(d) and not rep_map.get(j) in used_repeat:
                    #used_repeat.append(rep_map[j])
                    mol.append(j)
                else:
                    print 'warning incomplete molecule'#, idnum, dist, j, rep_map.get(j)
            
            repeat_mol = [rep_map[m] for m in mol]
            if len(mol) == len(moltypes):
                
                used_repeat.extend(repeat_mol)
                if mean_xyz:
                    x,y,z = search_df.loc[repeat_mol,['x','y','z']].mean().values
                else:
                    x,y,z = search_df.loc[i,['x','y','z']].values
                mol_data.append([molname,x,y,z,radius,color,transparency])
            #print i, repeat_mol
    
        df = atom_df.copy()
        if remove_atoms:
            df.drop(used_repeat, inplace=True)
        if mol_data:
            moldf = pd.DataFrame(mol_data,columns=['type','x','y','z','radius','color','transparency'])
            df = pd.concat([df,moldf])
        return df