# -*- coding: utf-8 -*-
"""
Created on Mon May 16 08:15:13 2016

@author: cjs14
"""
import math
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib import cm
from matplotlib.colors import Normalize
from six import string_types

from .shared import atom_data

class Atom_Manipulation(object):
    """ a class to manipulate atom data
    
    atom_df : pandas.DataFrame
        containing columns; x, y, z, type
    """
    def __init__(self, atom_df,undos=1):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns; x, y, z, type
        undos : int
            number of past dataframes to save
        """
        assert set(atom_df.columns).issuperset(['x','y','z','type'])
        assert undos > 0
        
        self._atom_df = atom_df.copy()
        #the one on the left is the newest
        self._old_atom_df = [None]*undos
        self._original_atom_df = atom_df.copy()

    @property
    def df(self):
        return self._atom_df.copy() 
    
    def _save_df(self):
        self._old_atom_df.pop() # remove the oldest (right)
        self._old_atom_df.insert(0,self._atom_df.copy()) # add newest left
        
    def undo_last(self):
        if self._old_atom_df[0] is not None:
               self._atom_df = self._old_atom_df.pop(0) # get newest left
               self._old_atom_df.append(None)  
        else:
            raise Exception('No previous dataframes')
        
    def revert_to_original(self):
        """ revert to original atom_df """
        self._save_df()
        self._atom_df = self._original_atom_df.copy()
        
    def change_variables(self, map_dict, vtype='type'):
        """ change particular variables according to the map_dict """
        self._save_df()
        self._atom_df.replace({vtype:map_dict}, inplace=True)

    def change_type_variable(self, atom_type, variable, value, type_col='type'):
        """ change particular variable for one atom type """
        self._save_df()
        if not hasattr(value, '__iter__'):
            # there is df.loc, df.at or df.set_value, not sure which ones best
            self._atom_df.set_value(self._atom_df[type_col]==atom_type, variable, value)
        else:
            self._atom_df[variable] = self._atom_df[variable].astype(object) 
            df = self._atom_df
            # ensure all indexes are unique
            df.index.name = 'old_index'
            df.reset_index(drop=False,inplace=True)
            for indx in df.loc[df[type_col]==atom_type].index:
                df.set_value(indx,variable,value)
            self._atom_df = df.set_index('old_index')                
            
    def apply_map(self, vmap, column, default=False, type_col='type'):
        """ change values in a column, according to a mapping of another column

        Properties
        ----------
        vmap : dict or str
           A dictionary mapping values, or a string associated with a column in
           the ipymd.shared.atom_data() dataframe (e.g. color and RVdW)
        column : str
            the column to change
        default : various
            the default value to put when the type key cannot be found,
            if False then the original value will not be overwritten
           
        """
        if isinstance(vmap, string_types):
            df = atom_data()
            vmap = df[vmap].dropna().to_dict()
        
        self._atom_df = self._atom_df.copy()
        if default is not False:
            self._atom_df[column] = default
        for key, val in vmap.iteritems():
            self.change_type_variable(key, column, val, type_col)
        
    def color_by_index(self, cmap='jet', minv=None, maxv=None):
        """change colors to map index values 
        
        cmap : string
            the colormap to apply, see available at http://matplotlib.org/examples/color/colormaps_reference.html
        minv, maxv : float
            optional min, max cmap value, otherwise take min, max value found in column
            
        """
        colormap = cm.get_cmap(cmap)
        var = self._atom_df.index
        minval = var.min() if minv is None else minv
        maxval = var.max() if maxv is None else maxv
        norm = Normalize(minval, maxval,clip=True)
        
        self._save_df()
        self._atom_df.color = [tuple(col[:3]) for col in colormap(norm(var),bytes=True)]

    def color_by_variable(self, colname, cmap='jet', minv=None, maxv=None):
        """change colors to map 
        
        colname : string
            a coloumn of the dataframe that contains numbers by which to color
        cmap : string
            the colormap to apply, see available at http://matplotlib.org/examples/color/colormaps_reference.html
        minv, maxv : float
            optional min, max cmap value, otherwise take min, max value found in column
            
        """
        colormap = cm.get_cmap(cmap)
        var = self._atom_df[colname]
        minval = var.min() if minv is None else minv
        maxval = var.max() if maxv is None else maxv
        norm = Normalize(minval, maxval,clip=True)
        
        self._save_df()
        self._atom_df.color = [tuple(col[:3]) for col in colormap(norm(var),bytes=True)]

    def color_by_categories(self, colname, cmap='jet', sort=True):
        """change colors to map 
        
        colname : string
            a column of the dataframe that contains categories by which to color
        cmap : string
            the colormap to apply, see available at http://matplotlib.org/examples/color/colormaps_reference.html
            
        """
        colormap = cm.get_cmap(cmap)        
        cats = self._atom_df[colname]   
        unique_cats = cats.unique()
        if sort:
            unique_cats = sorted(unique_cats)
        num_cats = float(cats.nunique())
        
        #potential way to always have a string the same color
        #for cat in unique_cats:
            #[(ord(c)-97)/122. for c in cat.lower() if ord(c)>=97 and ord(c)<=122]        
        
        color_dict = dict([(cat,colormap(i/num_cats,bytes=True)[:3]) for i,cat in enumerate(unique_cats)])
        
        self._save_df()
        self._atom_df.color = cats.map(color_dict)
                
    def filter_variables(self, values, vtype='type'):
        if isinstance(values, int):
            values = [values]
        if isinstance(values, float):
            values = [values]
        if isinstance(values, string_types):
            values = [values]
        
        self._save_df()
        self._atom_df = self._atom_df[self._atom_df[vtype].isin(values)]        
        
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
    
    def filter_inside_pts(self, points):
        """return only atoms inside the bounding shape of a set of points 

        points : np.array((N,3))        
        """ 
        inside = self._pnts_in_pointcloud(points, self._atom_df[['x','y','z']].values)
        self._save_df()
        self._atom_df = self._atom_df[inside]

    def filter_inside_box(self, vectors, origin=np.zeros(3)):
        """return only atoms inside box
        
        vectors : np.array((3,3))
            a, b, c vectors
        origin : np.array((1,3))
        
        """
        a,b,c = vectors + origin
        points = [origin, a, b, a+b, c, a+c, b+c, a+b+c]
        self.filter_inside_pts(points)
    
    def _rotate(self, vector, axis, theta):
        """rotate the vector v clockwise about the given axis vector 
        by theta degrees.
        
        e.g. df._rotate([0,1,0],[0,0,1],90) -> [1,0,0]
        
        vector : iterable or list of iterables
            vector to rotate [x,y,z] or [[x1,y1,z1],[x2,y2,z2]]
        axis : iterable
            axis to rotate around [x0,y0,z0] 
        theta : float
            rotation angle in degrees
        """
        theta = -1*theta
        
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
        
        #TODO think about using np.einsum or something more efficient?
        if len(np.asarray(vector).shape)==1:
            return np.dot(rotation_matrix, vector)
        else:
            new_vectors = []
            for v in vector:
               new_vectors.append(np.dot(rotation_matrix, v)) 
            return np.asarray(new_vectors)

    def filter_inside_hexagon(self, vectors, origin=np.zeros(3)):
        """return only atoms inside hexagonal prism
        
        vectors : np.array((2,3))
            a, c vectors
        origin : np.array((1,3))
        
        """
        a, c = vectors
        points = [self._rotate(a, c, angle) for angle in [0,60,120,180,240,300]]
        points += [p + c for p in points]
        points = np.array(points) + origin
        self.filter_inside_pts(points)

    def repeat_cell(self, vectors, repetitions=((0,1),(0,1),(0,1)),original_first=False):
        """ repeat atoms along vectors a, b, c 

        vectors : np.array((3,3))
            a,b,c vectors  
        original_first: bool
            if True, the original atoms will be first in the DataFrame
        """
        xreps,yreps,zreps = repetitions
        if isinstance(xreps, int):
            xreps = (0,xreps)
        if isinstance(yreps, int):
            yreps = (0,yreps)
        if isinstance(zreps, int):
            zreps = (0,zreps)

        dfs = []
        if original_first:
            dfs.append(self._atom_df.copy())            
            
        for i in range(xreps[0], xreps[1]+1):
            for j in range(yreps[0], yreps[1]+1):
                for k in range(zreps[0], zreps[1]+1):
                    if i==0 and j==0 and k==0 and original_first:
                        continue
                    atom_copy = self._atom_df.copy()
                    atom_copy[['x','y','z']] = (atom_copy[['x','y','z']]
                                + i*vectors[0]  + j*vectors[1] + k*vectors[2])
                    dfs.append(atom_copy)
        self._save_df()
        self._atom_df = pd.concat(dfs)
        #TODO check for identical atoms and warn
        
    def slice_x(self, minval=None, maxval=None):
        self._save_df()
        if minval is not None:
            self._atom_df = self._atom_df[self._atom_df['x']>=minval].copy()
        if maxval is not None:
            self._atom_df = self._atom_df[self._atom_df['x']<=maxval].copy()

    def slice_y(self, minval=None, maxval=None):
        self._save_df()
        if minval is not None:
            self._atom_df = self._atom_df[self._atom_df['y']>=minval].copy()
        if maxval is not None:
            self._atom_df = self._atom_df[self._atom_df['y']<=maxval].copy()

    def slice_z(self, minval=None, maxval=None):
        self._save_df()
        if minval is not None:
            self._atom_df = self._atom_df[self._atom_df['z']>=minval].copy()
        if maxval is not None:
            self._atom_df = self._atom_df[self._atom_df['z']<=maxval].copy()
                                    
    #TODO slice along arbitrary direction
                                    
    def translate_atoms(self, vector):
        """translate atoms by vector
        
        vector : list
            x, y, z translation
        
        """
        x,y,z = vector
        self._save_df()
        self._atom_df.x += x
        self._atom_df.y += y
        self._atom_df.z += z
    
    def rotate_atoms(self, angle, vector=[1,0,0]):
        """rotate the clockwise about the given axis vector 
        by theta degrees.
        
        e.g. for rotate_atoms(90,[0,0,1]); [0,1,0] -> [1,0,0]
        
        angle : float
            rotation angle in degrees
        vector : iterable
            vector to rotate around [x0,y0,z0] 

        """
        xyz = self._atom_df[['x','y','z']].values
        new_xyz = self._rotate(xyz, vector,angle)
        self._save_df()
        self._atom_df[['x','y','z']] = new_xyz
        
    def group_atoms_as_mols(self, atom_ids, name, remove_atoms=True, mean_xyz=True,
                            color='red',transparency=1.,radius=1.):
        """ combine atoms into a molecule
        atom_ids : list of lists
            list of dataframe indexes for each molecule
        name : string
            name of molecule
        remove_atoms : bool
            remove the grouped atoms from the dataframe
        mean_xyz : bool
            use the mean coordinate of atoms for molecule, otherwise use coordinate of first atom
            
        """
        # test no atoms in multiple molecules
        all_atoms = np.asarray(atom_ids).flatten()  
        assert len(set(all_atoms))==len(all_atoms), 'atoms in multiple molecules'

        mol_data = []
        for atoms in atom_ids:
            if mean_xyz:
                x,y,z = df.loc[atoms,['x','y','z']].mean().values
            else:
                x,y,z = df.loc[atoms[0],['x','y','z']].values
            mol_data.append([name,x,y,z,radius,color,transparency])
            
            if remove_atoms:
                df.drop(atoms, inplace=True)

        if mol_data:
            moldf = pd.DataFrame(mol_data,columns=['type','x','y','z','radius','color','transparency'])
            df = pd.concat([df,moldf])
        self._save_df()
        self._atom_df = df