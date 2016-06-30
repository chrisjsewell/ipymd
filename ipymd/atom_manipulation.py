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

from chemlab.db import ChemlabDB
#have to convert from nm to angstrom
vdw_dict = dict((k,v*10) for k,v in ChemlabDB().get("data", 'vdwdict').iteritems())

# Maps
default_atom_map = {
    "C": "gray",
    "O": "red",
    "H": "white",

    "N": "light_blue",
    "S": "gold",
    "Cl": "green",
    "B": " green",
    
    "P": "orange",
    "Fe": "orange",
    "Ba": "orange",

    "Na": "blue",
    "Mg": "forest_green",
    
    "Zn": "brown",
    "Cu": "brown",
    "Ni": "brown",
    "Br": "brown",

    "Ca": "dark_gray",
    "Mn": "dark_gray",
    "Al": "dark_gray",
    "Ti": "dark_gray",
    "Cr": "dark_gray",
    "Ag": "dark_gray",

    "F": " goldenrod",    
    "Si": "goldenrod",
    "Au": "goldenrod",
    
    "I": "purple",
        
    "Li": "fire_brick",
    "He": "pink",

    "Xx": "deep_pink",
}


light_atom_map = {
    "C": "gainsboro",
    "O": "light_salmon",
    "H": "snow",

"N": " pale_turquoise",
    "S": " light_goldenrod_yellow",
    "Cl": "pale_green",
    "B": " pale_green",
    
"P": "beige",
    "Fe": "beige",
    "Ba": "beige",

"Na": "lavender",
    "Mg": "aquamarine",
    
"Zn": "dark_salmon",
    "Cu": "dark_salmon",
    "Ni": "dark_salmon",
    "Br": "dark_salmon",

"Ca": "light_slate_gray",
    "Mn": "light_slate_gray",
    "Al": "light_slate_gray",
    "Ti": "light_slate_gray",
    "Cr": "light_slate_gray",
    "Ag": "light_slate_gray",

"F": " pale_goldenrod",    
    "Si": "pale_goldenrod",
    "Au": "pale_goldenrod",
    
"I": "lavender",
    
"Li": "light_coral",
    "He": "light_pink",

"Xx": "deep_pink",
}


class Atom_Manipulation(object):
    """ a class to manipulate atom data
    
    atom_df : pandas.DataFrame
        containing columns; x, y, z, type
    """
    def __init__(self, atom_df):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns; x, y, z, type
        """
        assert set(atom_df.columns).issuperset(['x','y','z','type'])
        
        self._atom_df_new = atom_df.copy()
        self._atom_df_old = None
        self._original_atom_df = atom_df.copy()

    @property
    def df(self):
        return self._atom_df_new.copy()    
    
    @property
    def _atom_df(self):
        return self._atom_df_new   

    @_atom_df.setter
    def _atom_df(self, atom_df):
        self._atom_df_old = self._atom_df_new
        self._atom_df_new = atom_df
    
    def undo_last(self):
        if self._atom_df_old is not None:
            self._atom_df_new = self._atom_df_old
            self._atom_df_old = None
            
        
    def revert_to_original(self):
        """ revert to original atom_df """
        self._atom_df = self._original_atom_df.copy()
        
    def change_variables(self, map_dict, vtype='type'):
        """ change particular variables """
        self._atom_df.replace({vtype:map_dict}, inplace=True)

    def change_type_variable(self, atom_type, variable, value):
        """ change particular variable for one atom type """
        self._atom_df.loc[self._atom_df['type']==atom_type, variable] = value

    def apply_colormap(self, colormap=default_atom_map):
        """
        colormap : dict
           A dictionary mapping atom types to colors, in str format 
           By default it is a default color scheme (light_atom_map also available).   
        """
        for key, val in colormap.iteritems():
            self.change_type_variable(key, 'color', val)
    
    def color_by_variable(self, variable, cmap='jet', minv=None, maxv=None):
        """change colors to map 
        
        variable : string
            a coloumn of the dataframe that contains numbers by which to color
        cmap : string
            the colormap to apply, see available at http://matplotlib.org/examples/color/colormaps_reference.html
        minv, maxv : float
            optional min, max cmap value, otherwise take min, max value found in column
            
        """
        colormap = cm.get_cmap(cmap)
        var = self._atom_df[variable]
        minval = var.min() if minv is None else minv
        maxval = var.max() if maxv is None else maxv
        norm = Normalize(minval, maxval,clip=True)
        
        self._atom_df.color = [tuple(col[:3]) for col in colormap(norm(var),bytes=True)]

    def apply_radiimap(self, radiimap=vdw_dict):
        """
        radii_map: dict
           A dictionary mapping atom types to radii. The default is the
           mapping contained in `chemlab.db.vdw.vdw_dict`
        
        """
        for key, val in radiimap.iteritems():
            self.change_type_variable(key, 'radius', val)
        
    def filter_variables(self, values, vtype='type'):
        if isinstance(values, int):
            values = [values]
        if isinstance(values, float):
            values = [values]
        if isinstance(values, basestring):
            values = [values]
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
        self._atom_df = pd.concat(dfs)
        #TODO check for identical atoms and warn
        
    def slice_x(self, minval=None, maxval=None):
        if minval is not None:
            self._atom_df = self._atom_df[self._atom_df['x']>=minval].copy()
        if maxval is not None:
            self._atom_df = self._atom_df[self._atom_df['x']<=maxval].copy()

    def slice_y(self, minval=None, maxval=None):
        if minval is not None:
            self._atom_df = self._atom_df[self._atom_df['y']>=minval].copy()
        if maxval is not None:
            self._atom_df = self._atom_df[self._atom_df['y']<=maxval].copy()

    def slice_z(self, minval=None, maxval=None):
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
        self._atom_df = self._atom_df.copy()
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
        self._atom_df = self._atom_df.copy()
        xyz = self._atom_df[['x','y','z']].values
        new_xyz = self._rotate(xyz, vector,angle)
        self._atom_df[['x','y','z']] = new_xyz
        

