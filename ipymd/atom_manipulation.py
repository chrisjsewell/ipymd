# -*- coding: utf-8 -*-
"""
Created on Mon May 16 08:15:13 2016

@author: cjs14
"""
import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
from matplotlib import cm
from matplotlib.colors import Normalize
from six import string_types

from .shared import atom_data
from .shared.transformations import transform_to_crystal, rotate_vectors

class Atom_Manipulation(object):
    """ a class to manipulate atom data
    
    atom_df : pandas.DataFrame
        containing columns; x, y, z, type
    """
    def __init__(self, atom_df, meta_series=None,undos=1):
        """ a class to manipulate atom data
        
        atom_df : pandas.DataFrame
            containing columns; x, y, z, type
        meta_series : pandas.Series
            containing columns; origin, a, b, c to define unit cell
            if none it will be constructed from the min/max x, y, z values
        undos : int
            number of past dataframes to save
        """
        assert set(atom_df.columns).issuperset(['x','y','z','type'])
        assert undos > 0
        
        if meta_series is None:
            meta_series = pd.Series([(atom_df.x.min(),atom_df.y.min(),atom_df.z.min()),
                                     (atom_df.x.max()-atom_df.x.min(),0.,0.),
                                     (0.,atom_df.y.max()-atom_df.y.min(),0.),
                                     (0.,0.,atom_df.z.max()-atom_df.z.min())],
                                     index=['origin','a','b','c'])
        else:
            assert set(meta_series.index).issuperset(['origin','a','b','c'])                
        
        self._atom_df = atom_df.copy()
        #the one on the left is the newest
        self._old_atom_df = [None]*undos
        self._original_atom_df = atom_df.copy()

        self._meta = meta_series.copy()
        self._old_meta = [None]*undos
        self._original_meta = meta_series.copy()

    @property
    def df(self):
        return self._atom_df.copy() 
    @property
    def meta(self):
        return self._meta.copy() 
    
    def _save(self):
        self._old_atom_df.pop() # remove the oldest (right)
        self._old_atom_df.insert(0,self._atom_df.copy()) # add newest left

        self._old_meta.pop() # remove the oldest (right)
        self._old_meta.insert(0,self._meta.copy()) # add newest left
        
    def undo_last(self):
        if self._old_atom_df[0] is not None:
               self._atom_df = self._old_atom_df.pop(0) # get newest left
               self._old_atom_df.append(None)  

               self._meta = self._old_meta.pop(0) # get newest left
               self._old_meta.append(None)  
        else:
            raise Exception('No previous dataframes')
        
    def revert_to_original(self):
        """ revert to original atom_df """
        self._save()
        self._atom_df = self._original_atom_df.copy()
        self._meta = self._original_meta.copy()
        
    def change_variables(self, map_dict, vtype='type'):
        """ change particular variables according to the map_dict """
        self._save()
        self._atom_df.replace({vtype:map_dict}, inplace=True)

    def change_type_variable(self, atom_type, variable, value, type_col='type'):
        """ change particular variable for one atom type """
        self._save()
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
        
        self._save()
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
        
        self._save()
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
        
        self._save()
        self._atom_df.color = cats.map(color_dict)
                
    def filter_variables(self, values, vtype='type'):
        if isinstance(values, int):
            values = [values]
        if isinstance(values, float):
            values = [values]
        if isinstance(values, string_types):
            values = [values]
        
        self._save()
        self._atom_df = self._atom_df[self._atom_df[vtype].isin(values)]        
            
#------------------------
# Geometric manipulation
            
    def repeat_cell(self,a=1,b=1,c=1,original_first=False):
        """ repeat atoms along a, b, c directions (and update unit cell)
    
        a : int or tuple
            repeats in 'a' direction, if tuple then defines repeats in -/+ direction
        b : int or tuple
            repeats in 'b' direction, if tuple then defines repeats in -/+ direction
        c : int or tuple
            repeats in 'c' direction, if tuple then defines repeats in -/+ direction
        original_first: bool
            if True, the original atoms will be first in the DataFrame
        """
        if isinstance(a, int):
            arange = range(0,abs(a)+1)
        else:
            arange = range(-abs(a[0]),abs(a[1])+1)
        if isinstance(b, int):
            brange = range(0,abs(b)+1)
        else:
            brange = range(-abs(b[0]),abs(b[1])+1)
        if isinstance(c, int):
            crange = range(0,abs(c)+1)
        else:
            crange = range(-abs(c[0]),abs(c[1])+1)
        
        avec = np.asarray(self._meta.a)
        bvec = np.asarray(self._meta.b)
        cvec = np.asarray(self._meta.c)        

        dfs = []
        if original_first:
            dfs.append(self._atom_df.copy())            
            
        for i in arange:
            for j in brange:
                for k in crange:
                    if i==0 and j==0 and k==0 and original_first:
                        continue
                    atom_copy = self._atom_df.copy()
                    atom_copy[['x','y','z']] = (atom_copy[['x','y','z']]
                                + i*avec  + j*bvec + k*cvec)
                    dfs.append(atom_copy)
        self._save()
        self._atom_df = pd.concat(dfs)
        #TODO check for identical atoms and warn
        
        origin = np.asarray(self._meta.origin)
        self._meta.origin = tuple(origin + arange[0]*avec + brange[0]*bvec + crange[0]*cvec)
        self._meta.a = tuple(avec*len(arange))
        self._meta.b = tuple(bvec*len(brange))
        self._meta.c = tuple(cvec*len(crange))
        
#    def _convert_basis_abc(self, xyz):
#        """ convert  x,y,z to a,b,c basis
#        
#        xyz : numpy.array((N,3))
#        """
#        origin = np.asarray(self._meta.origin)
#        avec = np.asarray(self._meta.a)
#        bvec = np.asarray(self._meta.b)
#        cvec = np.asarray(self._meta.c)
#        
#        lattparams_bb(np.asarray([self._meta.a,self._meta.b,self._meta.c]))
#        
#        # move to origin
#        xyz = xyz - origin
#        
#        # convert basis
#        basis = np.array([avec,bvec,cvec]).T
#        invbasis = np.linalg.inv(basis)
#        
#        abc = np.dot(invbasis,xyz.T).T
#        
#        return abc


#    def _convert_basis_xyz(self, abc):
#        """ convert  a,b,c to x,y,z  basis
#        
#        abc : numpy.array((N,3))
#        """
#        origin = np.asarray(self._meta.origin)
#        avec = np.asarray(self._meta.a)
#        bvec = np.asarray(self._meta.b)
#        cvec = np.asarray(self._meta.c)
#        
#        # convert basis
#        basis = np.array([avec,bvec,cvec]).T
#        xyz = np.dot(basis,abc.T).T
#        
#        # move from origin
#        xyz = xyz + origin
#        
#        return xyz
            
        
    def slice_fraction(self, amin=0, amax=1, bmin=0, bmax=1, cmin=0, cmax=1,
                       incl_max=False, update_uc=True,delta=0.01):
        """ slice along a,b,c directions (from origin) as fraction of vector length 

        incl_max : bool
            whether to slice < (False) <= (True) max values
        update_uc : bool
            update unit cell (a,b,c,origin) to match slice
        delta : float
            retain atoms within 'delta' fraction outside of slice plane)
        
        """
        self._save()
        
        abc = transform_to_crystal(self._atom_df[['x','y','z']].values,
                                 self._meta.a,self._meta.b,self._meta.c,
                                 self._meta.origin)

        if incl_max:
            mask = ((abc[:,0]>=amin-delta) & (abc[:,0]<=amax+delta) & 
                    (abc[:,1]>=bmin-delta) & (abc[:,1]<=bmax+delta) & 
                    (abc[:,2]>=cmin-delta) & (abc[:,2]<=cmax+delta)) 
        else:
            mask = ((abc[:,0]>=amin-delta) & (abc[:,0]<amax+delta) & 
                    (abc[:,1]>=bmin-delta) & (abc[:,1]<bmax+delta) & 
                    (abc[:,2]>=cmin-delta) & (abc[:,2]<cmax+delta))  

        self._atom_df = self._atom_df[mask] 
        
        if update_uc:
            origin = np.asarray(self._meta.origin)
            avec = np.asarray(self._meta.a)
            bvec = np.asarray(self._meta.b)
            cvec = np.asarray(self._meta.c)
            
            self._meta.origin = tuple(origin+amin*avec+bmin*bvec+cmin*cvec)
            self._meta.a = avec * (amax-amin)
            self._meta.b = bvec * (bmax-bmin)
            self._meta.c = cvec * (cmax-cmin)            
        
    def slice_absolute(self, amin=0, amax=None, bmin=0, bmax=None, cmin=0, cmax=None,
                       incl_max=False, update_uc=True, delta=0.01):
        """ slice along a,b,c directions (from origin) given absolute vector length 

        if amax, bmax or cmax is None, then will use the vector length    
        
        update_uc : bool
            update unit cell (a,b,c,origin) to match slice
        incl_max : bool
            whether to slice < (False) <= (True) max values
        delta : float
            retain atoms within 'delta' fraction outside of slice plane)
            
        """
        self._save()
        
        abc = transform_to_crystal(self._atom_df[['x','y','z']].values,
                                 self._meta.a,self._meta.b,self._meta.c,
                                 self._meta.origin)
        
        anorm = np.linalg.norm(self._meta.a)
        bnorm = np.linalg.norm(self._meta.b)
        cnorm = np.linalg.norm(self._meta.c)
        
        amax = anorm if amax is None else amax
        bmax = anorm if bmax is None else bmax
        cmax = anorm if cmax is None else cmax                            

        if incl_max:       
            mask = ((abc[:,0]>=(amin/anorm)-delta) & (abc[:,0]<=(amax/anorm)+delta) & 
                    (abc[:,1]>=(bmin/bnorm)-delta) & (abc[:,1]<=(bmax/bnorm)+delta) & 
                    (abc[:,2]>=(cmin/cnorm)-delta) & (abc[:,2]<=(cmax/cnorm)+delta)) 
        else:
            mask = ((abc[:,0]>=(amin/anorm)-delta) & (abc[:,0]<(amax/anorm)+delta) & 
                    (abc[:,1]>=(bmin/bnorm)-delta) & (abc[:,1]<(bmax/bnorm)+delta) & 
                    (abc[:,2]>=(cmin/cnorm)-delta) & (abc[:,2]<(cmax/cnorm)+delta)) 

        self._atom_df = self._atom_df[mask]                                          
                                    
        if update_uc:
            origin = np.asarray(self._meta.origin)
            avec = np.asarray(self._meta.a)
            bvec = np.asarray(self._meta.b)
            cvec = np.asarray(self._meta.c)
            
            self._meta.origin = tuple(origin+
                                      amin*avec/anorm + 
                                      bmin*bvec/bnorm + 
                                      cmin*cvec/cnorm)
            self._meta.a = avec * (amax-amin)/anorm
            self._meta.b = bvec * (bmax-bmin)/bnorm
            self._meta.c = cvec * (cmax-cmin)/cnorm            

    def translate(self, vector, update_uc=True):
        """translate atoms by vector
        
        vector : list
            x, y, z translation
        update_uc : bool
            update unit cell (a,b,c,origin) to match translation
        
        """
        x,y,z = vector
        self._save()
        self._atom_df.x += x
        self._atom_df.y += y
        self._atom_df.z += z
        
        if update_uc:
            self._meta.origin = np.array(self._meta.origin) + np.asarray(vector)
            
    def rotate(self, angle, vector=[1,0,0], update_uc=True):
        """rotate the clockwise about the given axis vector 
        by theta degrees.
        
        e.g. for rotate_atoms(90,[0,0,1]); [0,1,0] -> [1,0,0]
        
        angle : float
            rotation angle in degrees
        vector : iterable
            vector to rotate around [x0,y0,z0] 
        update_uc : bool
            update unit cell (a,b,c,origin) to match rotation

        """
        coords = self._atom_df[['x','y','z']].values
        new_coords = rotate_vectors(coords,vector,angle)
        self._save()
        self._atom_df[['x','y','z']] = new_coords
        
        if update_uc:
            self._meta.a = tuple(rotate_vectors(self._meta.a,vector,angle)[0])
            self._meta.b = tuple(rotate_vectors(self._meta.b,vector,angle)[0])
            self._meta.c = tuple(rotate_vectors(self._meta.c,vector,angle)[0])
        
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
        self._save()
        inside = self._pnts_in_pointcloud(points, self._atom_df[['x','y','z']].values)
        self._save()
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
        df = self._atom_df.copy()
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
        self._save()
        self._atom_df = df