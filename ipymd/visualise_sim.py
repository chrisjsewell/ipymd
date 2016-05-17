# -*- coding: utf-8 -*-
"""
Created on Sun May  1 23:47:03 2016

@author: cjs14
"""
from io import BytesIO

import numpy as np

from chemlab.graphics.qtviewer import QtViewer
from chemlab.graphics.renderers.line import LineRenderer
#from chemlab.graphics.postprocessing import (FXAAEffect, GammaCorrectionEffect, 
#                                             OutlineEffect, SSAOEffect)
from chemlab.graphics.colors import get as str_to_colour
from chemlab.graphics import colors as chemlab_colors
from chemlab.db import ChemlabDB
from chemlab.graphics.transformations import rotation_matrix
def orbit_z(self, angle):
    # Subtract pivot point
    self.position -= self.pivot        
    # Rotate
    rot = rotation_matrix(-angle, self.c)[:3,:3]
    self.position = np.dot(rot, self.position)        
    # Add again the pivot point
    self.position += self.pivot
    
    self.a = np.dot(rot, self.a)
    self.b = np.dot(rot, self.b)
    self.c = np.dot(rot, self.c)     
from chemlab.graphics.camera import Camera
Camera.orbit_z = orbit_z

from IPython.display import Image as ipy_Image
from PIL import Image, ImageChops

# in order to set atoms as transparent
from .chemlab_patch.atom import AtomRenderer
from .chemlab_patch.triangle import TriangleRenderer
from .chemlab_patch.box import BoxRenderer
from .chemlab_patch.hexagon import HexagonRenderer

class Visualise_Sim(object):
    """ 
    A class to visualise atom data    
    """
    _unit_dict = {'real':{'distance':0.1}}
    
    def __init__(self, colormap=None, radiimap=None, units='real'):
        """
        colormap: dict, should contain the 'Xx' key,value pair
           A dictionary mapping atom types to colors. By default it is the color
           scheme provided by `chemlab.graphics.colors.default_atom_map`. The 'Xx'
           symbol value is taken as the default color.        
        radii_map: dict, should contain the 'Xx' key,value pair.
           A dictionary mapping atom types to radii. The default is the
           mapping contained in `chemlab.db.vdw.vdw_dict`        
        For units *real*, these are the units:
        
            mass = grams/mole
            distance = Angstroms
            time = femtoseconds
            energy = Kcal/mole
            velocity = Angstroms/femtosecond
            force = Kcal/mole-Angstrom
            torque = Kcal/mole
            temperature = Kelvin
            pressure = atmospheres
            dynamic viscosity = Poise
            charge = multiple of electron charge (1.0 is a proton)
            dipole = charge*Angstroms
            electric field = volts/Angstrom
            density = gram/cm^dim
        """
        assert units=='real', 'currently only supports real units'
        self._units = units
        self._atomcolors = None
        self.change_atom_colormap()
        self._atomradii = None  
        self.change_atom_radiimap()
        
        # rendered objects
        self._atoms = []
        self._boxes = []
        self._hexagons = []
        self._axes = None
        self._triangles = []
        
    def remove_all_objects(self):
        self._atoms = []
        self._boxes = []
        self._hexagons = []
        self._axes = None
        self._triangles = []                
        
    def change_atom_colormap(self, colormap=None, colorstrs=False):
        """
        colormap : dict, should contain the 'Xx' key,value pair
           A dictionary mapping atom types to colors, in RGBA format (0-255) 
           By default it is the color scheme provided by 
           `chemlab.graphics.colors.default_atom_map`. The 'Xx' symbol value is 
           taken as the default color.   
        colorstrs : bool
            if True, colors should be strings matching colors in `chemlab.graphics.colors`
        """
        if colormap is None:
            self._atomcolors = chemlab_colors.default_atom_map
        else:
            assert colormap.has_key('Xx'), "colormap should contain an 'Xx' default key"
            if colorstrs:
                colormap = dict([[k,str_to_colour(v)] for k,v in colormap.iteritems()])
                if None in colormap.values():
                    raise ValueError('one or more colors not found')
            self._atomcolors = colormap

    def change_atom_radiimap(self, radiimap=None):
        """
        radii_map: dict, should contain the 'Xx' key,value pair.
           A dictionary mapping atom types to radii. The default is the
           mapping contained in `chemlab.db.vdw.vdw_dict`        
        """
        if radiimap is None:
            self._atomradii = ChemlabDB().get("data", 'vdwdict')
        else:
            assert radiimap.has_key('Xx'), "radiimap should contain an 'Xx' default key"
            self._atomradii = radiimap
        
    def add_atoms(self, atoms_df, type_map={}, spheres=True, alpha=1., shading='phong'):
        """ add atoms to visualisation

        atoms_df : pandas.DataFrame
            a table of atom data, must contain columns;  xs, yx, zs and type
        type_map : dict
            mapping of types
        spheres : bool
            whether the atoms are rendered as spheres or points
        alpha : float
            how transparent the atoms are (if spheres) 0 to 1
        shading : str
            phong or toon
        """
        assert set(['xs','ys','zs','type']).issubset(set(atoms_df.columns))
        
        r_array = np.array(atoms_df[['xs','ys','zs']])
        r_array = self._unit_conversion(r_array, 'distance')
        
        type_array = atoms_df['type'].map(lambda x: type_map.get(x,x)).tolist()
        
        backend = 'impostors' if spheres else 'points'
        
        assert alpha <= 1. and alpha > 0., 'alpha must be between 0 and 1'

        self._atoms.append([r_array.copy(), type_array, backend, alpha, shading])   
        
    def remove_atoms(self, n=1):
        """ remove the last n sets of atoms to be added """
        assert len(self._atoms) >= n
        self._atoms = self._atoms[:-n]
        
    def add_box(self, vectors, origin=np.zeros(3), color='black', width=1):
        """ add wireframed box to visualisation
        
       vectors : np.ndarray((3,3), dtype=float)
          The three vectors representing the sides of the box.
       origin : np.ndarray((3,3), dtype=float), default to zero
          The origin of the box.
        color : str
          the color of the wireframe, in chemlab colors
          
        """
        vectors = self._unit_conversion(vectors.copy(), 'distance')
        origin = self._unit_conversion(origin.copy(), 'distance')
        color = str_to_colour(color)
        
        self._boxes.append([vectors, origin, color, width])

    def remove_boxes(self, n=1):
        """ remove the last n boxes to be added """
        assert len(self._boxes) >= n
        self._boxes = self._boxes[:-n]

    def add_hexagon(self, vectors, origin=np.zeros(3), color='black', width=1):
        """ add wireframed hexagonal prism to visualisation
        
       vectors : np.ndarray((2,3), dtype=float)
          The three vectors representing the orthogonal, a,c lattices.
       origin : np.ndarray((3,3), dtype=float), default to zero
          The origin of the box.
        color : str
          the color of the wireframe, in chemlab colors
          
        """
        vectors = self._unit_conversion(vectors.copy(), 'distance')
        origin = self._unit_conversion(origin.copy(), 'distance')
        color = str_to_colour(color)
        
        self._hexagons.append([vectors, origin, color, width])

    def remove_hexagons(self, n=1):
        """ remove the last n boxes to be added """
        assert len(self._hexagons) >= n
        self._hexagons = self._hexagons[:-n]

    def add_axes(self, axes=np.array([[1,0,0],[0,1,0],[0,0,1]]), 
                 length=1., offset=(-1.2,0.2), colors=('red','green','blue'),
                 width=1.5):
        """ add axes 

        axes : np.array(3,3)
            to turn off axes, set to None
        axes_offset : tuple
            x, y offset from top top-left atom
        
        """
        self._axes = [length*axes/np.linalg.norm(axes, axis=1), width,
                      offset, [str_to_colour(col) for col in colors]]

    def add_plane(self, vectors, origin=np.zeros(3), rev_normal=False, 
                  color='red', alpha=1.):
        """ add flat plane to visualisation
        
       vectors : np.ndarray((2,3), dtype=float)
          The three vectors representing the edges of the plane.
       origin : np.ndarray((3,3), dtype=float), default to zero
          The origin of the plane.
       rev_normal : bool
           whether to reverse direction of normal (for lighting calculations)
       color : str
          the color of the plane, in chemlab colors
          
        """
        vectors = self._unit_conversion(vectors.copy(), 'distance')
        origin = self._unit_conversion(origin.copy(), 'distance')

        c = str_to_colour(color)
        colors = np.array([c,c,c,c,c,c])
        assert alpha <= 1. and alpha > 0., 'alpha must be between 0 and 1'
        
        n = np.cross(vectors[0], vectors[1])
        if rev_normal:
            n = -n
        n = n / np.linalg.norm(n)
        normals = np.array([n,n,n,n,n,n])
        
        vertices = np.array([origin,vectors[0]+origin,vectors[1]+origin,
                             vectors[0]+origin,vectors[0]+vectors[1]+origin,vectors[1]+origin])
        
        self._triangles.append([vertices, normals, colors, alpha])
    
    def remove_planes(self, n=1):
        """ remove the last n planes to be added """
        assert len(self._triangles) >= n
        self._triangles = self._triangles[:-n]

    def _unit_conversion(self, values, measure):
        """ 
        values : np.array 
        measure : str       
        """
        if not self._unit_dict.has_key(self._units):
            raise NotImplementedError
        if not self._unit_dict[self._units].has_key(measure):
            raise NotImplementedError

        return values * self._unit_dict[self._units][measure]

    def _trim_image(self, im):
        """
        a simple solution to trim whitespace on the image
        
        1. It gets the border colour from the top left pixel, using getpixel, 
        so you don't need to pass the colour.
        2. Subtracts a scalar from the differenced image, 
        this is a quick way of saturating all values under 100, 100, 100 to zero. 
        So is a neat way to remove any 'wobble' resulting from compression.
        """
        bg = Image.new(im.mode, im.size, im.getpixel((0,0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        if bbox:
            return im.crop(bbox)

    def _concat_images_horizontal(self, images, gap=10, background='white'):
        """ concatentate one or more PIL images horizontally 

        Parameters
        ----------
        images : PIL.Image list
            the images to concatenate
        gap : int
            the pixel gap between images
        background : PIL.ImageColor
            background color (as supported by PIL.ImageColor)
        """
        if len(images) == 1: return images[0]
        
        total_width = sum([img.size[0] for img in images]) + len(images)*gap
        max_height = max([img.size[1] for img in images])
        
        final_img = Image.new("RGBA", (total_width, max_height), color=background)
        
        horizontal_position = 0
        for img in images:
            final_img.paste(img, (horizontal_position, 0))
            horizontal_position += img.size[0] + gap
        
        return final_img

    def _concat_images_vertical(self, images, gap=10, background='white'):
        """ concatentate one or more PIL images vertically 

        Parameters
        ----------
        images : PIL.Image list
            the images to concatenate
        gap : int
            the pixel gap between images
        """
        if len(images) == 1: return images[0]
        
        total_width = max([img.size[0] for img in images]) 
        max_height = sum([img.size[1] for img in images]) + len(images)*gap
        
        final_img = Image.new("RGBA", (total_width, max_height), color=background)
        
        vertical_position = 0
        for img in images:
            final_img.paste(img, (0, vertical_position))
            vertical_position += img.size[1] + gap
        
        return final_img

    def get_image(self, xrot=0, yrot=0, zrot=0, fov=10., width=400, height=400):
        """ get image of atom configuration
        
        requires atoms to have, at least variables xs, yx, zs and type
        
        Parameters
        ----------
        rotx: rotation about x (degrees)
        roty: rotation about y (degrees)
        rotz: rotation about z (degrees)
        (start x-axis horizontal, y-axis vertical)
        
        Return
        ------
        image : PIL.Image

        """        
        # an array of all points in the image (used to calculate axes position)
        all_array = None

        # initialize graphic engine
        v = QtViewer()
        w = v.widget
        #TODO could add option to change background color, but then need inversion of other colors
        w.background_color = str_to_colour('white')
        w.initializeGL()
        w.camera.fov = fov

        ## ADD RENDERERS
        ## ----------------------------
         
        #atoms renderer
        for r_array, type_array, backend, alpha, shading in self._atoms:
            if alpha < 1.:
                transparent = True
                cols=np.array(self._atomcolors.values())
                cols[:,3] = int(alpha*255)
                colormap = dict(zip(self._atomcolors.keys(), cols.tolist()))
            else:
                colormap = self._atomcolors
                transparent = False
            v.add_renderer(AtomRenderer, r_array, type_array,
                        color_scheme=colormap, radii_map=self._atomradii,
                        backend=backend, transparent=transparent, shading=shading) 
            all_array = r_array if all_array is None else np.concatenate([all_array,r_array])
        
        #boxes render
        for box_vectors, box_origin, box_color, box_width in self._boxes:
            v.add_renderer(BoxRenderer,box_vectors,box_origin,
                           color=box_color, width=box_width)           
            #TODO account for other corners of box?
            b_array = box_vectors + box_origin                           
            all_array = b_array if all_array is None else np.concatenate([all_array,b_array])

        #hexagonal prism render
        for hex_vectors, hex_origin, hex_color, hex_width in self._hexagons:
            v.add_renderer(HexagonRenderer,hex_vectors,hex_origin,
                           color=hex_color, width=hex_width)           
            #TODO account for other vertices of hexagon?
            h_array = hex_vectors + hex_origin                           
            all_array = h_array if all_array is None else np.concatenate([all_array,h_array])

        #surfaces render
        for vertices, normals, colors, alpha in self._triangles:
            if alpha < 1.:
                transparent = True
                colors[:,3] = int(alpha*255)
            else:
                transparent = False
            v.add_renderer(TriangleRenderer,vertices, normals, colors, 
                           transparent=transparent)           
            all_array = vertices if all_array is None else np.concatenate([all_array,vertices])

        ## ----------------------------

        if all_array is None:
            # Cleanup
            w.close()
            v.clear()
            raise Exception('nothing available to render')

        # transfrom coordinate system
        w.camera.orbit_x(xrot*np.pi/180.)
        w.camera.orbit_y(yrot*np.pi/180.)
        w.camera.orbit_z(zrot*np.pi/180.)
        
        # axes renderer
        # TODO option to show a,b,c instead of x,y,z
        if self._axes is not None:
            axes, axes_width, axes_offset, axes_colors = self._axes          
            
            # find top-left coordinate after transformations and 
            # convert to original coordinate system
            
            ones = np.ones((all_array.shape[0],1))
            trans_array = np.apply_along_axis(w.camera.matrix.dot,1,np.concatenate((all_array,ones),axis=1))[:,[0,1,2]]
            t_top_left = [trans_array[:,0].min() + axes_offset[0], 
                          trans_array[:,1].max() + axes_offset[1], 
                          trans_array[:,2].min(), 1]

            x0, y0, z0 = np.linalg.inv(w.camera.matrix).dot(t_top_left)[0:3]
            origin = [x0, y0, z0]
           
            all_array = np.concatenate([all_array, [origin]])
            
            vectors = axes + origin
            for vector, color in zip(vectors, axes_colors):
                # for some reason it won't render if theres not a 'dummy' 2nd element
                startends = [[origin, vector],[origin, vector]]                      
                colors = [[color, color],[color, color]]
                #TODO add as arrows instead of lines 
                v.add_renderer(LineRenderer, startends, colors, width=axes_width)
                #TODO add x,y,z labels (look at chemlab.graphics.__init__?)

        w.camera.autozoom(all_array)

        # convert scene to image
        image = w.toimage(width, height)
        image = self._trim_image(image)

        # Cleanup
        w.close()
        v.clear()
       
        return image
        
    def visualise(self, images, columns=1): 
        """ visualise image(s) in IPython
        
        Parameters
        ----------
        images : list or single PIL.Image
        columns : int
            number of image columns 
        
        Returns
        -------
        image : IPython.display.Image

        """        
        try:
            img_iter = iter(images)
        except TypeError:
            img_iter = iter([images])
        
        img_columns = []
        img_rows = []
        for image in img_iter: 
            if len(img_columns) < columns:
                img_columns.append(image)
            else:
                img_rows.append(self._concat_images_horizontal(img_columns))
                img_columns = [image]
        if img_columns:
            img_rows.append(self._concat_images_horizontal(img_columns))
        image = self._concat_images_vertical(img_rows)
        
        b = BytesIO()
        image.save(b, format='png')
        data = b.getvalue()
        del b
        
        return ipy_Image(data=data)

    def basic_vis(self, atoms_df=None, sim_box=None, type_map={}, 
                  spheres=True, alpha=1, 
                  xrot=0, yrot=0, zrot=0, fov=10.,
                  axes=np.array([[1,0,0],[0,1,0],[0,0,1]]), axes_length=1., axes_offset=(-1.2,0.2),
                  width=400, height=400):
        """ basic visualisation shortcut
        
        invoking add_atoms, add_box (if sim_box), add_axes, get_image and visualise functions
        
        """
        if not atoms_df is None:
            self.add_atoms(atoms_df, type_map=type_map, spheres=spheres, alpha=alpha)
        if not sim_box is None:
            self.add_box(sim_box[0], sim_box[1])
        if not axes is None:
            self.add_axes(axes=axes, length=axes_length, offset=axes_offset)
            
        image = self.get_image(xrot=xrot, yrot=yrot, zrot=zrot, fov=fov, 
                               width=width, height=height)
        
        # cleanup
        self._atoms.pop()
        if not sim_box is None: self._boxes.pop()

        return self.visualise(image)