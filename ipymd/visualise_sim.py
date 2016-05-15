# -*- coding: utf-8 -*-
"""
Created on Sun May  1 23:47:03 2016

@author: cjs14
"""
from io import BytesIO

import numpy as np

from chemlab.graphics.qtviewer import QtViewer
from chemlab.graphics.renderers.atom import AtomRenderer
from chemlab.graphics.renderers.box import BoxRenderer
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

    def get_image(self, atoms_df, sim_box=None, type_dict={}, spheres=True, 
                  xrot=0, yrot=0, zrot=0, fov=10., 
                  show_axes=True, axes_offset=(-0.2,0.2), axes_length=1,
                  width=400, height=400):
        """ get image of atom configuration
        
        requires atoms to have, at least variables xs, yx, zs and type
        
        Parameters
        ----------
        atom_df : pandas.DataFrame
            a dataframe of variables for each atom
        sim_bounds : np.array((6,3), dtype=float)
            list of coordinates for simulation bounds [x0,y0,z0,a,b,c]
        sphere : render spheres, otherwise points
        rotx: rotation about x (degrees)
        roty: rotation about y (degrees)
        rotz: rotation about z (degrees)
        (start x-axis horizontal, y-axis vertical)
        
        Return
        ------
        image : PIL.Image

        """
        assert set(['xs','ys','zs','type']).issubset(set(atoms_df.columns))
        r_array = np.array([s[['xs','ys','zs']] for i,s in atoms_df.iterrows()])
        r_array = self._unit_conversion(r_array, 'distance')
        all_array = r_array
        
        type_array = [type_dict.get(s['type'], s['type']) for i,s in atoms_df.iterrows()]

        # initialize graphic engine
        v = QtViewer()
        w = v.widget
        #TODO could add option to change background, but should also have simulation box change
        w.background_color = str_to_colour('white')
        w.initializeGL()
        w.camera.fov = fov

        ## add renderers
        ## ---------------
        rends = []

        #simulation bounding box render
        if not sim_box is None:
            sim_box = self._unit_conversion(sim_box, 'distance')

            # move r_array so origin is at (0,0,0)
            r_array[:,0] = r_array[:,0] - sim_box[0,0]           
            r_array[:,1] = r_array[:,1] - sim_box[1,1]            
            r_array[:,2] = r_array[:,2] - sim_box[2,2]           
            
            vectors = sim_box[3:6]
            rends.append(v.add_renderer(BoxRenderer, vectors))
            
            #TODO account for other corners of box?
            all_array = np.concatenate([r_array,vectors])
         
        #atoms renderer
        backend = 'impostors' if spheres else 'points'
        rends.append(v.add_renderer(AtomRenderer, r_array, type_array,
                    color_scheme=self._atomcolors, radii_map=self._atomradii,
                    backend=backend))            
        
        # transfrom coordinate system
        w.camera.orbit_x(xrot*np.pi/180.)
        w.camera.orbit_y(yrot*np.pi/180.)
        w.camera.orbit_z(zrot*np.pi/180.)
        
        # axes renderer
        # TODO option to show a,b,c instead of x,y,z
        if show_axes:
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
            
            vectors = [[x0+axes_length,y0,z0], 
                       [x0,y0+axes_length,z0], 
                       [x0,y0,z0+axes_length]]
            # colors consistent with ovito
            colors = [str_to_colour('red'),str_to_colour('green'),str_to_colour('blue')]
            for vector, color in zip(vectors, colors):
                # for some reason it won't render if theres not a 'dummy' 2nd element
                startends = [[origin, vector],[origin, vector]]                      
                colors = [[color, color],[color, color]]
                #TODO add as arrows instead of lines 
                rends.append(v.add_renderer(LineRenderer, startends, colors))
                #TODO add x,y,z labels (look at chemlab.graphics.__init__?)

        w.camera.autozoom(all_array)

        # convert scene to image
        image = w.toimage(width, height)
        image = self._trim_image(image)

        # Cleanup
        for r in rends:
            del r
        del v
        del w
        
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
        return ipy_Image(data=data)
