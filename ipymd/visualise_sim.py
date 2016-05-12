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

from IPython.display import Image as ipy_Image

class Visualise_Sim(object):
    """ 

    For style *real*, these are the units:
    
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
    _unit_dict = {'real':{'distance':0.1}}
    
    def __init__(self, units='real'):
        assert units=='real', 'currently only supports real units'
        self._units = units
        
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
        
    def visualise(self, atoms_df, type_dict={}, bounds=None, xrot=0, yrot=0, fov=10.,
                  width=400, height=400):
        """ 
        rotx: rotation about x 
        roty: rotation about y
        (start x-axis horizontal, y-axis vertical)
        bounds: np.array((3,2), dtype=float)
        """
        assert set(['xs','ys','zs','type']).issubset(set(atoms_df.columns))
        r_array = np.array([s[['xs','ys','zs']] for i,s in atoms_df.iterrows()])
        r_array = self._unit_conversion(r_array, 'distance')
        
        type_array = [type_dict.get(s['type'], 'Xx') for i,s in atoms_df.iterrows()]

        # initialize graphic engine
        v = QtViewer()
        w = v.widget 
        w.initializeGL()
        w.camera.fov = fov

        rends = []
        if not bounds is None:
            bounds = self._unit_conversion(bounds, 'distance')
            x0, y0, z0 = bounds[:,0]
            x1, y1, z1 = bounds[:,1]
            
            # move r_array so origin is at (0,0,0)
            r_array[:,0] = r_array[:,0] - x0           
            r_array[:,1] = r_array[:,1] - y0           
            r_array[:,2] = r_array[:,2] - z0           
            
            vectors = np.array([[x1-x0,0,0],[0,y1-y0,0],[0,0,z1-z0]])
            rends.append(v.add_renderer(BoxRenderer, vectors))
            
        rends.append(v.add_renderer(AtomRenderer, r_array, type_array))
        
        if not bounds is None:
            w.camera.autozoom(np.concatenate([r_array,vectors]))
        else:
            w.camera.autozoom(r_array)
        
        w.camera.orbit_x(xrot*np.pi/180.)
        w.camera.orbit_y(yrot*np.pi/180.)
        
        image = w.toimage(width, height)
        b = BytesIO()
        image.save(b, format='png')
        data = b.getvalue()

        # Cleanup
        for r in rends:
            del r
        del v
        del w
        
        return ipy_Image(data=data)