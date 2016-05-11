# -*- coding: utf-8 -*-
"""
Created on Sun May  1 23:47:03 2016

@author: cjs14
"""
from io import BytesIO

import numpy as np

from chemlab.graphics.qtviewer import QtViewer
from chemlab.graphics.renderers.atom import AtomRenderer

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
    _unit_dict = {'real':{'distance':10.}}
    
    def __init__(self, units='real'):
        assert units=='real', 'currently only supports real units'
        self._units = units
        
    def _unit_conversion(self, values, measure):
        
        if not self._unit_dict.has_key(self._units):
            raise NotImplementedError
        if not self._unit_dict[self._units].has_key(measure):
            raise NotImplementedError

        return values * self._unit_dict[self._units][measure]
        
    def visualise(self, atoms_df, type_dict={}, xrot=0, yrot=0, width=400, height=400):

        assert set(['xs','ys','zs','type']).issubset(set(atoms_df.columns))

        # initialize graphic engine
        v = QtViewer()
        w = v.widget 
        w.initializeGL()

        r_array = np.array([np.array([s['xs'], s['ys'], s['zs']]) for i,s in atoms_df.iterrows()])
        r_array = self._unit_conversion(r_array, 'distance')
        
        type_array = [type_dict.get(s['type'], 'Xx') for i,s in atoms_df.iterrows()]

        renderer = AtomRenderer
        r = v.add_renderer(renderer, r_array, type_array)

        w.camera.autozoom(r_array)
        w.camera.orbit_x(xrot*np.pi/180.)
        w.camera.orbit_y(yrot*np.pi/180.)
        
        image = w.toimage(width, height)
        b = BytesIO()
        image.save(b, format='png')
        data = b.getvalue()

        # Cleanup
        del r
        del v
        del w
        
        return ipy_Image(data=data)