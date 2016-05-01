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
    """ """
    def __init__(self):
        pass
    def visualise(self, atoms_df, width=400, height=400):
        # initialize graphic engine
        v = QtViewer()
        w = v.widget 
        w.initializeGL()

        r_array = np.array([np.array([s['xs'], s['ys'], s['zs']]) for i,s in atoms_df.iterrows()])
        type_array = ['Fe' for i,s in atoms_df.iterrows()]

        renderer = AtomRenderer
        v.add_renderer(renderer, r_array, type_array)

        w.camera.autozoom(r_array)
        w.camera.orbit_x(45*np.pi/180.)
        w.camera.orbit_y(45*np.pi/180.)
        
        image = w.toimage(width, height)
        b = BytesIO()
        image.save(b, format='png')
        data = b.getvalue()

        # Cleanup
        del v
        del w
        
        return ipy_Image(data=data)