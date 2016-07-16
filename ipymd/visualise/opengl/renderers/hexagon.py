# -*- coding: utf-8 -*-
"""
Created on Mon May 16 12:41:12 2016

@author: cjs14
"""

from __future__ import division
import pkgutil
import math
import numpy as np
from OpenGL.GL import *

from .base import ShaderBaseRenderer
from ....shared.colors import black

class HexagonRenderer(ShaderBaseRenderer):
    '''Used to render one wireframed hexagonal prism.
    
       **Parameters**
        
       widget:
          The parent QChemlabWidget
       vectors: np.ndarray((2,3), dtype=float)
          The two vectors representing the orthogonal a,c crystal vectors.
       origin: np.ndarray((3,), dtype=float), default to zero
          The origin of the box.
       color: 4 int tuple
          r,g,b,a color in the range [0,255]
       width: float
          width of wireframe lines

    '''

    def __init__(self, widget, vectors, origin=np.zeros(3), color=black, width=1.5):
        vert = pkgutil.get_data("ipymd.visualise.opengl.renderers.opengl_shaders",
                                "default_persp.vert")
        frag = pkgutil.get_data("ipymd.visualise.opengl.renderers.opengl_shaders",
                                "no_light.frag")

        self.color = color
        super(HexagonRenderer, self).__init__(widget, vert, frag)
        self.origin = origin
        self.vectors = vectors
        self.width = width
        

    def _rotate(self, v, axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta degrees.
        """
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
        return np.dot(rotation_matrix, v)
        
    def draw_vertices(self):
        # We need 12 vertices to draw a hexagonal prism, top and bottom
        a, c = self.vectors

        b1 = a
        b2 = self._rotate(a, c, 60)
        b3 = self._rotate(a, c, 120)
        b4 = self._rotate(a, c, 180)
        b5 = self._rotate(a, c, 240)
        b6 = self._rotate(a, c, 300)
                
        t1 = b1 + c
        t2 = b2 + c
        t3 = b3 + c
        t4 = b4 + c
        t5 = b5 + c
        t6 = b6 + c

        lines_vertices = np.array([
            b1, b2, b2, b3, b3, b4, b4, b5, b5, b6, b6, b1,
            t1, t2, t2, t3, t3, t4, t4, t5, t5, t6, t6, t1,
            b1, t1,
            b2, t2,
            b3, t3,
            b4, t4,
            b5, t5,
            b6, t6,
        ]) 
        
        lines_vertices += self.origin
        r, g, b, a = self.color
        
        glColor4f(r/255, g/255, b/255, a/255)
        glLineWidth(self.width)
        glBegin(GL_LINES)
        for i in lines_vertices:
            glVertex3f(*i)
        glEnd()
        glLineWidth(1)
    
    def update(self, vectors):
        """Update the box vectors.
        """

        self.vectors = vectors
