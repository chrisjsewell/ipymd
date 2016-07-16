import OpenGL.GL as ogl
from OpenGL.raw import GL
from OpenGL.arrays import ArrayDatatype as ADT

import numpy as np

class VertexBuffer(object):

  def __init__(self, data, usage):
    self.buffer = ogl.GLuint(0)
    self.buffer = ogl.glGenBuffers(1)
    self.usage = usage
    self.data = data
    
    # Add a little warning
    if data.dtype == np.float:
      Warning('This float array is 64 bit')
    
    ogl.glBindBuffer(ogl.GL_ARRAY_BUFFER, self.buffer)
    ogl.glBufferData(ogl.GL_ARRAY_BUFFER, ADT.arrayByteCount(data),
                 ADT.voidDataPointer(data), usage)
    ogl.glBindBuffer(ogl.GL_ARRAY_BUFFER, 0)

  def __del__(self):
    # this causes an error otherwise
    if bool(ogl.glDeleteBuffers):
      ogl.glDeleteBuffers(1, ogl.GLuint(self.buffer))
    else:
      return

  def bind(self):
    ogl.glBindBuffer(ogl.GL_ARRAY_BUFFER, self.buffer)
    
  def unbind(self):
    ogl.glBindBuffer(ogl.GL_ARRAY_BUFFER, 0)

  def set_data(self, data):
    self.bind()
    ogl.glBufferData(ogl.GL_ARRAY_BUFFER, ADT.arrayByteCount(data), ADT.voidDataPointer(data), self.usage)
    self.unbind()
      
  def bind_colors(self, size, type, stride=0):
    self.bind()
    ogl.glColorPointer(size, type, stride, None)

  def bind_edgeflags(self, stride=0):
    self.bind()
    ogl.glEdgeFlagPointer(stride, None)

  def bind_indexes(self, type, stride=0):
    self.bind()
    ogl.glIndexPointer(type, stride, None)

  def bind_normals(self, type, stride=0):
    self.bind()
    ogl.glNormalPointer(type, stride, None)

  def bind_texcoords(self, size, type, stride=0):
    self.bind()
    ogl.glTexCoordPointer(size, type, stride, None)

  def bind_vertexes(self, size, type, stride=0):
    self.bind()
    ogl.glVertexPointer(size, type, stride, None)
    
  def bind_attrib(self, attribute, size, type, normalized=ogl.GL_FALSE, stride=0):
    self.bind()
    ogl.glVertexAttribPointer(attribute, size, type, normalized, stride, None)