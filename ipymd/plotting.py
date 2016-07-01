# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:45:06 2016

@author: cjs14
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from io import BytesIO
from IPython import get_ipython, display

class Plotting(object):
    """ a class to deal with data plotting """
    def __init__(self,nrows=1,ncols=1,figsize=(5,4)):
        """ a class to deal with data plotting 
        """
        #ensure IPython shows matplotlib in inline mode
        ipython = get_ipython()
        ipython.run_line_magic('matplotlib', 'inline')

        if nrows==0 and ncols==0:
            self._fig = plt.figure(figsize=figsize)
            self._axes = []
        else:
            self._fig, axes = plt.subplots(nrows,ncols,squeeze=False,figsize=figsize)
            self._axes = axes.flatten()
        plt.close()
    
    def _get_mplfigure(self):
        return self._fig

    def _set_mplfigure(self, fig):
        self._fig = fig
        self._axes = fig.get_axes()
            
    figure = property(_get_mplfigure, _set_mplfigure)
    
    def _get_axes(self,):
        return self._axes
    
    axes = property(_get_axes)
    
    def display_plot(self, tight_layout=False):
        """ display plot in IPython 

        if tight_layout is True it may crop anything outside axes
        
        """
        ipython = get_ipython()
        current_config = ipython.run_line_magic('config', "InlineBackend.print_figure_kwargs")
        new_config = current_config.copy()        
        if tight_layout:
            new_config['bbox_inches'] = 'tight'
        else:
            new_config['bbox_inches'] = None

        ipython.run_line_magic('config', 
            'InlineBackend.print_figure_kwargs = {0}'.format(new_config))

        display.display(self._fig)

        ipython.run_line_magic('config', 
            'InlineBackend.print_figure_kwargs = {0}'.format(current_config))        

    def get_image(self,size=300,dpi=300, tight_layout=False):
        """return as PIL image

        if tight_layout is True it may crop anything outside axes
        
        """
        bbox_inches = 'tight' if tight_layout else None
        buf = BytesIO()
        self._fig.savefig(buf, dpi=dpi,format='png',bbox_inches=bbox_inches)
        buf.seek(0)
        img = Image.open(buf)
        img.thumbnail((int(size),int(size)),Image.ANTIALIAS)
        buf.close()
        return img
    
    def resize_axes(self,width=0.8,height=0.8,left=0.1,bottom=0.1, axes=0):
        """ resiaze axes, for instance to fit object outside of it """
        self._axes[axes].set_position([left,bottom,width,height])

    def add_image(self, image, axes=0, interpolation="bicubic", no_axis=True):

        self._axes[axes].imshow(image, interpolation="bicubic")
        
        if no_axis:
            self._axes[axes].get_xaxis().set_visible(False)
            self._axes[axes].get_yaxis().set_visible(False)
            self._axes[axes].set_frame_on(False)
    
    # from from http://matplotlib.org/examples/pylab_examples/demo_annotation_box.html
    # TODO resize axes automatically
    def add_image_annotation(self, img, xy=(0,0), arrow_xy=None, axes=0,
                  zoom=1, xytype='axes points', arrow_xytype='data',
                  arrowprops=dict(facecolor='black',
                                  arrowstyle="simple",
                                  connectionstyle="arc3,rad=0.2",
                                  alpha=0.4)):
        """ add an image to the plot 
        
        coordtype:
        
        ====================  ====================================================
        argument              coordinate system
        ====================  ====================================================
          'figure points'     points from the lower left corner of the figure
          'figure pixels'     pixels from the lower left corner of the figure
          'figure fraction'   0,0 is lower left of figure and 1,1 is upper right
          'axes points'       points from lower left corner of axes
          'axes pixels'       pixels from lower left corner of axes
          'axes fraction'     0,0 is lower left of axes and 1,1 is upper right
          'data'              use the axes data coordinate system
        ====================  ====================================================
        
        
        for arrowprops see http://matplotlib.org/users/annotations_guide.html#annotating-with-arrow
        """
        imagebox = OffsetImage(img, zoom=zoom)
        
        if arrow_xy is None:
            arrow_xy = (0,0)
            arrow_xytype='data'
            arrowprops={}
    
        ab = AnnotationBbox(imagebox, xy=arrow_xy,
                            xybox=xy,
                            xycoords=arrow_xytype,
                            boxcoords=xytype,#"offset points",
                            pad=0.5,
                            arrowprops=arrowprops,
                            )    
        self._axes[axes].add_artist(ab)

#TODO convert colors to r,g,b before creating set (in case mixed strings and rgb)
def create_legend_image(labels, colors, size=100, ncol=1, title=None, frameon=False,
                        sort_labels=True, colbytes=False, dpi=300, **kwargs):
    """ create a standalone image of a legend 

    labels : list
        a list of text labels
    colors : list
        a list of colors, as (r,g,b) or names from matplotlib
    colbytes : bool
        whether colors are in byte format (1-255) or not (0-1)
    kwargs
        additional arguments for legend matplotlib.legend
    """
    patches = []
    names=[]
    for label,color in set(zip(labels,colors)):
        if not isinstance(color,basestring) and colbytes:
            color = [i/255. for i in color]
        
        names.append(label)
        patches.append(mpl.patches.Patch(color=color))
    
    if sort_labels:    
        patches = [patch for (name,patch) in sorted(zip(names,patches))]
        names = sorted(names)
                
    plot = Plotting(0,0,figsize=(0.1,0.1))
    plot.figure.legend(patches,names,loc='center',frameon=frameon,ncol=ncol,title=title, **kwargs)
        
    return plot.get_image(size,tight_layout=True,dpi=dpi)   

#http://matplotlib.org/examples/api/colorbar_only.html
def create_colormap_image(vmin,vmax,title='',cmap='jet',ticks=2,length=2, width=0.4,
                          horizontal=True, size=100, dpi=300):

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if horizontal:
        plot = Plotting(figsize=(length,width))
    else:
        plot = Plotting(figsize=(width,length))
    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    orientation = 'horizontal' if horizontal else 'vertical'
    cb1 = mpl.colorbar.ColorbarBase(plot.axes[0], cmap=cmap,
                                    norm=norm,
                                    orientation=orientation)
    cb1.set_ticks(np.linspace(vmin,vmax,num=ticks))
    if title:
        cb1.set_label(title)
    return plot.get_image(size,tight_layout=True,dpi=dpi)     
