# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:45:06 2016

@author: cjs14
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from io import BytesIO
from IPython import get_ipython, display

from ._xkcdify import _XKCDify
    
class Plotter(object):
    """ a class to deal with data plotting """
    def __init__(self,nrows=1,ncols=1,figsize=(5,4)):
        """ a class to deal with data plotting 
        
        Attributes
        ----------
        figure : matplotlib.figure
            the figure
        axes : list or single matplotlib.axes
            if more than one then returns a list (ordered in reading direction),
            else returns one instance
            
        """
        #ensure IPython shows matplotlib in inline mode
        ipython = get_ipython()
        if ipython is not None:
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
        if len(self._axes) == 1:
            return self._axes[0]
        else:
            return self._axes
    
    axes = property(_get_axes)
    
    def display_plot(self, tight_layout=False):
        """ display plot in IPython 

        if tight_layout is True it may crop anything outside axes
        
        """
        ipython = get_ipython()
        if ipython is None:
            raise RuntimeError('not in an ipython shell')
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

    def apply_xkcd_style(self,axes=0,mag=1.0,
            f1=50, f2=0.01, f3=15,
            bgcolor='w',
            xaxis_loc=None, yaxis_loc=None,
            xaxis_arrow='+', yaxis_arrow='+',
            ax_extend=0.1):
        """ apply the xkcd style to one or more axes i.e. for schematic plots.
        This should be done after the axes is finalised (i.e everthing plotted)
        
        Parameters
        ----------
        axes : int or list of ints
            the axes to be modified.
        mag : float
            the magnitude of the distortion
        f1, f2, f3 : int, float, int
            filtering parameters.  f1 gives the size of the window, f2 gives
            the high-frequency cutoff, f3 gives the size of the filter
        xaxis_loc, yaxis_loc : float
            The locations to draw the x and y axes.  If not specified, they
            will be drawn from the bottom left of the plot
        xaxis_arrow, yaxis_arrow : str
            where to draw arrows on the x/y axes.  Options are '+', '-', '+-', or ''
        ax_extend : float
            How far (fractionally) to extend the drawn axes beyond the original
            axes limits
        """
        if len(self._axes)==1:
            expand_axes=True
        else:
            expand_axes=False
        for ax_id in np.array(axes, ndmin=1):
            self._axes[axes] = _XKCDify(self._axes[axes],mag,f1,f2,f3,bgcolor,
            xaxis_loc,yaxis_loc,xaxis_arrow,yaxis_arrow, ax_extend, expand_axes)
        
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

