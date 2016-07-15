# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 02:53:47 2016

@author: cjs14
"""
import numpy as np
import matplotlib as mpl
from chemlab.graphics.colors import get as str_to_colour

from plotter import Plotter

#TODO convert colors to r,g,b before creating set (in case mixed strings and rgb)
def create_legend_image(labels, colors, size=100, ncol=1, title=None, frameon=False,
                        sort_labels=True, colbytes=False, chemlabcols=True, dpi=300, **kwargs):
    """ create a standalone image of a legend 

    labels : list
        a list of text labels
    colors : list
        a list of colors, as (r,g,b) or names from matplotlib
    colbytes : bool
        whether colors are in byte format (1-255) or not (0-1)
    chemlabcols : bool
        using color names defined by chemlab, otherwise matplotlib   
    kwargs
        additional arguments for legend matplotlib.legend
    """
    patches = []
    names=[]
    for label,color in set(zip(labels,colors)):
        if not isinstance(color,basestring) and colbytes:
            color = [i/255. for i in color]
        if isinstance(color,basestring) and chemlabcols:
            color = str_to_colour(color)
            color = [i/255. for i in color]
        
        names.append(label)
        patches.append(mpl.patches.Patch(color=color))
    
    if sort_labels:    
        patches = [patch for (name,patch) in sorted(zip(names,patches))]
        names = sorted(names)
                
    plot = Plotter(0,0,figsize=(0.1,0.1))
    plot.figure.legend(patches,names,loc='center',frameon=frameon,ncol=ncol,title=title, **kwargs)
        
    return plot.get_image(size,tight_layout=True,dpi=dpi)   

#http://matplotlib.org/examples/api/colorbar_only.html
def create_colormap_image(vmin,vmax,title='',cmap='jet',ticks=2,length=2, width=0.4,
                          horizontal=True, size=100, dpi=300):

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    if horizontal:
        plot = Plotter(figsize=(length,width))
    else:
        plot = Plotter(figsize=(width,length))
    # ColorbarBase derives from ScalarMappable and puts a colorbar
    # in a specified axes, so it has everything needed for a
    # standalone colorbar.  There are many more kwargs, but the
    # following gives a basic continuous colorbar with ticks
    # and labels.
    orientation = 'horizontal' if horizontal else 'vertical'
    cb1 = mpl.colorbar.ColorbarBase(plot.axes, cmap=cmap,
                                    norm=norm,
                                    orientation=orientation)
    cb1.set_ticks(np.linspace(vmin,vmax,num=ticks))
    if title:
        cb1.set_label(title)
    return plot.get_image(size,tight_layout=True,dpi=dpi)     
