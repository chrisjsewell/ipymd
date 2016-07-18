# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 16:45:06 2016

@author: cjs14
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from PIL import Image
from io import BytesIO
from IPython import get_ipython, display

from IPython.display import HTML
from .JSAnimation.IPython_display import display_animation

def style(style):
    """A context manager to apply matplotlib style settings from a style specification.
    
    Popular styles include; default, ggplot, xkcd, and are used in the the following manner:
    
    ::
    
        with ipymd.plotting.style('default'):
            plot = ipymd.plotting.Plotter()
            plot.display_plot()

    Parameters
    ----------
    style : str, dict, or list
        A style specification. Valid options are:

        +------+-------------------------------------------------------------+
        | str  | The name of a style or a path/URL to a style file. For a    |
        |      | list of available style names, see `style.available`.       |
        +------+-------------------------------------------------------------+
        | dict | Dictionary with valid key/value pairs for                   |
        |      | `matplotlib.rcParams`.                                      |
        +------+-------------------------------------------------------------+
        | list | A list of style specifiers (str or dict) applied from first |
        |      | to last in the list.                                        |
        +------+-------------------------------------------------------------+

    """
    if style=='xkcd':
        return plt.xkcd()
    else: 
        return plt.style.context(style)    
    
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
    
    def _get_axes(self):
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
            self._fig.show()
        else:        
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
        
    def add_image(self, image, axes=0, interpolation="bicubic", hide_axes=True, 
                  width=1., height=1.,origin=(0.,0.), **kwargs):
        """add image to axes
        """        
        x0,y0=origin
        self._axes[axes].imshow(image, interpolation="bicubic", extent=(x0,x0+width,y0,y0+height) ,**kwargs)
        
        if hide_axes:
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
    
def animation_line(x_iter, y_iter, interval=20, xlim=(0,1),ylim=(0,1),
                 incl_controls=True, plot=None,ax=0,**plot_kwargs):
    """create an animation of multiple x,y data sets
    
    x_iter : iterable
        any iterable of x data sets, e.g. [[1,2,3],[4,5,6]] 
    y_iter : iterable
        an iterable of y data sets, e.g. [[1,2,3],[4,5,6]]
    interval : int
        draws a new frame every *interval* milliseconds  
    xlim : tuple
        the x_limits for the axes (ignored if using existing plotter)
    ylim : tuple
        the y_limits for the axes (ignored if using existing plotter)
    incl_controls : bool
        include Javascript play controls
    plot : ipymd.plotting.Plotter
        an existing plotter object
    ax : int
        the id number of the axes on which to plot (if using existing plotter) 
    plot_kwargs : various
        key word arguments to pass to plot method, e.g. marker='o', color='b', ...
    
    Returns
    -------
    html : IPython.core.display.HTML
        a html object
        
    Notes
    -----
    x_iter and y_iter can be yield functions such as:
    
    ::
    
        def y_iter(x_iter):
            for xs in x_iter:
                yield [i**2 for i in xs]
    
    This means that the values do not have to be necessarily pre-computed.
        
    """
    if plot is None:
        plotter = Plotter()
    else:
        plotter = plot
    if isinstance(plotter.axes, list):
        ax = plotter.axes[ax]
    else:
        ax = plotter.axes
    if plot is None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    xiter = iter(x_iter)
    yiter = iter(y_iter)
    
    xy_plot, = ax.plot([], [], animated=True, **plot_kwargs)
    
    # initialization function: plot the background of each frame
    def init():
        xy_plot.set_data([], [])
        return (xy_plot,)
    
    # animation function. This is called sequentially
    def animate(i):
        try:
            xs = xiter.next()
            ys = yiter.next()
        except StopIteration:
            xs = []
            ys = []
        xy_plot.set_data(xs,ys)
        return (xy_plot,)            

    anim = animation.FuncAnimation(plotter.figure, animate, init_func=init,
                                   frames=x_iter, interval=interval, blit=True)
    
    if incl_controls:
        html =  display_animation(anim)
    else:
        html = HTML(anim.to_html5_video())
    
    # cleanup
    xy_plot.remove()
    
    return html

def animation_scatter(x_iter, y_iter, interval=20, xlim=(0,1),ylim=(0,1),
                 incl_controls=True, plot=None,ax=0,**plot_kwargs):
    """create an animation of multiple x,y data sets
    
    x_iter : iterable
        any iterable of x data sets, e.g. [[1,2,3],[4,5,6]] 
    y_iter : iterable
        an iterable of y data sets, e.g. [[1,2,3],[4,5,6]]
    interval : int
        draws a new frame every *interval* milliseconds  
    xlim : tuple
        the x_limits for the axes (ignored if using existing plotter)
    ylim : tuple
        the y_limits for the axes (ignored if using existing plotter)
    incl_controls : bool
        include Javascript play controls
    plot : ipymd.plotting.Plotter
        an existing plotter object
    ax : int
        the id number of the axes on which to plot (if using existing plotter) 
    plot_kwargs : various
        key word arguments to pass to plot method, e.g. marker='o', color='b', ...
    
    Returns
    -------
    html : IPython.core.display.HTML
        a html object
        
    Notes
    -----
    x_iter and y_iter can be yield functions such as:
    
    ::
    
        def y_iter(x_iter):
            for xs in x_iter:
                yield [i**2 for i in xs]
    
    This means that the values do not have to be necessarily pre-computed.
        
    """
    if plot is None:
        plotter = Plotter()
    else:
        plotter = plot
    if isinstance(plotter.axes, list):
        ax = plotter.axes[ax]
    else:
        ax = plotter.axes
    if plot is None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    xiter = iter(x_iter)
    yiter = iter(y_iter)
    
    xy_plot = ax.scatter([], [], animated=True, **plot_kwargs)
    
    # initialization function: plot the background of each frame
    def init():
        xy_plot.set_offsets([])
        return (xy_plot,)
    
    # animation function. This is called sequentially
    def animate(i):
        try:
            xs = xiter.next()
            ys = yiter.next()
        except StopIteration:
            xs = []
            ys = []
        xy_plot.set_offsets(zip(xs,ys))
        return (xy_plot,)            

    anim = animation.FuncAnimation(plotter.figure, animate, init_func=init,
                                   frames=x_iter, interval=interval, blit=True)
    
    if incl_controls:
        html =  display_animation(anim)
    else:
        html = HTML(anim.to_html5_video())
    
    # cleanup
    xy_plot.remove()
    
    return html

def animation_contourf(x_iter, y_iter, z_iter, interval=20, 
                       xlim=(0,1),ylim=(0,1),zlim=(0,1.), 
                        cmap='viridis', cbar=True,
                        incl_controls=True, plot=None,ax=0,**plot_kwargs):
    """create an animation of multiple x,y data sets
    
    x_iter : iterable
        any iterable of x data sets, e.g. [[1,2,3],[4,5,6]] 
    y_iter : iterable
        an iterable of y data sets, e.g. [[1,2,3],[4,5,6]]
    y_iter : iterable
        an iterable of z(x,y) data sets, each set must be of shape (len(x), len(y))
    interval : int
        draws a new frame every *interval* milliseconds  
    xlim : tuple
        the x_limits for the axes (ignored if using existing plotter)
    ylim : tuple
        the y_limits for the axes (ignored if using existing plotter)
    zlim : tuple
        the z_limits for the colormap
    cmap : str or matplotlib.cm
        the colormap to use (see http://matplotlib.org/examples/color/colormaps_reference.html)
    incl_controls : bool
        include Javascript play controls
    plot : ipymd.plotting.Plotter
        an existing plotter object
    ax : int
        the id number of the axes on which to plot (if using existing plotter) 
    plot_kwargs : various
        key word arguments to pass to plot method, e.g. marker='o', color='b', ...
    
    Returns
    -------
    html : IPython.core.display.HTML
        a html object
        
    Notes
    -----
    x_iter and y_iter can be yield functions such as:
    
    ::
    
        def y_iter(x_iter):
            for xs in x_iter:
                yield [i**2 for i in xs]
    
    This means that the values do not have to be necessarily pre-computed.
        
    """
    if plot is None:
        plotter = Plotter()
    else:
        plotter = plot
    if isinstance(plotter.axes, list):
        ax = plotter.axes[ax]
    else:
        ax = plotter.axes
    if plot is None:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    
    xiter = iter(x_iter)
    yiter = iter(y_iter)
    ziter = iter(z_iter)    
    
    zmin, zmax = zlim
    
    # has to be in a list to pass between nested functions
    c_plot = [ax.contourf([0,1], [0,1], [[0,1],[0,1]],
                          vmin=zmin,vmax=zmax,cmap=cmap, 
                          animate=True,**plot_kwargs)]
    if cbar:
        plt.colorbar(c_plot[0], ax=ax)
        plt.close()
    
    
    # initialization function: plot the background of each frame
    def init():
        for coll in c_plot[0].collections:
            ax.collections.remove(coll)
        c_plot[0] = ax.contourf([0,1], [0,1], [[0,0],[0,0]],
                                vmin=zmin,vmax=zmax,cmap=cmap, 
                                animate=True,**plot_kwargs)        
        return c_plot[0].collections
    
    # animation function. This is called sequentially
    def animate(i):
        try:
            xs = xiter.next()
            ys = yiter.next()
            zs = ziter.next()
        except StopIteration:
            xs = []
            ys = []
            zs = []
            
        X, Y = np.meshgrid(xs, ys)
        
        for coll in c_plot[0].collections:
            ax.collections.remove(coll)
        
        c_plot[0] = ax.contourf(X,Y,zs,vmin=zmin,vmax=zmax,cmap=cmap,
                                animate=True,**plot_kwargs)        
        return c_plot[0].collections

    anim = animation.FuncAnimation(plotter.figure, animate, init_func=init,
                                   frames=x_iter, interval=interval, blit=True)
    
    if incl_controls:
        html =  display_animation(anim)
    else:
        html = HTML(anim.to_html5_video())
    
    # cleanup
    for coll in c_plot[0].collections:
        ax.collections.remove(coll)
    
    return html


