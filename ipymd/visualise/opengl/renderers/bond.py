from .base import AbstractRenderer
from .line import LineRenderer
from .cylinder_imp import CylinderImpostorRenderer
import numpy as np

class BondRenderer(AbstractRenderer):
    '''
        Render chemical bonds as cylinders or lines.

        **Parameters**
    
        widget:
           The parent QChemlabWidget
        r_array: np.ndarray((NATOMS, 3), dtype=float)
            The coordinate array
        bonds: np.ndarray((NBONDS, 2), dtype=int)
            An array of integer pairs that represent the bonds.
        colors_start: np.ndarray((NBONDS, 1), dtype=int)
            An array of colors that represent the bond color.
        colors_end: np.ndarray((NBONDS, 1), dtype=int)
            An array of colors that represent the bond color.
        radius: float, default=0.02
            The radius of the bonds
        style: "cylinders" | "lines"
            Whether to render the bonds as cylinders or lines.
    
    '''
    def __init__(self, widget, starts, ends, colors_start, colors_end, radii,
                 backend="impostors", shading='phong', transparent=True,
                 linewidth=5):
        #super(BondRenderer, self).__init__(widget)
        
        
        self.radii = np.asarray(radii)
        self.colors_start = np.array(colors_start, 'uint8')
        self.colors_end = np.array(colors_end, 'uint8')
        
        if backend == 'lines':
            # need to duplicate color for start and end of the line
            cols = np.empty((self.colors_start.shape[0],2,4),dtype='uint8')
            cols[:,0,:] = self.colors_start
            cols[:,1,:] = self.colors_end
            
            self.cr1 = LineRenderer(widget, np.array(zip(starts,ends)),
                                    cols,width=linewidth)
            self.cr2 = None
        elif backend == 'impostors':
            
            middles = (starts + ends)/2
            
            bounds = np.empty((len(starts), 2, 3))
            bounds[:, 0, :] = starts
            bounds[:, 1, :] = middles
            
            self.cr1 = CylinderImpostorRenderer(widget, bounds, self.radii,
                                                self.colors_start, shading=shading,
                                                transparent=transparent)

            bounds = np.empty((len(starts), 2, 3))
            bounds[:, 0, :] = middles
            bounds[:, 1, :] = ends

            self.cr2 = CylinderImpostorRenderer(widget, bounds, self.radii,
                                                self.colors_end, shading=shading,
                                                transparent=transparent)
        else:
            raise Exception("Available backends: lines, impostors")
                    
    def draw(self):
        self.cr1.draw()
        if self.cr2 is not None:
            self.cr2.draw()

    #TODO bond/cylinder update functions
    def update_positions(self, bounds):
        if bounds.size == 0:
            return
        
        self.cr1.update_bounds(bounds)

    def update_colors(self, colors_start,colors_end):
        self.colors_start = np.array(colors_start, 'uint8')
        self.colors_end = np.array(colors_end, 'uint8')
        
        self.cr1.update_colors(self.colors_start)
