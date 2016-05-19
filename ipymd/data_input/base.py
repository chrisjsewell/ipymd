# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:49:20 2016

@author: cjs14
"""
import itertools

_col_dict = {
'blues': ['light_steel_blue', 'powder_blue', 'light_blue', 'sky_blue', 
    'light_sky_blue', 'deep_sky_blue', 'dodger_blue', 'cornflower_blue', 'steel_blue', 
    'royal_blue', 'blue', 'medium_blue', 'dark_blue', 'navy', 'midnight_blue'], 
'oranges': ['orange_red', 'tomato', 'coral', 'dark_orange', 'orange', 'gold'], 
'browns': ['cornsilk', 'blanched_almond', 'bisque', 'navajo_white', 'wheat', 
    'burly_wood', 'tan', 'rosy_brown', 'sandy_brown', 'goldenrod', 'dark_goldenrod', 
    'peru', 'chocolate', 'saddle_brown', 'sienna', 'brown', 'maroon'], 
'greens': ['dark_olive_green', 'olive', 
    'olive_drab', 'yellow_green', 'lime_green', 'lime', 'lawn_green', 'chartreuse', 
    'green_yellow', 'spring_green', 'medium_spring_green', 'light_green', 'pale_green', 
    'dark_sea_green', 'medium_sea_green', 'sea_green', 'forest_green', 'green', 
    'dark_green'], 
'purples': ['lavender', 'thistle', 'plum', 'violet', 'orchid', 'fuchsia', 'magenta', 
    'medium_orchid', 'medium_purple', 'blue_violet', 'dark_violet', 'dark_orchid', 
    'dark_magenta', 'purple', 'indigo', 'dark_slate_blue', 'slate_blue', 
    'medium_slate_blue'], 
'yellows': ['yellow', 'light_yellow', 'lemon_chiffon', 
    'light_goldenrod_yellow', 'papaya_whip', 'moccasin', 'peach_puff', 'pale_goldenrod', 
    'khaki', 'dark_khaki'], 
'pinks': ['pink', 'light_pink', 'hot_pink', 'deep_pink', 
    'pale_violet_red', 'medium_violet_red'], 
'cyans': ['medium_aquamarine', 'aqua', 'cyan', 'light_cyan', 'pale_turquoise', 
    'aquamarine', 'turquoise', 'medium_turquoise', 'dark_turquoise', 'light_sea_green', 
    'cadet_blue', 'dark_cyan', 'teal'], 
'reds': ['light_salmon', 'salmon', 'dark_salmon', 
         'light_coral', 'indian_red', 'crimson', 'fire_brick', 'dark_red', 'red']}

class DataInput(object):
    """
    Data is divided into two levels; sytem and atom
    
    system is simply a table containing variables (columns) for each timestep (rows)
    atom is a series of tables, one for each timestep, 
        containing variables (columns) for each atom (rows)
    systems bounding box is also supplied for each time step
    """
    def __init__(self):
        raise NotImplemented
    def get_system_data(self, step=None):
        """ return pandas.DataFrame or pandas.Series (for single step) """
        raise NotImplemented
    def get_atom_data(self, step):
        """ return pandas.DataFrame """
        raise NotImplemented
    def get_simulation_box(self, step):
        """ return list of coordinates (np.array(3)) for [a,b,c] & origin"""
        raise NotImplemented
        
    def count_timesteps(self):
        """ return int of total number of timesteps """
        raise NotImplemented

    def _add_radii(self, atom_df):
        atom_df['radius'] = 1.
        
    #TODO cycle through specific colors in each list
    def _add_colors(self, atom_df):
        """ add colors to atom_df, with different color for each atom type """

        atom_df['transparency'] = 1.
        
        col_keys = _col_dict.keys()
        col_cycle = itertools.cycle(col_keys)
        for typ in atom_df['type'].unique():            
            atom_df.loc[atom_df['type']==typ,'color'] = _col_dict[col_cycle.next()][0]
        
    def _skiplines(self, f, num=1):
        """ skip line(s) in an open file """
        for n in range(num):
            line = next(f)
        return line