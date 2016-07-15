# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 01:56:51 2016

@author: cjs14
"""

col_dict = {
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

def available_colors():
    return col_dict.copy()
