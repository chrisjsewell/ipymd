# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:46:22 2016

@author: cjs14
"""
# so tab completion works
from . import data_input
from . import data_output
from .visualise import visualise_sim
from . import atom_manipulation
from . import atom_analysis
from . import plotting 
from . import shared
from .shared.colors import available_colors
from .shared import get_data_path 

from ._version import __version__

def version():
    return __version__

