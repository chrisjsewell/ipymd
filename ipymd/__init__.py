# -*- coding: utf-8 -*-
"""
Created on Sun May  1 22:46:22 2016

@author: cjs14
"""
# so tab completion works
from . import data_input
from . import visualise_sim
from . import atom_manipulation
from . import atom_analysis

import os
import inspect
from . import test_data
def get_test_path(data, check_exists=False):
    """return a directory path to the test data

    data : str or list of str
        file name or list of sub-directories and file name (e.g. ['lammps','data.txt'])   
    """
    basepath = os.path.dirname(os.path.abspath(inspect.getfile(test_data)))
    
    if isinstance(data, basestring): data = [data]
    
    dirpath = os.path.join(basepath, *data)
    
    if check_exists:
        assert os.path.exists(dirpath), '{0} does not exist'.format(dirpath)
    
    return dirpath
