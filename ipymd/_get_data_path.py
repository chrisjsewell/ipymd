# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 03:21:03 2016

@author: cjs14
"""
import os
import inspect
from . import test_data

def get_data_path(data, check_exists=False, module=test_data):
    """return a directory path to data within a module

    data : str or list of str
        file name or list of sub-directories and file name (e.g. ['lammps','data.txt'])   
    """
    basepath = os.path.dirname(os.path.abspath(inspect.getfile(module)))
    
    if isinstance(data, basestring): data = [data]
    
    dirpath = os.path.join(basepath, *data)
    
    if check_exists:
        assert os.path.exists(dirpath), '{0} does not exist'.format(dirpath)
    
    return dirpath
