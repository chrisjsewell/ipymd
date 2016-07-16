# shared resources

import os
import inspect
import pandas as pd

from .. import test_data
from . import atomdata

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
    
def atom_data():
    """return a dataframe of atomic data
    """
    path = get_data_path('element.txt',module=atomdata)
    df = pd.read_csv(path,comment='#')
    df.set_index('Symb',inplace=True)
    
    red = df.Red*255
    green = df.Green*255
    blue = df.Blue*255
    
    df['color'] = zip(red.values.astype(int),
                      green.values.astype(int),
                      blue.values.astype(int))
                      
    df.drop(['Red','Green','Blue'],axis=1,inplace=True)
    
    return df
