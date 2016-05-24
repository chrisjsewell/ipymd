# -*- coding: utf-8 -*-
"""
Created on Mon May 23 17:55:14 2016

@author: cjs14
"""
import os
import numpy as np
import datetime

from ._version import __version__

class Data_Output(object):
    """
    
    """
    def __init__(self, atom_df, abc, origin=np.zeros(3)):
        self._atom_df = atom_df.copy()
        self._abc = np.array(abc)
        self._origin = np.array(origin)
    
    def save_xyz(self, outpath='out.xyz', overwrite=False,
                 header=''):

        if os.path.exists(outpath) and not overwrite:
            raise IOError('file already exists; {0}'.format(outpath))
        
        raise NotImplementedError
        
    def save_gromacs(self, outpath='out.gro', overwrite=False,
                     header=''):

        if os.path.exists(outpath) and not overwrite:
            raise IOError('file already exists; {0}'.format(outpath))
        
        raise NotImplementedError

    def save_lammps(self, outpath='out.lammps', overwrite=False,
                    atom_type='atomic', header=''):
        """ to adhere to http://lammps.sandia.gov/doc/read_data.html?highlight=read_data 
        
        Example
        -------
        In [1]: import pandas as pd
        In [2]: df = pd.DataFrame([['Fe',2,3,4,1],
                                  ['Cr',2,3,3,-1],
                                  ['Fe',4,3,1,1]],columns=['type','xs','ys','zs','q'])
        In [3]: from ipymd.data_output import Data_Output as data_out
        In [4]: data = data_out(df, [[1,0,0],[0,1,0],[0,0,1]])
        In [5]: data.save_lammps('test.lammps', atom_type='charge', overwrite=True,
                                header='my header')
        In [6]: cat test.lammps
        # This file was created by ipymd (v0.0.1) on 2016-05-23 20:51:16 
        # type map: {'Cr': 2, 'Fe': 1} 
        # my header 
        
        3 atoms 
        2 atom types 
        
        # simulation box boundaries
        0.0000 1.0000 xlo xhi 
        0.0000 1.0000 ylo yhi 
        0.0000 1.0000 zlo zhi 
        0.0000 0.0000 0.0000 xy xz yz 
        
        Atoms 
        
        1 1 1.0000 2.0000 3.0000 4.0000 
        2 2 -1.0000 2.0000 3.0000 3.0000 
        3 1 1.0000 4.0000 3.0000 1.0000 

        """
        if os.path.exists(outpath) and not overwrite:
            raise IOError('file already exists; {0}'.format(outpath))
        
        assert atom_type in ['atomic', 'charge']
        
        xlo, ylo, zlo = self._origin

        a, b, c = self._abc
        xhi = a[0] + xlo
        xy = b[0]
        yhi = b[1] + ylo
        xz = c[0]
        yz = c[1]
        zhi = c[2] + zlo        
        
        xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
        xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
        ylo_bound = ylo + min(0.0,yz)
        yhi_bound = yhi + max(0.0,yz)
        zlo_bound = zlo
        zhi_bound = zhi
        
        num_atoms = self._atom_df.shape[0] 
        types = self._atom_df['type'].unique()
        num_types = len(types)
        type_map = dict(zip(types, [i+1 for i in range(len(types))]))

        with open(outpath, 'w+') as f:
            
            # header comments 
            f.write('# This file was created by ipymd (v{0}) on {1} \n'.format(
                    __version__, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            f.write('# type map: {0} \n'.format(type_map))
            f.write('# {0} \n'.format(header))
            
            # header
            f.write('\n')
            f.write('{0} atoms \n'.format(num_atoms))
            f.write('{0} atom types \n'.format(num_types))
            f.write('\n')
            
            f.write('# simulation box boundaries\n')
            f.write('{0:.4f} {1:.4f} xlo xhi \n'.format(xlo_bound, xhi_bound))
            f.write('{0:.4f} {1:.4f} ylo yhi \n'.format(ylo_bound, yhi_bound))
            f.write('{0:.4f} {1:.4f} zlo zhi \n'.format(zlo_bound, zhi_bound))
            f.write('{0:.4f} {1:.4f} {1:.4f} xy xz yz \n'.format(xy, xz, yz))            
            f.write('\n')
            
            # body
            f.write('Atoms \n')
            f.write('\n')
            
            for i, (ix, s) in enumerate(self._atom_df.iterrows()):
                if atom_type == 'atomic':
                    f.write('{0} {1} {3:.4f} {4:.4f} {5:.4f} \n'.format(
                    i+1, type_map[s.type], *s[['xs','ys','zs']].values))
                elif atom_type == 'charge':
                    f.write('{0} {1} {2:.4f} {3:.4f} {4:.4f} {5:.4f} \n'.format(
                    i+1, type_map[s.type], *s[['q','xs','ys','zs']].values))
                        


            