# -*- coding: utf-8 -*-
"""
Created on Wed May 18 22:24:10 2016

@author: cjs14

adapted from from http://physics.bu.edu/~erikl/research/tools/crystals/read_cif.py
"""
import os
import math

import numpy as np
import pandas as pd
    
from .base import DataInput

class CIF(DataInput):
    """ Build a crystal from  a Crystallographic Information File (.cif)

    """ 
    def __init__(self, file_path, ignore_overlaps=False):
        """ Build a crystal from  a Crystallographic Information File (.cif)

        here is a typical example of a CIF file:
            
             _cell_length_a 4.916
             _cell_length_b 4.916
             _cell_length_c 5.4054
             _cell_angle_alpha 90
             _cell_angle_beta 90
             _cell_angle_gamma 120
             _cell_volume 113.131
             _exptl_crystal_density_diffrn      2.646
             _symmetry_space_group_name_H-M 'P 32 2 1'
             loop_
             _space_group_symop_operation_xyz
               'x,y,z'
               'y,x,2/3-z'
               '-y,x-y,2/3+z'
               '-x,-x+y,1/3-z'
               '-x+y,-x,1/3+z'
               'x-y,-y,-z'
             loop_
             _atom_site_label
             _atom_site_fract_x
             _atom_site_fract_y
             _atom_site_fract_z
             Si   0.46970   0.00000   0.00000
             O   0.41350   0.26690   0.11910
    
        Parameters
        -----------
        file_path : str
            path to file
            
        """  
        data = self._read_cif_file(file_path)
        atoms_df, vectors = self._convert_cif_data(data, ignore_overlaps)
        self._add_colors(atoms_df)
        self._add_radii(atoms_df)
    
        self._atoms_df = atoms_df
        self._vectors = vectors
    
    def get_atom_data(self):
        """ return atom data """
        return self._atoms_df.copy()

    def get_simulation_box(self):
        """ return list of coordinates origin & [a,b,c] """
        return self._vectors

    #TODO read and deal with occupancy values (and overlapping atoms)        
    def _read_cif_file(self, file_path):
        """
        
        Read CIF file, and extract the necessary info in the form of a dictionary.
        E.g., the value of "_cell_volume" can be found with data['_cell_volume'].
        
        """
        assert os.path.exists(file_path), '{0} does not exist'.format(file_path)
        
        data = {}
        
        # Open the CIF file.
        with open(file_path, 'r') as f:
    
            reading_sym_ops = False
            
            atom_headers = []
    
            # Read lines one by one.
            for line in f:
    
                # Split into columns.
                cols = line.split()
                if (len(cols) == 0): continue
    
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Identify the keyword.  Here are the simply "key: value" ones.
                if (cols[0] == '_cell_length_a'):
                    data['_cell_length_a'] = float(cols[1])
    
                elif (cols[0] == '_cell_length_b'):
                    data['_cell_length_b'] = float(cols[1])
    
                elif (cols[0] == '_cell_length_c'):
                    data['_cell_length_c'] = float(cols[1])
    
                elif (cols[0] == '_cell_angle_alpha'):
                    data['_cell_angle_alpha'] = float(cols[1])
    
                elif (cols[0] == '_cell_angle_beta'):
                    data['_cell_angle_beta'] = float(cols[1])
    
                elif (cols[0] == '_cell_angle_gamma'):
                    data['_cell_angle_gamma'] = float(cols[1])
    
                elif (cols[0] == '_cell_volume'):
                    data['_cell_volume'] = float(cols[1])
    
    
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Extract the symmetry operations.  This will be a list of
                # strings such as:
                #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
                elif (cols[0] == '_space_group_symop_operation_xyz'):
                    reading_sym_ops = True
                    data['_space_group_symop_operation_xyz'] = []
    
                elif (reading_sym_ops):
                    
                    # Add the operation if the string is between single quotes.
                    # Otherwise it's a sign we are done with the list.
                    if (cols[0][0] == '\''  and  cols[0][-1] == '\''):
                        data['_space_group_symop_operation_xyz'].append(cols[0][1:-1])
                    else:
                        reading_sym_ops = False
                        # Note: it's safe to ignore this line completely.
    
    
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                # Search for the keyword "_atom_site_label" which indicates the
                # start of the atom_site data.
                elif (cols[0] == '_atom_site_label') and not atom_headers:
                    
                    while len(cols) < 4:
                        atom_headers.append(cols[0])
                        data[cols[0]] = []
                        line = self._skiplines(f)
                        cols = line.split()
                    
                    # Stop reading atom sites if we found a line with fewer
                    # columns, and which does not start with '_atom_site_'.
                    while len(cols) == len(atom_headers):
                        for i, name in enumerate(atom_headers):
                            data[name].append(cols[i])
                        line = self._skiplines(f)
                        cols = line.split()
        
        # Return the extracted data.
        return data
        
    def _extract_element(self, label):
        """Converts an "_atom_type_label" into an element name. """
    
        elem2 = ['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti',
                 'Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
                 'Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
                 'Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
                 'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt',
                 'Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa',
                 'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']
    
        if (label[0:2] in elem2):
            return label[0:2]
    
        elem1 = ['H','B','C','N','O','F','P','S','K','V','Y','I','W','U']
    
        if (label[0] in elem1):
            return label[0]
    
        raise Exception('WARNING: could not convert "%s" into element name!' % label)
        return label
    
    def _convert_cif_data(self, data, ignore_overlaps=False):
        
        radians = math.radians
        cos, sin = math.cos, math.sin
        sqrt = math.sqrt
                
        # Extract lengths and angles from the CIF data.
        La = float(data['_cell_length_a'])
        Lb = float(data['_cell_length_b'])
        Lc = float(data['_cell_length_c'])
        alpha = radians(float(data['_cell_angle_alpha']))
        beta = radians(float(data['_cell_angle_beta']))
        gamma = radians(float(data['_cell_angle_gamma']))
        volume = float(data['_cell_volume'])
        
        # Extract the symmetry operations.  This will be a list of strings such as:
        #    ['x,y,z', 'y,x,2/3-z', '-y,x-y,2/3+z', '-x,-x+y,1/3-z', ... ]
        ops = data['_space_group_symop_operation_xyz']
        
        # For proper evaluation, we need to convert "2/3" into "2./3", etc. to prevent
        # integer division which would turn e.g. 2/3 into 0.
        for i in range(len(ops)):
            ops[i] = ops[i].replace("/", "./")
        
        
        # Get the atom labels and coordinates.
        labels = data['_atom_site_label']
        fX = [ float(s) for s in data['_atom_site_fract_x'] ]
        fY = [ float(s) for s in data['_atom_site_fract_y'] ]
        fZ = [ float(s) for s in data['_atom_site_fract_z'] ]
        if data.has_key('_atom_site_occupancy'):
            occ = [ float(s) for s in data['_atom_site_occupancy'] ]
        else:
            occ = [ 1.0 for _ in data['_atom_site_label'] ]
        
        # Create a list of 4-tuples, where each tuple is an atom:
        #   [ ('Si', 0.4697, 0.0, 0.0),  ('O', 0.4135, 0.2669, 0.1191),  ... ]
        atoms = [ (labels[i], fX[i], fY[i], fZ[i], occ[i]) for i in range(len(labels)) ]
        
        # Make sure that all atoms lie within the unit cell.  Also convert names such
        # as 'Oa1' into 'O'.
        for i in range(len(atoms)):
            (name,xn,yn,zn,oc) = atoms[i]
            xn = (xn + 10.0) % 1.0
            yn = (yn + 10.0) % 1.0
            zn = (zn + 10.0) % 1.0
            name = self._extract_element(name)
            atoms[i] = (name,xn,yn,zn,oc)
            
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Use symmetry operations to create the unit cell.
        
        
        # The CIF file consists of a few atom positions plus several "symmetry
        # operations" that indicate the other atom positions within the unit cell.  So
        # using these operations, create copies of the atoms until no new copies can be
        # made.
        
        
        # Two atoms are on top of each other if they are less than "eps" away.
        eps = 0.01  # in Angstrom
        
        
        # For each atom, apply each symmetry operation to create a new atom.
        overlap_atoms = []
        imax = len(atoms)
        i=0
        while (i < imax):
        
            label,x,y,z,oc = atoms[i]
        
            for op in ops:
        
                # Python is awesome: calling e.g. eval('x,y,1./2+z') will convert the
                # string into a 3-tuple using the current values for x,y,z!
                xn,yn,zn = eval(op)
        
                # Make sure that the new atom lies within the unit cell.
                xn = (xn + 10.0) % 1.0
                yn = (yn + 10.0) % 1.0
                zn = (zn + 10.0) % 1.0
        
                # Check if the new position is actually new, or the same as a previous
                # atom.
                new_atom = True
                for at in atoms:
                    if (abs(at[1]-xn) < eps  and  abs(at[2]-yn) < eps  and  abs(at[3]-zn) < eps):
                        new_atom = False
        
                        # Check that this is the same atom type.
                        if (at[0] != label):
                            if at[4] + oc != 1:
                                raise Exception('overlapping atoms do not have occupancy summing to unity')
                            overlap_atoms.append((label,xn,yn,zn,oc))
                                
                # If the atom is new, add it to the list!
                if (new_atom):
                    atoms.append( (label,xn,yn,zn,oc) )  # add a 4-tuple
        
        
            # Update the loop iterator.
            i = i + 1
            imax = len(atoms)

        if overlap_atoms and not ignore_overlaps:
            raise Exception('invalid CIF file: atom of type %s overlaps with atom of type %s' % (at[0],label))
        
        atoms+=overlap_atoms

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Convert the fractional coordinates into real coordinates.
        
        
        # The primitive vectors a,b,c are such that 
        #
        #   cos(alpha) = b.c / |b||c|
        #   cos(beta)  = a.c / |a||c|
        #   cos(gamma) = a.b / |a||b|
        #
        # with the convention
        #
        #   a = La*xhat
        #   b = bx*xhat + by*yhat
        #   c = cx*xhat + cy*yhat + cz*zhat
        #
        cosa = cos(alpha)
        #sina = sin(alpha)
        cosb = cos(beta)
        #sinb = sin(beta)
        cosg = cos(gamma)
        sing = sin(gamma)
        
        cosa2 = cosa * cosa
        cosb2 = cosb * cosb
        sing2 = sing * sing
        
        ax = La
        
        bx = Lb * cosg
        by = Lb * sing
        
        cx = Lc * cosb
        cy = Lc * (cosa - cosg*cosb) / sing
        cz = Lc * sqrt( 1 - (cosa2 + cosb2 - 2*cosg*cosb*cosa) / sing2 )
        
        # Use the volume to check if we did the vectors right.
        V = ax*by*cz
        if ( abs(V - volume) > 0.1):
            raise Exception('volume does not match that calculated from primitive vectors')
        
        a = np.array([ax,0,0])
        b = np.array([bx,by,0])
        c = np.array([cx,cy,cz])
    
        for i in range(len(atoms)):
        
            # Get label and fractional coordinates.
            label,xf,yf,zf,oc = atoms[i]
        
            xa = xf * ax  # contribution of a-vector to the x-coordinate of this atom
            #ya = 0       # a-vector has no y-component, so does not affect y of atom
            #za = 0       # a-vector has no z-component, so does not affect z of atom
            
            xb = yf * bx  # contribution of b-vector to the x-coordinate of this atom
            yb = yf * by  # contribution of b-vector to the y-coordinate of this atom
            #zb = 0       # b-vector has no z-component, so does not affect z of atom
        
            xc = zf * cx  # contribution of c-vector to the x-coordinate of this atom
            yc = zf * cy  # contribution of c-vector to the y-coordinate of this atom
            zc = zf * cz  # contribution of c-vector to the z-coordinate of this atom
        
            # Add all contributions.
            xn = xa + xb + xc
            yn = yb + yc
            zn = zc
        
            atoms[i] = (label, xn, yn, zn,oc)
        
        df = pd.DataFrame(atoms, columns=['type', 'xs', 'ys', 'zs','occupancy'])
    
        return df, (np.array([a,b,c]),np.array([0.,0.,0.]))   
        
        
