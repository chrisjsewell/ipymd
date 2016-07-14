# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 20:18:04 2016

Derived from:  LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov
   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.
   https://github.com/lammps/lammps/tree/lammps-icms/src/USER-DIFFRACTION
   
This package contains the commands neeed to calculate x-ray and
electron diffraction intensities based on kinematic diffraction 
theory. Detailed discription of the computation can be found in the
following works:

Coleman, S.P., Spearot, D.E., Capolungo, L.  (2013) Virtual 
diffraction analysis of Ni [010] symmetric tilt grain boundaries, 
Modelling and Simulation in Materials Science and Engineering, 21 
055020. doi:10.1088/0965-0393/21/5/055020

Coleman, S.P., Sichani, M.M., Spearot, D.E.  (2014) A computational 
algorithm to produce virtual x-ray and electron diffraction patterns 
from atomistic simulations, JOM, 66 (3), 408-416. 
doi:10.1007/s11837-013-0829-3 

Coleman, S.P., Pamidighantam, S. Van Moer, M., Wang, Y., Koesterke, L. 
Spearot D.E (2014) Performance improvement and workflow development 
of virtual diffraction calculations, XSEDE14, 
doi:10.1145/2616498.2616552


@author: cjs14
"""



