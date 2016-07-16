# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 20:18:04 2016

Derived from:  
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
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


@author: chris sewell
"""
#TODO Averaging over (thermal) phase space
#TODO Analysis of triclinic cells
#TODO Calculation of structure factor coefficients from atom charge & type (rather than pre-defining ionic state)
#TODO fast fourier transform?
#TODO parallelization of fourier summation (http://ipyparallel.readthedocs.io/en/latest/intro.html#getting-started)

import math
import numpy as np
import pandas as pd

from ..shared import get_data_path
from . import data
from . import basic
from .. import plotting

def _set_thetas(min2theta=1.,max2theta=179.):
    """ set min and max angles to assess

    Properties
    ----------
    min2theta : float
        minimum 2 theta range to explore (degrees)
    max2theta : float
        maximum 2 theta range to explore (degrees)
    
    Returns
    -------
    min_theta : float
        minimum theta range to explore (radians)
    max_theta : float
        maximum theta range to explore (radians)

    """
    # Process angles
    min_theta = math.radians(min2theta) / 2.
    max_theta = math.radians(max2theta) / 2.
    return min_theta, max_theta

def _compute_rmesh(sim_abc, wlambda, min_theta, max_theta,
                  rspace=[1,1,1], periodic=[True,True,True], manual=False):
    """Compute full reciprocal lattice mesh
    
    Properties
    ----------
    sim_abc : numpy.array((3,3))
        a,b,c cell vectors (length units)
    wlambda : float
        radiation wavelength (length units)
        x-rays usually in the range 0.1 to 100 Angstroms
    min_theta : float
        minimum theta range to explore (radians)
    max_theta : float
        maximum theta range to explore (radians)
    rspace : list of floats
        parameters to multiply the spacing of the reciprocal lattice nodes 
        in the h, k, and l directions respectively
    periodic : list of bools
        whether periodic boundary in the h, k, and l directions respectively
    manual : bool
        use manual spacing of reciprocal lattice points based on the values of the c parameters 
        (good for comparing diffraction results from multiple simulations, but small c required).
        
    Returns
    -------
    rmesh : np.array((N,3))
        mesh of k points defining reciprocal lattice
    
    """
    # get cell parameters
    a,b,c,alpha,beta,gamma = basic.lattparams_bb(sim_abc)
    cell_lenghts = (a,b,c)
    if alpha!=90 or beta!=90 or gamma!=90:
        raise ValueError("Compute XRD does not work with triclinic structures")
    
    # maximum reciprocal lattice vector |K|, 
    # calculated from Bragg's law 
    Kmax = 2 * math.sin(max_theta) / wlambda 
        
    # Calculate spacing between reciprocal lattice points
    # Using distance based on periodic repeating distance
    if manual:
        inv_cell_lengths = [1.0,1.0,1.0]
    else:
        if not np.any(periodic):
            raise ValueError("Compute XRD must have at least one periodic boundary unless manual spacing specified")
        
        # calc inverse dimension for periodic directions
        inv_cell_lengths = [1./cell_lenghts[i] if periodic[i] else np.nan for i in range(3)]
        ave_inv = np.nanmean(inv_cell_lengths)

        # Use the average inverse dimensions for non-periodic directions
        inv_cell_lengths = [icl if periodic[i] else ave_inv for i,icl in enumerate(inv_cell_lengths)]
        
    # resolution (i.e. spacing) of reciprocal space points
    dK = [inv_cell_lengths[i] * rspace[i] for i in range(3)]
    # maximum integer value for K points in each dimension
    Knmax = [math.ceil(Kmax / dK[i]) for i in range(3)]
    
    # create the full reprocal lattice indices grid
    rmesh = np.mgrid[-Knmax[0]:Knmax[0]+1:1, 
               -Knmax[1]:Knmax[1]+1:1, 
               -Knmax[1]:Knmax[1]+1:1].reshape(3,-1).T
    
    # resize reciprocal mesh to correct spacing
    for i in range(3):
        rmesh[:,i] *= dK[i]
    
    return rmesh

def _restrict_rmesh(rmesh, wlambda, min_theta, max_theta):
    """filter mesh points to only those in Eswald's sphere 
    (and angular limits)
    
    Parameters
    ----------
    rmesh : np.array((N,3))
        mesh of k points defining reciprocal lattice
    wlambda : float
        radiation wavelength (length units)
        x-rays usually in the range 0.1 to 100 Angstroms
    min_theta : float
        minimum theta range to explore (radians)
    max_theta : float
        maximum theta range to explore (radians)

    Returns
    -------
    rmesh_sphere : np.array((N,3))
        mesh of k points defining reciprocal lattice, 
        retricted to Eswald's shere (and angular limits)
    k_mods : np.array((N,1))
         modulus for each k-point
    thetas : np.array((N,1))
        angles for each k-point (radians)

    """
    # calculate the length (squared) of each mesh vector
    K_sqr = np.sum(np.square(rmesh),axis=1)
    # select only mesh points within the Eswald sphere radius
    radius_mask = K_sqr * wlambda**2 <= 2**2 # i.e. (2sin(pi))^2
    # calculate the angle of each remaining mesh vector
    K = np.sqrt(K_sqr[radius_mask])
    theta = np.arcsin(wlambda * np.sqrt(K_sqr[radius_mask]) * 0.5)
    # select only mesh points within the angular limits
    angle_mask = np.logical_and(theta <= max_theta, theta >= min_theta)
    # return remaining mesh points
    return rmesh[radius_mask][angle_mask], K[angle_mask], theta[angle_mask]

def _get_sf_coeffs():
    datapath = get_data_path('xray_scattering_factors_coefficients.csv',
                             module=data)
    return pd.read_csv(datapath, index_col=0,comment='#')

def _calc_struct_factors(atoms_df,rmesh_sphere,wlambda,k_mods):
    """ calculate atomic scattering factors, fj, 
    for each atom at each reciprocal lattice point

    Parameters
    ----------
    atoms_df : pandas.DataFrame
        a dataframe of info for each atom, including columns; x,y,z,type
    rmesh_sphere : np.array((N,3))
        mesh of k points defining reciprocal lattice, 
        retricted to Eswald's shere (and angular limits)
    wlambda : float
        radiation wavelength (length units)
        x-rays usually in the range 0.1 to 100 Angstroms
    k_mods : np.array((N,1))
         modulus for each k-point, only required for calclating Lorentz-polarization factor

    Returns
    -------
    intensities : dict(np.array((N,1)))
         structure factor for each k-point (values), for each atom type (keys)

    """
    
    # calculate (|K|/2)^2
    K_2_sqr = (0.5*k_mods)**2

    # get the structure factor coefficients
    sf_coeffs_df = _get_sf_coeffs()
    
    struct_factors = {}
    for atype in atoms_df.type.unique():
        sfs = sf_coeffs_df.loc[atype]
        struct_factors[atype] = 0
        struct_factors[atype] += sfs.A1 * np.exp(-sfs.B1*K_2_sqr)
        struct_factors[atype] += sfs.A2 * np.exp(-sfs.B2*K_2_sqr)
        struct_factors[atype] += sfs.A3 * np.exp(-sfs.B3*K_2_sqr)
        struct_factors[atype] += sfs.A4 * np.exp(-sfs.B4*K_2_sqr)
        struct_factors[atype] += sfs.C
    
    return struct_factors

def _calc_intensities(atoms_df, rmesh_sphere, wlambda, struct_factors,
                     thetas=None,k_mods=None,use_Lp=True):
    """ calculate diffraction intensities for each atom at each reciprocal lattice point

    Parameters
    ----------
    atoms_df : pandas.DataFrame
        a dataframe of info for each atom, including columns; x,y,z,type
    rmesh_sphere : np.array((N,3))
        mesh of k points defining reciprocal lattice, 
        retricted to Eswald's shere (and angular limits)
    wlambda : float
        radiation wavelength (length units)
        x-rays usually in the range 0.1 to 100 Angstroms
    k_mods : np.array((N,1))
         modulus for each k-point, only required for calclating Lorentz-polarization factor
    thetas : np.array((N,1))
        angles for each k-point (radians), only required for calclating Lorentz-polarization factor
    use_Lp : bool
        switch to apply Lorentz-polarization factor

    Returns
    -------
    intensities : np.array((N,1))
         intensity for each k-point

    """
    # compute F(K)
    F = np.zeros(rmesh_sphere.shape[0]) + 0*1j
    for xyz,atype in zip(atoms_df[['x','y','z']].values, atoms_df.type):
        inner_dot = 2 * np.pi * np.dot(xyz,rmesh_sphere.T)
        F += struct_factors[atype] * (np.cos(inner_dot) + 1j*np.sin(inner_dot))
    # compute Lp(theta)
    if use_Lp:
        sin_thetas = 0.5*k_mods*wlambda
        Lp = (1+np.cos(2*thetas)**2)/(np.cos(thetas)*sin_thetas**2)
    else:
        Lp = 1.
    # calculate intensities
    return Lp*F*np.conjugate(F)/float(atoms_df.shape[0])

def compute_xrd(atoms_df, sim_abc,wlambda, min2theta=1.,max2theta=179.,
                rspace=[1,1,1], periodic=[True,True,True], manual=False, lp=True):
    r"""Compute predicted x-ray diffraction intensities for a given wavelength
    
    Properties
    ----------
    atoms_df : pandas.DataFrame
        a dataframe of info for each atom, including columns; x,y,z,type
    sim_abc : numpy.array((3,3))
        a,b,c cell vectors (length units)
    wlambda : float
        radiation wavelength (length units)
        x-rays usually in the range 0.1 to 100 Angstroms
    min2theta : float
        minimum 2 theta range to explore (degrees)
    max2theta : float
        maximum 2 theta range to explore (degrees)
    rspace : list of floats
        parameters to adjust the spacing of the reciprocal lattice nodes 
        in the h, k, and l directions respectively
    periodic : list of bools
        whether periodic boundary in the h, k, and l directions respectively
    manual : bool
        use manual spacing of reciprocal lattice points based on the values of the c parameters 
        (good for comparing diffraction results from multiple simulations, but small c required).
    lp : bool
        switch to apply Lorentz-polarization factor

    Returns
    -------
    2thetas : np.array((N,1))
        2theta angles for each k-point (degrees)
    intensities : np.array((N,1))
         intensity for each k-point
    
    Notes
    -----
    This is an implementation of the virtual x-ray diffraction pattern algorithm
    by Coleman *et al*. [ref1]_
    
    The algorithm proceeds in the following manner:

    1. Define a crystal structure by position (x,y,z) and atom/ion type.
    2. Define the x-ray wavelength to use
    3. Compute the full reciprocal lattice mesh
    4. Filter reciprocal lattice points by those in the Eswald's sphere
    5. Compute the structure factor at each reciprocal lattice point, for each atom type
    6. Compute the x-ray diffraction intensity at each reciprocal lattice point
    7. Group and sum intensities by angle

    The reciprocal k-point modulii are calculated from Bragg's law:    
    
    .. math::
    
        \left| {\mathbf{K}} \right| = \frac{1}{{d_{\text{hkl}} }} = \frac{2\sin \left( \theta \right)}{\lambda }
        
    and are restricted to within the Eswald's sphere, as illustrated:
        
    .. image:: ../images/xrd_mesh.jpg
    
    The atomic scattering factors, fj, accounts for the reduction in 
    diffraction intensity due to Compton scattering, with coefficients based on 
    the electron density around the atomic nuclei.
    
    .. math::
    
        f_j \left(\frac{\sin \theta}{\lambda}\right)
         = \left[ \sum\limits^4_i a_i \exp \left( -b_i \frac{\sin^2 \theta}{\lambda^2} \right)\right] + c 
         = \left[ \sum\limits^4_i a_i \exp \left( -b_i \left(\frac{\left| {\mathbf{K}} \right|}{2}\right)^2 \right)\right] + c
    
    The relative diffraction intensity from x-rays is computed at each 
    reciprocal lattice point through:
    
    .. math::
    
        I_x(\mathbf{K}) = Lp(\theta) \frac{F(\mathbf{K})F^*(\mathbf{K})}{N}
        
    such that:
    
    .. math::
    
        F ({\mathbf{K}} )= 
        \sum\limits_{j = 1}^{N} {f_{j}.e^{\left( {2\pi i \, {\mathbf{K}} \cdot {\mathbf{r}}_{j} } \right)}} 
        = \sum\limits_{j = 1}^{N} {f_j.\left[ \cos \left( 2\pi \mathbf{K} \cdot \mathbf{r}_j \right) + i \sin \left( 2\pi \mathbf{K} \cdot \mathbf{r}_j \right) \right]}

    and the Lorentz-polarization factor is:

    .. math::

        Lp(\theta) = \frac{1+\cos^2 (2\theta)}{\cos(\theta)\sin^2(\theta)}    
    
    References
    ----------
    
    .. [ref1] 1.Coleman, S. P., Sichani, M. M. & Spearot, D. E. 
        A Computational Algorithm to Produce Virtual X-ray and Electron Diffraction 
        Patterns from Atomistic Simulations. JOM 66, 408–416 (2014).


    """
    min_theta, max_theta = _set_thetas(min2theta,max2theta)
    rmesh = _compute_rmesh(sim_abc,wlambda,min_theta, max_theta,rspace, periodic, manual)
    rmesh_sphere, k_mods, thetas = _restrict_rmesh(rmesh,wlambda,min_theta, max_theta)
    struct_factors = _calc_struct_factors(atoms_df,rmesh_sphere,wlambda,k_mods)
    I = _calc_intensities(atoms_df,rmesh_sphere,wlambda,struct_factors,thetas,k_mods,use_Lp=lp)
    
    return np.degrees(2*thetas), I

#TODO new numpy (v1.11.0) has auto bin size selection
def plot_xrd_hist(ang2thetas, intensities, bins=180*100, wlambda=None,barwidth=None):
    """ create histogram plot of xrd spectrum
    
    Properties
    ----------
    barwidth : float or None
        if None then the barwidht will be the bin width
    
    Returns
    -------
    plot : ipymd.plotting.Plotting
        a plot object
    
    """
    I_hist, theta_edges = np.histogram(ang2thetas,bins=bins,
                                       weights=np.real(intensities),density=True)
    if barwidth is None:
        bin_width = (ang2thetas.max() - ang2thetas.min())/bins
    else:
        bin_width = barwidth
    theta_left = theta_edges[:-1]
    zero_mask = I_hist!=0
    
    plot = plotting.Plotter()
    plot.axes.bar(theta_left[zero_mask],I_hist[zero_mask],bin_width,
                     label=r'$\lambda = {0}$'.format(wlambda))
    plot.axes.set_xlabel(r'Scatteting Angle ($2 \theta$)')
    plot.axes.set_ylabel('Relative Intensity')
    plot.axes.set_yticklabels([])
    plot.axes.grid(True) 
    if wlambda is not None:
        plot.axes.legend(loc='upper right',framealpha=0.5)
    return plot
    