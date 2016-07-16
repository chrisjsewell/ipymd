.. PyGauss documentation master file, created by
   sphinx-quickstart on Sun Jun 14 01:13:38 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Molecular Dynamics Analysis for IPython (ipyMD)!
===================================

This package aims to provide a means of producing **reusable** analysis of Molecular Dynamics (MD) output in the IPython Notebook. 

.. figure:: images/small_sim.png
    :width: 200px
    :align: center
    :height: 100px
    :alt: alternate text
    :figclass: align-center

    Analysis of the atomic coordination of Na, wrt Cl, for an NaCl nano-crystal.

There are many programs for 3D visualisation of MD output (my favourite being `Ovito <http://www.ovito.org>`__). However, there lacks a means to produce a more thorough, documented analysis of the data. IPython Notebooks are ideal for this type of analysis and so the objective of `ipymd` is to produce a Python package that can be used in conjuction with programmes like Ovito, to produce documented and reuseable analysis.  

The aim of `ipymd` is to produce IPython Notebooks that include:

- Static images of the simulations
- Analysis of simulation data, including graphical plots

It has been created with the goal to be:

- Easy to use
- Easy to extend
`chemlab <http://chemlab.readthedocs.io/>`__
It builds primarily on the `chemlab <http://chemlab.readthedocs.io/>`__ package, that is an API layer on top of OpenGL. Data is parsed in standard formats, such as [pandas](http://pandas.pydata.org/) dataframes, which are easy to create and use independantly from this package, in order to extend its functionality.  

+--------------------------+-----------------------------------------------+
|**Author**                | Chris Sewell                                  |
+--------------------------+-----------------------------------------------+
|**Project Page**          | https://github.com/chrisjsewell/ipymd         |
+--------------------------+-----------------------------------------------+

Contents
--------

.. toctree::
   :maxdepth: 3

   install
   tutorial
   package_api

License
-------

ipymd is released under the GNU GPL3 or GNU LGPL license, if the PyQt parts are omitted (in ipymd.visualise.opengl) and the ipymd.data_input.spacegroup package is omitted,
ipyMD is released under the `GNU GPLv3
<http://www.gnu.org/licenses/gpl.html>`_ and its main developer is
Chris Sewell.
