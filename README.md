# IPython Molecular Dynamics Package (ipymd)
Analysis of Molecular Dynamics output in the IPython Notebook

This package aims to provide a means of producing **reusable** analysis of Molecular Dynamics (MD) output in the IPython Notebook. 

There are many programs for 3D visualisation of MD output (my favourite being [Ovito](http://www.ovito.org/index.php)). However, there lacks a means to produce a more thorough, documented analysis of the data. IPython Notebooks are ideal for this type of analysis and so the objective of `ipymd` is to produce a Python package that can be used in conjuction with programmes like Ovito, to produce documented and reuseable analysis.  

The aim of `ipymd` is to produce IPython Notebooks that include:

- Static images of the simulations
- Plots of simulation data

It will build primarily on the chemlab package, that is an API layer on top of OpenGL.   
