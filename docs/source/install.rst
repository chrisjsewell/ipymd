Quick Start
-----------------

.. image:: https://anaconda.org/cjs14/ipymd/badges/version.svg   
    :target: https://anaconda.org/cjs14/ipymd

`Anaconda <https://www.continuum.io/>`__ is recommended to create a
Python environment within which to use ipymd:

::

    conda create -n ipymd -c cjs14 ipymd
    source activate ipymd
    jupyter notebook

Currently the conda package is only available for OSX. For other operating systems, 
or to use the latest version from Github, the following environment should work:

::

    conda create -n ipymd python=2.7.11=0 numpy scipy matplotlib pandas ipython ipython-notebook pillow pyopengl pyqt six

If there are any issues, see the known working package dependancies list: 
https://github.com/chrisjsewell/ipymd/blob/master/working_dependencies_list_osx.txt 