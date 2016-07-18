# -*- coding: utf-8 -*-
"""
Created on Thu May 07 14:34:45 2015

@author: chris
"""
#!/usr/bin/env python

from distutils.core import setup

import os
def readme(file, git_path, img_folder):
    if os.path.exists(file):
        descript = open(file).read()
        descript = descript.replace('image:: ', 
                                    'image:: {0}/raw/master/{1}/'.format(git_path, img_folder))
        return descript
    return ''

import re
def version(path):
    verstrline = open(path, "rt").read()
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        raise RuntimeError("Unable to find version string in %s." % (path,))
    return verstr

setup(name='ipymd',
      version=version("ipymd/_version.py"),
      description='Analysis of Molecular Dynamics output in the IPython Notebooks',
      keywords = "chemistry physics molecular dynamics",
      long_description=readme('setup_README.rst',
                              'https://github.com/chrisjsewell/ipymd',
                              'docs/source/images'),
      author='Chris Sewell',
      author_email='chrisj_sewell@hotmail.com',
      url='https://http://ipymd.readthedocs.io',
      license = "GPL3",
      platforms = ["Any."],
      packages=['ipymd', 
                'ipymd.atom_analysis',
                'ipymd.atom_analysis.data',
                'ipymd.data_input',
                'ipymd.data_input.spacegroup',
                'ipymd.plotting',
                'ipymd.plotting.JSAnimation',
                'ipymd.plotting.JSAnimation.icons',
                'ipymd.shared',
                'ipymd.shared.atomdata',
                'ipymd.shared.fonts',
                'ipymd.test_data',
                'ipymd.test_data.atom_dump',
                'ipymd.visualise',
                'ipymd.visualise.opengl',
                'ipymd.visualise.opengl.postprocessing',
                'ipymd.visualise.opengl.postprocessing.shaders',
                'ipymd.visualise.opengl.renderers',
                'ipymd.visualise.opengl.renderers.opengl_shaders',
                ],
      package_data={'': ['*.rst', '*.txt'],
                    'ipymd.atom_analysis':['*.png','*.jpg'],
                    'ipymd.atom_analysis.data':['*.csv'],
                    'ipymd.data_input.spacegroup': ['*.csv','*.dat'],
                    'ipymd.plotting.JSAnimation.icons' : ['*.png'],
                    'ipymd.shared.atomdata':['*.txt'],
                    'ipymd.shared.fonts': ['*.ttf'],
                    'ipymd.test_data': ['*.dump', '*.cif','*.data'],
                    'ipymd.test_data.atom_dump': ['*.dump'],
                    'ipymd.visualise.opengl.postprocessing.shaders' : ['*.frag','*.vert'],
                    'ipymd.visualise.opengl.renderers.opengl_shaders' : ['*.frag','*.vert'],
                    },
      install_requires=[
                          "python=2.7.11=0",
                          "numpy",
                          "scipy",
                          "matplotlib",
                          "pandas",
                          "ipython",
                          "ipython-notebook",
                          "pillow",
                          "pyopengl",
                          "pyqt",
                          "six"
                       ],               
     )
