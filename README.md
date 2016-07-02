# IPython Molecular Dynamics Analysis (ipyMD)

This package aims to provide a means of producing **reusable** analysis of Molecular Dynamics (MD) output in the IPython Notebook (http://ipymd.rtfd.io).


![png](docs/source/images/small_sim.png)


There are many programs for 3D visualisation of MD output (my favourite being [Ovito](http://www.ovito.org/index.php)). However, there lacks a means to produce a more thorough, documented analysis of the data. IPython Notebooks are ideal for this type of analysis and so the objective of `ipymd` is to produce a Python package that can be used in conjuction with programmes like Ovito, to produce documented and reuseable analysis.  

The aim of `ipymd` is to produce IPython Notebooks that include:

- Static images of the simulations
- Analysis of simulation data, including graphical plots

It has been created with the goal to be:

- Easy to use
- Easy to extend

It builds primarily on the [chemlab](http://chemlab.readthedocs.io/en/latest/) package, that is an API layer on top of OpenGL. Data is parsed in standard formats, such as [pandas](http://pandas.pydata.org/) dataframes, which are easy to create and use independantly from this package, in order to extend its functionality.  



## User Tutorial

### Instillation of Dependant Packages

[Anaconda](https://www.continuum.io/) is recommended to create a Python environment within which to use ipymd:

    conda create -n ipymd -c cjs14 ipymd
    source activate ipymd
    ipython notebook

Currently the conda package is only available for OSX. For other operating systems, see the conda build file for guidance on package dependancies: https://github.com/chrisjsewell/ipymd/blob/master/conda_recipe/meta.yaml

In the IPython Notebook, the ipymd package can then be imported.

```python
import ipymd
print ipymd.version()
```

    0.1.0


### Basic Atom Creation and Visualisation

The input for a basic atomic visualisation, is a [pandas](http://pandas.pydata.org/) Dataframe that specifies the coordinates, size and color of each atom in the following manner:


```python
import pandas as pd
df = pd.DataFrame(
        [[2,3,4,1,[0, 0, 255],1],
         [1,3,3,1,'orange',1],
         [4,3,1,1,'blue',1]],
        columns=['x','y','z','radius','color','transparency'])
```

Distances are measured in Angstroms, and colors can be defined in [r,g,b] format (0 to 255) or as a string defined in `available_colors`.


```python
print ipymd.available_colors()['reds']
```

    ['light_salmon', 'salmon', 'dark_salmon', 'light_coral', 'indian_red', 'crimson', 'fire_brick', 'dark_red', 'red']


The `Visualise_Sim` class can then be used to setup a visualisation, which is returned in the form of a `PIL` image.


```python
vis = ipymd.visualise_sim.Visualise_Sim()
vis.add_atoms(df)
img1 = vis.get_image(size=400,quality=5)
img1
```




    <PIL.Image._ImageCrop image mode=RGBA size=134x54 at 0x1581DBA70>



To convert this into an image viewable in IPython, simply parse it to the `visualise` function.


```python
vis.visualise(img1)
```




![png](docs/source/images/output_12_0.png)



Extending this basic procedure, additional objects can be added to the visualisation, the viewpoint can be rotated and multiple images can be output at once, as shown in the following example:


```python
vis.add_axes(length=0.2, offset=(-0.3,0))
vis.add_box([[5,0,0],[0,5,0],[0,0,5]])
vis.add_plane([[5,0,0],[0,5,2]],alpha=0.3)
vis.add_hexagon([[1,0,0],[0,0,.5]],[0,0,2],color='green')

img2 = vis.get_image(xrot=45, yrot=45)
#img2 = vis.draw_colorlegend(img2,1,2)
vis.visualise([img1,img2])
```




![png](docs/source/images/output_14_0.png)



Images can also be added to plots, with the `Plotting` class:


```python
plot = ipymd.plotting.Plotting()
plot.axes[0].scatter([0,0.5,1.2],[0,0.5,1])
plot.axes[0].grid(True)
plot.add_image_annotation(img2,(250,100),(0.5,0.5),zoom=0.5)
plot.resize_axes(width=0.5)
plot.display_plot(tight_layout=False)
```


![png](docs/source/images/output_16_0.png)



### Atom Creation From Other Sources

The `ipymd.data_input` module includes a number of classes to automate the intial creation of the atoms Dataframe, from various sources. Most classes will return a sub-class of `DataInput`, with a `get_atoms` method to return the atoms Dataframe and a `get_simulation_box` method to return the vertexes and origin of the simulation box.

#### Crystal Parameters

This class allows atoms to be created in ordered crystal, as defined by their space group and crystal parameters:


```python
data = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.4, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])

sim_abc, sim_origin = data.get_simulation_box()
print sim_origin
print sim_abc
atoms_df = data.get_atom_data()
atoms_df.head(2)
```

    [ 0.  0.  0.]
    [[  2.70000000e+01   0.00000000e+00   0.00000000e+00]
     [  1.65327318e-15   2.70000000e+01   0.00000000e+00]
     [  1.65327318e-15   1.65327318e-15   2.70000000e+01]]





<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>type</th>
      <th>x</th>
      <th>y</th>
      <th>z</th>
      <th>transparency</th>
      <th>color</th>
      <th>radius</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>Na</td>
      <td>0.000000e+00</td>
      <td>0.0</td>
      <td>0.0</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>Na</td>
      <td>3.306546e-16</td>
      <td>2.7</td>
      <td>2.7</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>




```python
vis2 = ipymd.visualise_sim.Visualise_Sim()
vis2.add_axes()
vis2.add_box(sim_abc, sim_origin)
vis2.add_atoms(atoms_df)
images = [vis2.get_image(xrot=xrot,yrot=45) for xrot in [0,45]]
vis2.visualise(images, columns=2)
```




![png](docs/source/images/output_22_0.png)



A dataframe is available which lists the alternative names for each space group:


```python
df = ipymd.data_input.crystal.get_spacegroup_df()
df.loc[[1,194,225]]
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>System_type</th>
      <th>Point group</th>
      <th>Short_name</th>
      <th>Full_name</th>
      <th>Schoenflies</th>
      <th>Fedorov</th>
      <th>Shubnikov</th>
      <th>Fibrifold</th>
    </tr>
    <tr>
      <th>Number</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>1</th>
      <td>triclinic</td>
      <td>1</td>
      <td>P1</td>
      <td>P 1</td>
      <td>$C_1^1$</td>
      <td>1s</td>
      <td>$(a/b/c)\cdot 1$</td>
      <td>-</td>
    </tr>
    <tr>
      <th>194</th>
      <td>hexagonal</td>
      <td>6/m 2/m 2/m</td>
      <td>P63/mmc</td>
      <td>P 63/m 2/m 2/c</td>
      <td>$D_{6h}^4$</td>
      <td>88a</td>
      <td>$(c:(a/a))\cdot m:6_3\cdot m$</td>
      <td>-</td>
    </tr>
    <tr>
      <th>225</th>
      <td>cubic</td>
      <td>4/m 3 2/m</td>
      <td>Fm3m</td>
      <td>F 4/m 3 2/m</td>
      <td>$O_h^5$</td>
      <td>73s</td>
      <td>$\left ( \frac{a+c}{2}/\frac{b+c}{2}/\frac{a+b...</td>
      <td>$2^{-}:2$</td>
    </tr>
  </tbody>
</table>
</div>



#### Crystallographic Information Files

.cif files are a common means to store crystallographic data and can be loaded as follows:


```python
cif_path = ipymd.get_test_path('example_crystal.cif')
data = ipymd.data_input.cif.CIF(cif_path)
vis = ipymd.visualise_sim.Visualise_Sim()
vis.basic_vis(data.get_atom_data(), data.get_simulation_box(),
              xrot=45,yrot=45)
```




![png](docs/source/images/output_27_0.png)



NB: at present, fractional occupancies of lattice sites are returned in the atom Dataframe, but cannot be visualised as such. It is intended that eventually occupancy will be visualised by partial spheres. 


```python
data.get_atom_data().head(1)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>type</th>
      <th>x</th>
      <th>y</th>
      <th>z</th>
      <th>occupancy</th>
      <th>transparency</th>
      <th>color</th>
      <th>radius</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>Fe</td>
      <td>4.363536</td>
      <td>2.40065</td>
      <td>22.642804</td>
      <td>1</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
    </tr>
  </tbody>
</table>
</div>



#### Lammps Input Data

The input data for LAMMPS simulations (supplied to `read_data`) can be input. Note that the `get_atom_data` method requires that the atom_style is defined, in order to define what each data column refers to.


```python
lammps_path = ipymd.get_test_path('lammps_input.data')
data = ipymd.data_input.lammps.LAMMPS_Input(lammps_path)

vis = ipymd.visualise_sim.Visualise_Sim()
vis.basic_vis(data.get_atom_data(atom_style='charge'), data.get_simulation_box(),xrot=45,yrot=45)
```




![png](docs/source/images/output_32_0.png)



#### Lammps Output Data

Output data can be read in the form of a single file or, it is advisable for efficiency, that a single file is output for each timestep, where `*` is used to define the variable section of the filename. The `get_atoms` and `get_simulation_box` methods not take a variable to define which timestep is returned.


```python
lammps_path = ipymd.get_test_path('atom_onefile.dump')
data = ipymd.data_input.lammps.LAMMPS_Output(lammps_path)

vis = ipymd.visualise_sim.Visualise_Sim()
vis.basic_vis(data.get_atom_data(98), data.get_simulation_box(98),
              spheres=True,xrot=45,yrot=45,quality=5)
```




![png](docs/source/images/output_35_0.png)




```python
lammps_path = ipymd.get_test_path(['atom_dump','atoms_*.dump'])
data = ipymd.data_input.lammps.LAMMPS_Output(lammps_path)
print data.count_timesteps()

vis = ipymd.visualise_sim.Visualise_Sim()
vis.basic_vis(data.get_atom_data(98), data.get_simulation_box(98),
              spheres=False,xrot=90,yrot=0)
```

    98





![png](docs/source/images/output_36_1.png)



### Atom Manipulation

The atoms Dataframe is already very easy to manipulate using the standard [pandas](http://pandas.pydata.org/) methods. But an `Atom_Manipulation` class has also been created to carry out standard atom manipulations, such as setting variables dependant on atom type or altering the geometry, as shown in this example:


```python
data = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.4, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])

manipulate_atoms = ipymd.atom_manipulation.Atom_Manipulation

new_df = manipulate_atoms(data.get_atom_data())

new_df.apply_radiimap({'Na':1.5, 'Cl':1})
new_df.apply_colormap({'Na':'blue','Cl':'green'})
new_df.change_type_variable('Na', 'transparency', 0.5)
new_df.slice_z(10,20)

vis2 = ipymd.visualise_sim.Visualise_Sim()
vis2.add_box(*data.get_simulation_box())
vis2.add_axes(offset=(-1.3,-0.7))
vis2.add_atoms(new_df.df, spheres=True)

img1 = vis2.get_image(xrot=45,yrot=45)

vis2.remove_atoms()
new_df.repeat_cell(data.get_simulation_box()[0],((-1,1),(-1,1),(-1,1)))
new_df.color_by_variable('z')
vis2.add_atoms(new_df.df, spheres=True)
img2 = vis2.get_image(xrot=90,yrot=0)

img3 = ipymd.plotting.create_colormap_image(new_df.df.z.min(), new_df.df.z.max(),
                                            horizontal=True,title='z position',size=150)

vis2.visualise([img1,img2, (280,1), img3], columns=2)
```




![png](docs/source/images/output_39_0.png)



NB: the default radii map is by atom  Van der Waals radii (`ipymd.atom_manipulation.vdw_dict`) and the default color map is by the same as in chemlab (`ipymd.atom_manipulation.default_atom_map`).

### Geometric Analysis

Given the simple and flexible form of the atomic data and visualisation, it is now easier to add more complex geometric analysis. These analyses are being contained in the `Atom_Analysis` class, and some initial examples are detailed below:

#### Atomic Coordination

The two examples below show computation of the coordination of Na, w.r.t Cl, in a simple NaCl crystal (which should be 6). The first does not include a consideration of the repeating boundary conditions, and so outer atoms have a lower coordination number. But the latter computation provides a method which takes this into consideration, by repeating the Cl lattice in each direction before computation.


```python
data = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.4, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])
df = data.get_atom_data()
df['coord'] = 0

analysis = ipymd.atom_analysis.Atom_Analysis()

coord_df = manipulate_atoms(df)
coord_df.filter_variables('Na')

lattice_df = manipulate_atoms(df)
lattice_df.filter_variables('Cl')

na_coord = analysis.calc_coordination(coord_df.df,lattice_df.df)

df.loc[df['type'] == 'Na','coord'] = na_coord

new_df = manipulate_atoms(df)
new_df.filter_variables('Na')
new_df.color_by_variable('coord',minv=3,maxv=7)

vis = ipymd.visualise_sim.Visualise_Sim()
vis.add_box(*data.get_simulation_box())
vis.add_atoms(new_df.df)

img = vis.get_image(xrot=45,yrot=45)

img2 = ipymd.plotting.create_legend_image(new_df.df.coord,new_df.df.color, title='Na Coordination',size=150,colbytes=True)

vis.visualise([img,img2],columns=2)
```

    //anaconda/envs/ipymd/lib/python2.7/site-packages/pandas/core/generic.py:2177: SettingWithCopyWarning: 
    A value is trying to be set on a copy of a slice from a DataFrame.
    Try using .loc[row_indexer,col_indexer] = value instead
    
    See the the caveats in the documentation: http://pandas.pydata.org/pandas-docs/stable/indexing.html#indexing-view-versus-copy
      self[name] = value





![png](docs/source/images/output_45_1.png)




```python
data = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.4, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])
df = data.get_atom_data()

df = analysis.calc_type_coordination(df, 'Na','Cl',repeat_vectors=data.get_simulation_box()[0])

new_df = manipulate_atoms(df)
new_df.filter_variables('Na')
new_df.color_by_variable('coord_Na_Cl',minv=3,maxv=7)

vis = ipymd.visualise_sim.Visualise_Sim()
vis.add_box(*data.get_simulation_box())
vis.add_atoms(new_df.df)

img = vis.get_image(xrot=45,yrot=45)

img2 = ipymd.plotting.create_legend_image(new_df.df.coord_Na_Cl,new_df.df.color, title='Na Coordination',size=150,colbytes=True)

vis.visualise([img,img2],columns=2)
```




![png](docs/source/images/output_46_0.png)



#### Atomic Structure Comparison

`compare_to_lattice` takes each atomic coordinate in df1 and computes the distance to the nearest atom (i.e. lattice site) in df2:


```python
import numpy as np
data1 = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.4, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])
df1 = data1.get_atom_data()

print ('Average distance to nearest atom (identical)', 
       np.mean(analysis.compare_to_lattice(df1,df1)))

data2 = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], ['Na', 'Cl'], 
    225, cellpar=[5.41, 5.4, 5.4, 90, 90, 90], 
    repetitions=[5, 5, 5])
df2 = data2.get_atom_data()

print ('Average distance to nearest atom (different)', 
       np.mean(analysis.compare_to_lattice(df1,df2)))

```

    ('Average distance to nearest atom (identical)', 0.0)
    ('Average distance to nearest atom (different)', 0.022499999999999343)


#### Common Neighbour Analysis (CNA)

CNA ([Honeycutt and Andersen, J. Phys. Chem. 91, 4950](http://dx.doi.org/10.1021/j100303a014)) is an algorithm to compute a signature for pairs of atoms, which is designed to characterize the local structural environment. Typically, CNA is used as an effective filtering method to classify atoms in crystalline systems ([Faken and Jonsson, Comput. Mater. Sci. 2, 279](http://dx.doi.org/10.1016/0927-0256(94%2990109-0)), with the goal to get a precise understanding of which atoms are associated with which phases, and which are associated with defects.

Common signatures for nearest neighbours are:

- FCC = 12 x 4,2,1
- HCP = 6 x 4,2,1 & 6 x 4,2,2
- BCC = 6 x 6,6,6 & 8 x 4,4,4
- Diamond = 12 x 5,4,3 & 4 x 6,6,3

which are tested below:


```python
data = ipymd.data_input.crystal.Crystal(
    [[0.0, 0.0, 0.0]], ['Al'], 
    225, cellpar=[4.05, 4.05, 4.05, 90, 90, 90], 
    repetitions=[5, 5, 5])
fcc_vector = data.get_simulation_box()[0]
fcc_df = data.get_atom_data()

data = ipymd.data_input.crystal.Crystal(
    [[0.33333,0.66667,0.25000]], ['Mg'], 
    194, cellpar=[3.21, 3.21, 5.21, 90, 90, 120], 
    repetitions=[5,5,5])
hcp_vector = data.get_simulation_box()[0]
hcp_df = data.get_atom_data()

data = ipymd.data_input.crystal.Crystal(
    [[0,0,0]], ['Fe'], 
    229, cellpar=[2.866, 2.866, 2.866, 90, 90, 90], 
    repetitions=[5,5,5])
bcc_vector = data.get_simulation_box()[0]
bcc_df = data.get_atom_data()

data = ipymd.data_input.crystal.Crystal(
    [[0,0,0]], ['C'], 
    227, cellpar=[3.57, 3.57, 3.57, 90, 90, 90], 
    repetitions=[2,2,2])
diamond_vector = data.get_simulation_box()[0]
diamond_df = data.get_atom_data()
```


```python
analysis= ipymd.atom_analysis.Atom_Analysis()
print analysis.cna_sum(fcc_df,repeat_vectors=fcc_vector)
print analysis.cna_sum(hcp_df,repeat_vectors=hcp_vector)
print analysis.cna_sum(bcc_df,repeat_vectors=bcc_vector)
print analysis.cna_sum(diamond_df,upper_bound=10,max_neighbours=16,repeat_vectors=diamond_vector)
```

    Counter({'4,2,1': 6000})
    Counter({'4,2,2': 1500, '4,2,1': 1500})
    Counter({'6,6,6': 2000, '4,4,4': 1500})
    Counter({'5,4,3': 768, '6,6,3': 256})


For each atom, the CNA for each nearest-neighbour can be output:


```python
analysis.common_neighbour_analysis(hcp_df,repeat_vectors=hcp_vector).head(5)
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>type</th>
      <th>x</th>
      <th>y</th>
      <th>z</th>
      <th>transparency</th>
      <th>color</th>
      <th>radius</th>
      <th>cna</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>Mg</td>
      <td>-0.000016</td>
      <td>1.853304</td>
      <td>1.3025</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
      <td>{u'4,2,2': 6, u'4,2,1': 6}</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>Mg</td>
      <td>1.605016</td>
      <td>0.926638</td>
      <td>3.9075</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
      <td>{u'4,2,2': 6, u'4,2,1': 6}</td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>Mg</td>
      <td>-0.000016</td>
      <td>1.853304</td>
      <td>6.5125</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
      <td>{u'4,2,2': 6, u'4,2,1': 6}</td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>Mg</td>
      <td>1.605016</td>
      <td>0.926638</td>
      <td>9.1175</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
      <td>{u'4,2,2': 6, u'4,2,1': 6}</td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>Mg</td>
      <td>-0.000016</td>
      <td>1.853304</td>
      <td>11.7225</td>
      <td>1</td>
      <td>light_salmon</td>
      <td>1</td>
      <td>{u'4,2,2': 6, u'4,2,1': 6}</td>
    </tr>
  </tbody>
</table>
</div>



This can be used to produce a plot identifying likely structure of an unknown structure:


```python
lammps_path = ipymd.get_test_path('thermalized_troilite.dump')
data = ipymd.data_input.lammps.LAMMPS_Output(lammps_path)
df = data.get_atom_data(0)
df = df[df.type==1]
plt = analysis.cna_plot(df,repeat_vectors=data.get_simulation_box(0)[0])
```


![png](docs/source/images/output_57_0.png)


A visualisation of the probable local character of eac atom can also be created. Note the *accuracy* parameter in the `cna_categories` method allows for more robust fitting to the ideal signatures:


```python
lammps_path = ipymd.get_test_path('thermalized_troilite.dump')
data = ipymd.data_input.lammps.LAMMPS_Output(lammps_path)

df = data.get_atom_data()
df = df[df.type==1]
df = analysis.cna_categories(df,repeat_vectors=data.get_simulation_box()[0],accuracy=0.7)
manip = ipymd.atom_manipulation.Atom_Manipulation(df)
manip.color_by_categories('cna')
#manip.apply_colormap({'Other':'blue','FCC':'green','HCP':'red'}, type_col='cna')
manip.change_type_variable('Other','transparency',0.5,type_col='cna')
atom_df = manip.df

vis = ipymd.visualise_sim.Visualise_Sim()
vis.add_box(*data.get_simulation_box())
vis.add_atoms(atom_df)

img = vis.get_image(xrot=45,yrot=45)

img2 = ipymd.plotting.create_legend_image(atom_df.cna,atom_df.color, 
                title='CNA Category\nof Fe Sublattice',size=150,colbytes=True)

vis.visualise([img,img2],columns=2)
```




![png](docs/source/images/output_59_0.png)



### System Analysis

Within the `LAMMPS_Output` class there is also the option to read in a systems data file, with a log of global variables for each simulation timestep.


```python
data = ipymd.data_input.lammps.LAMMPS_Output(
    sys_path=ipymd.get_test_path('system.dump'))
```


```python
sys_data = data.get_system_data()
sys_data.tail()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>time</th>
      <th>natoms</th>
      <th>a</th>
      <th>b</th>
      <th>vol</th>
      <th>press</th>
      <th>temp</th>
      <th>peng</th>
      <th>keng</th>
      <th>teng</th>
      <th>enth</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>94</th>
      <td>18800</td>
      <td>5880</td>
      <td>4.000919</td>
      <td>4.000920</td>
      <td>106961.684112</td>
      <td>3984.476256</td>
      <td>20.821004</td>
      <td>-576758.385663</td>
      <td>364.871078</td>
      <td>-576393.514586</td>
      <td>-570178.024689</td>
    </tr>
    <tr>
      <th>95</th>
      <td>19000</td>
      <td>5880</td>
      <td>4.002667</td>
      <td>4.002667</td>
      <td>107055.104770</td>
      <td>3457.822223</td>
      <td>20.543535</td>
      <td>-576743.788115</td>
      <td>360.008650</td>
      <td>-576383.779465</td>
      <td>-570985.120033</td>
    </tr>
    <tr>
      <th>96</th>
      <td>19200</td>
      <td>5880</td>
      <td>3.995696</td>
      <td>3.995696</td>
      <td>106682.545808</td>
      <td>1880.020740</td>
      <td>14.949236</td>
      <td>-576635.072775</td>
      <td>261.973143</td>
      <td>-576373.099632</td>
      <td>-573448.059018</td>
    </tr>
    <tr>
      <th>97</th>
      <td>19400</td>
      <td>5880</td>
      <td>3.990131</td>
      <td>3.990134</td>
      <td>106385.666301</td>
      <td>499.341106</td>
      <td>21.315529</td>
      <td>-576735.710976</td>
      <td>373.537223</td>
      <td>-576362.173753</td>
      <td>-575587.433134</td>
    </tr>
    <tr>
      <th>98</th>
      <td>19600</td>
      <td>5880</td>
      <td>3.979611</td>
      <td>3.979609</td>
      <td>105825.294988</td>
      <td>-2300.315883</td>
      <td>17.161320</td>
      <td>-576652.409368</td>
      <td>300.738098</td>
      <td>-576351.671270</td>
      <td>-579901.871560</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = sys_data.plot('time','temp')
ax.set_xlabel('Time (fs)')
ax.set_ylabel('Temperature (K)');
ax.grid()
```


![png](docs/source/images/output_64_0.png)

