# IPython Molecular Dynamics Package (ipymd)
Analysis of Molecular Dynamics output in the IPython Notebook

This package aims to provide a means of producing **reusable** analysis of Molecular Dynamics (MD) output in the IPython Notebook. 

There are many programs for 3D visualisation of MD output (my favourite being [Ovito](http://www.ovito.org/index.php)). However, there lacks a means to produce a more thorough, documented analysis of the data. IPython Notebooks are ideal for this type of analysis and so the objective of `ipymd` is to produce a Python package that can be used in conjuction with programmes like Ovito, to produce documented and reuseable analysis.  

The aim of `ipymd` is to produce IPython Notebooks that include:

- Static images of the simulations
- Plots of simulation data

It will build primarily on the chemlab package, that is an API layer on top of OpenGL.   

## Example Output

```python
%matplotlib inline
from ipymd.md_data import LAMMPS_Data
from ipymd.visualise_sim import Visualise_Sim
```


```python
data = LAMMPS_Data(
    sys_path='ipymd/lammps_test_data/system.dump',
    atom_path='ipymd/lammps_test_data/atom.dump')
```


```python
sys_data = data.get_system_data_all()
sys_data.head()
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
      <th>1</th>
      <td>200</td>
      <td>5880</td>
      <td>3.997601</td>
      <td>3.997601</td>
      <td>106784.302378</td>
      <td>2568.163297</td>
      <td>6.616167</td>
      <td>-576911.132565</td>
      <td>115.942920</td>
      <td>-576795.189644</td>
      <td>-572795.687453</td>
    </tr>
    <tr>
      <th>2</th>
      <td>400</td>
      <td>5880</td>
      <td>3.993996</td>
      <td>3.993997</td>
      <td>106591.809864</td>
      <td>1560.503603</td>
      <td>8.739034</td>
      <td>-576962.187377</td>
      <td>153.144425</td>
      <td>-576809.042952</td>
      <td>-574383.189834</td>
    </tr>
    <tr>
      <th>3</th>
      <td>600</td>
      <td>5880</td>
      <td>3.989881</td>
      <td>3.989882</td>
      <td>106372.285789</td>
      <td>364.540620</td>
      <td>8.262727</td>
      <td>-576965.242403</td>
      <td>144.797535</td>
      <td>-576820.444868</td>
      <td>-576254.921821</td>
    </tr>
    <tr>
      <th>4</th>
      <td>800</td>
      <td>5880</td>
      <td>3.986465</td>
      <td>3.986465</td>
      <td>106190.203677</td>
      <td>-586.959616</td>
      <td>7.597382</td>
      <td>-576960.911674</td>
      <td>133.137903</td>
      <td>-576827.773772</td>
      <td>-577736.783571</td>
    </tr>
    <tr>
      <th>5</th>
      <td>1000</td>
      <td>5880</td>
      <td>3.988207</td>
      <td>3.988208</td>
      <td>106283.055838</td>
      <td>128.391396</td>
      <td>4.990469</td>
      <td>-576921.379605</td>
      <td>87.453896</td>
      <td>-576833.925708</td>
      <td>-576634.915276</td>
    </tr>
  </tbody>
</table>
</div>




```python
ax = sys_data.plot('time','temp')
ax.set_xlabel('Time (fs)')
ax.set_ylabel('Temperature (K)');
```


![png](output_3_0.png)



```python
atom_data = data.get_atom_data(0)
atom_data.head()
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>type</th>
      <th>xs</th>
      <th>ys</th>
      <th>zs</th>
      <th>mass</th>
      <th>q</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>259</td>
      <td>1</td>
      <td>0.000000</td>
      <td>1.071430</td>
      <td>2.08333</td>
      <td>55.845</td>
      <td>2.231530e-12</td>
    </tr>
    <tr>
      <th>1</th>
      <td>267</td>
      <td>1</td>
      <td>0.000000</td>
      <td>0.357143</td>
      <td>2.08333</td>
      <td>55.845</td>
      <td>2.295430e-12</td>
    </tr>
    <tr>
      <th>2</th>
      <td>269</td>
      <td>1</td>
      <td>0.357143</td>
      <td>0.714286</td>
      <td>2.08333</td>
      <td>55.845</td>
      <td>2.219790e-12</td>
    </tr>
    <tr>
      <th>3</th>
      <td>271</td>
      <td>1</td>
      <td>0.714286</td>
      <td>1.071430</td>
      <td>2.08333</td>
      <td>55.845</td>
      <td>2.151480e-12</td>
    </tr>
    <tr>
      <th>4</th>
      <td>279</td>
      <td>1</td>
      <td>0.357143</td>
      <td>0.000000</td>
      <td>2.08333</td>
      <td>55.845</td>
      <td>2.308750e-12</td>
    </tr>
  </tbody>
</table>
</div>




```python
vis = Visualise_Sim()
vis.visualise(atom_data,xrot=45,yrot=45)
```




![png](output_5_0.png)



