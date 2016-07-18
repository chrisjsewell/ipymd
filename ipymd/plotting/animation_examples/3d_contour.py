# -*- coding: utf-8 -*-
"""
Created on Sun Jul 17 21:56:26 2016

@author: cjs14

http://scicomp.stackexchange.com/questions/7030/plotting-a-2d-animated-data-surface-on-matplotlib
"""

from mpl_toolkits.mplot3d import axes3d
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
from ipymd.plotting.JSAnimation.IPython_display import display_animation

def generate(X, Y, phi):
    R = 1 - np.sqrt(X**2 + Y**2)
    return np.cos(2 * np.pi * X + phi) * R

fig = plt.figure()
ax = axes3d.Axes3D(fig)
#plt.close()

xs = np.linspace(-1, 1, 50)
ys = np.linspace(-1, 1, 50)
X, Y = np.meshgrid(xs, ys)
Z = generate(X, Y, 0.0)
wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
ax.set_zlim(-1,1)

def update(i, ax, fig):
    ax.cla()
    phi = i * 360 / 2 / np.pi / 100
    Z = generate(X, Y, phi)
    wframe = ax.plot_wireframe(X, Y, Z, rstride=2, cstride=2)
    ax.set_zlim(-1,1)
    return wframe,

ani = animation.FuncAnimation(fig, update, 
        frames=xrange(100), 
        fargs=(ax, fig), interval=100)
display_animation(ani)