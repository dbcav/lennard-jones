# -*- coding: utf-8 -*-
"""
@author: david
"""

import LJ_mol as LJ
from LJ_mol import molecule as molecule
import random
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import LJ_grid 
from LJ_grid import grid


########Test pairlist
##initialize inputs
xlow = -1
xup = 1
ylow = -1 
yup = 1
n_part=9
sigma=1
eps=0.5
#positions = [[random.uniform(xlow,xup), random.uniform(ylow,yup)] for i in range(0,n_part)]
#positions = [[ -0.99,-0.99],[-0.94,-.82], [-0.55,-0.73], [0.5,-0.5],[0.82,-0.92], 
#             [0.75,0.75]],[-0.79,-0.92],[-0.90,-0.85],[0.75,0.33]]

positions = [[-0.75,-0.75], [-0.75,0],[-0.75,0.75],[0,-0.75],[0,0],[0,0.75],
             [0.75,-0.75],[0.75,0],[0.75,0.75]]
mols = [molecule(pos=positions[i], sigma=sigma, eps=eps) for i in range(len(positions))]
grid_x_steps = 3
grid_y_steps= 3

grid = grid(xlow=xlow, xup=xup, ylow=ylow, yup=yup, mols=mols, grid_x_steps=grid_x_steps,
            grid_y_steps=grid_y_steps)



grid.calc_pairlist()
print(grid.pairlist)
print(grid.debugpair)
print(len(grid.debugpair))
for i in range(len(positions)):
   plt.scatter(positions[i][0],positions[i][1])
plt.xlim((-1,1))
plt.ylim((-1,1))
