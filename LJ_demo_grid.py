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

##initialize inputs
xlow = -1
xup = 1
ylow = -1 
yup = 1
n_part=5
sigma=1
eps=0.5
positions = [[random.uniform(xlow,xup), random.uniform(ylow,yup)] for i in range(0,n_part)]
#positions = [[]]
mols = [molecule(pos=positions[i], sigma=sigma, eps=eps) for i in range(len(positions))]
grid_x_steps = 3
grid_y_steps= 3

##construct grid
grid = grid(xlow=xlow, xup=xup, ylow=ylow, yup=yup, mols=mols, grid_x_steps=grid_x_steps,
            grid_y_steps=grid_y_steps)

print(grid.gridPartNum)


print("find neighbors test: "+str(grid.findNeighbors(4)))
grid.calc_pairlist()
print(grid.pairlist)
for i in range(len(positions)):
   plt.scatter(positions[i][0],positions[i][1])
plt.xlim((-1,1))
plt.ylim((-1,1))