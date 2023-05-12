
from LJ_mol import molecule as molecule
import random
import numpy as np
import matplotlib.pyplot as plt
from LJ_grid import grid
import time

random.seed(10)
xlow = -1
xup = 1
ylow = -1 
yup = 1
n_part=2
sigma=0.5
eps=0.00000000000005

positions=[[0,0.05],[0,-0.05], [0.05,0.0]]

mols = [molecule(pos=positions[i], sigma=sigma, eps=eps) for i in range(len(positions))]
grid_x_steps = 2
grid_y_steps=2

grid = grid(xlow=xlow, xup=xup, ylow=ylow, yup=yup, mols=mols, grid_x_steps=grid_x_steps,
            grid_y_steps=grid_y_steps, cutoffdist=5)
grid.calc_pairlist()
#print("DEBUGLIST "+str(grid.debugpair))

for i in range(len(positions)):
   plt.scatter(positions[i][0],positions[i][1])
plt.xlim((-1,1))
plt.ylim((-1,1))


grid.updateForce()


V = [grid.mols[i].force for i in range(len(grid.mols))]
print(V)
origin = np.array([[grid.mols[i].pos[0] for i in range(len(grid.mols))],[grid.mols[i].pos[1] for i in range(len(grid.mols))]]) # origin point
xcoord = [V[i][0] for i in range(len(V))]
ycoord = [V[i][1] for i in range(len(V))]
plt.quiver(*origin,xcoord, ycoord, color=['r','b','g'], alpha=0.5, scale=1)


plt.gca().set_aspect('equal')
plt.show()


#print(grid.debugpair)
#print(grid.mols[0].force)





