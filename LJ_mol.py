# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 19:50:46 2023

@author: david
"""

import numpy as np
class molecule:
    ## molecule class
    def __init__(self, sigma, eps, pos):
        self.pos= pos
        self.sigma=sigma
        self.eps=eps
        self.grid_x=-1
        self.grid_y=-1
        self.force=np.zeros(len(pos))
    def get_sigma(self):
        return self.sigma
    def get_eps(self):
        return self.eps
    def get_pos(self):
        return self.pos
    def cell_id(self):
        return self.cell_id
    def get_serialNo(self):
        return self.serialNo
    def updateForce(self, newForce):
        for i in range(len(self.force)):
            self.force[i]+=newForce[i]
        return  


### Lorentz-Berthelot combination rule for finding parameters epsilon and sigma
### see https://en.wikipedia.org/wiki/Combining_rules
def LJ(mol1, mol2):
    eps = np.sqrt(mol1.get_eps()*mol2.get_eps())
    sig = (mol1.get_sigma()+mol2.get_sigma())/2
    dist = np.sqrt(np.sum([(mol1.get_pos()[i] - mol2.get_pos()[i])**2 for i in range(len(mol1.get_pos()))]))
    sixth = (sig/dist)**6
    return eps*(sixth**2 - 2*sixth)

def vLJ(mol1, vpos,vsig):
    ##calculates the LJ-potential of a molecule with a "virtual particle"
    ##could be used in calculating a "potential field"
    dist = np.sqrt(np.sum([(mol1.get_pos()[i] -vpos[i])**2 for i in range(len(mol1.get_pos()))]))
    sig = (vsig+mol1.get_sigma())/2
    return np.sqrt(mol1.get_eps())*((sig/dist)**12 - 2*(sig/dist)**6)

def LJForce(mol1, mol2):
    eps = np.sqrt(mol1.get_eps()*mol2.get_eps())
    sig = (mol1.get_sigma()+mol2.get_sigma())/2
    r = [mol1.get_pos()[i]-mol2.get_pos()[i] for i in range(len(mol1.get_pos()))]
    dist = np.sqrt(np.sum([(r[i])**2 for i in range(len(mol1.get_pos()))]))
    rhat = [x/dist for x in r]
    sixth = (sig/dist)**6
    LJForce = [4*eps*(12*sixth**2/dist -6*sixth/dist)*j for j in rhat]
    return LJForce
