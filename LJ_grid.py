# -*- coding: utf-8 -*-
"""
@author: david
"""

import numpy as np

from LJ_mol import LJForce
class grid():
    
    def __init__(self, xlow, xup, ylow, yup, mols, grid_x_steps, grid_y_steps, cutoffdist=None, printflag=False):
        self.xlow=xlow
        self.xup=xup
        self.ylow=ylow
        self.yup=yup
        self.mols=mols
        self.printflag=printflag
        self.cutoffdist = cutoffdist if cutoffdist is not None else max((xup-xlow)/grid_x_steps, (yup-ylow)/grid_y_steps)
        self.grid_x_steps=grid_x_steps
        self.grid_y_steps=grid_y_steps
        self.gridx = np.linspace(xlow,xup,grid_x_steps)
        self.gridy = np.linspace(ylow,yup,grid_y_steps)
        self.n_part = len(mols)
        self.gridPartNum = np.zeros((grid_x_steps, grid_y_steps))
        self.pairlist = []
        x_step_size = (xup-xlow)/grid_x_steps
        y_step_size = (yup-ylow)/grid_y_steps
        self.cellNo = grid_x_steps*grid_y_steps
        ##position of which the first particle of i-cell should be stored in GridList
        self.gridPointer= np.zeros(self.cellNo)

        #GridList
        self.gridList = [0 for i in range(self.n_part)]
        self.debugpair=[]
        
        ##Assign cell IDs to each molecule
        ##initialize the gridPointer array
        serialcounter=0
        for mol in mols:
            mol.serialNo = serialcounter
            serialcounter+=1
            x = mol.get_pos()[0]
            y = mol.get_pos()[1]
            x_cell=0
            y_cell =0
            x_cell_lb=xlow
            x_cell_up = xlow+x_step_size
            for i in range(0,grid_x_steps):
                if  x_cell_lb< x and x_cell_up > x:
                    x_cell = i
                    mol.grid_x=i
                    break
                else:
                    x_cell_lb+= x_step_size
                    x_cell_up += x_step_size
                    
                    
            y_cell_lb=ylow
            y_cell_up=ylow+y_step_size
            for i in range(0,grid_y_steps):
                if y_cell_lb<y and y_cell_up > y:
                    y_cell=i
                    mol.grid_y=i
                    break
                else:
                    y_cell_lb+= y_step_size
                    y_cell_up += y_step_size
                    
            mol.cell_id = y_cell*(grid_x_steps)+(x_cell)
            self.gridPointer[mol.cell_id]+=1
            self.gridPartNum[x_cell,y_cell]+=1
            
            
        ##modify gridPointer to correspond to correct indices
        gridPointerCumSum = np.cumsum(self.gridPointer)
        self.gridPointer2 = [0]
        for i in range(0,len(self.gridPointer)-1):
            self.gridPointer2.append(int(gridPointerCumSum[i]))
        
        ##we'll need gridpointer2 again, so use a temp to initialize gridList
        tempGridPtr=self.gridPointer2[:]
        
        
        ##construct gridlist    
        for mol in mols:
            self.gridList[int(tempGridPtr[mol.cell_id])]=mol
            tempGridPtr[mol.cell_id]+=1
            
            #self.gridList[int(self.gridPointer2[mol.cell_id])]=mol
            #self.gridPointer2[mol.cell_id]+=1
            
    def calc_pairlist(self):
        ##Algorithm 1 of Watanabe-Suzuki-Ito
        for i in range(self.cellNo):
            if i%100==0:
                if self.printflag: print('cell i= '+str(i))
           
            si = int(self.gridPointer2[i])
            ei = int(si+self.gridPointer[i])
            Lc = self.gridList[si:ei]
           
            neighbors = self.findNeighbors(i)
            neighbors = [n for n in neighbors if n>i]
            
            for j in neighbors:
                sj=self.gridPointer2[j]
                ej=int(sj+self.gridPointer[j])
                Lc.extend(self.gridList[sj:ej])
               
            for k in range(0,int(self.gridPointer[i])):
                for ell in range(k+1, len(Lc)):
                    m=Lc[k]
                    n=Lc[ell]
                    dist = np.sqrt(np.sum([(m.pos[i]-n.pos[i])**2 for i in range(len(m.pos))]))
                    #print("COMPARING :"+str(m.get_serialNo())+" AND "+str(n.get_serialNo() ))
                    ##update pairlist IF THESE PARTICLES INTERACT
                    ##CHECK DIST VS RADIUS OF INFLUENCE HERE THEN APPEND TO PAIRLIST
                    if dist < self.cutoffdist: 
                        #print("ADDING PAIR: "+str(m.get_serialNo())+" AND "+str(n.get_serialNo()))
                        self.debugpair.append((m.get_serialNo(),n.get_serialNo()))
                        self.pairlist.append((m,n))
        if self.printflag: print("pairlist calculated")
        return 
    
    def findNeighbors(self, cellNo):
        cellCoord = (cellNo % self.grid_x_steps, int(np.floor(cellNo/self.grid_x_steps)))
        neighbors=[]
        for  i in [-1,0,1]:
            for j in [-1,0,1]:
                if i==0 and j==0:
                    continue
                if cellCoord[0]+i<0 or cellCoord[0]+i>=self.grid_x_steps:
                    continue
                if cellCoord[1]+j<0 or cellCoord[1]+j >= self.grid_y_steps:
                    continue
                neighbors.append((cellCoord[1]+j)*(self.grid_x_steps)+(cellCoord[0]+i))
        return neighbors
    
    
    def updateForce(self):
        for pair in self.pairlist:
            current_force = LJForce(pair[0],pair[1])
            pair[0].updateForce(current_force)
            pair[1].updateForce([-1*a for a in current_force])
        return
    
    
    def calcLJ(self):
        for pair in self.pairlist:
            continue
        return
            
                        
                    
        