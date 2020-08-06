#!/usr/bin/env python
# coding: utf-8



import os, sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as mpl
import pandas as pd


sys.path.append("bin/")

import ZernikeFunc as ZF
import SurfaceFunc as SF


import time 


SAVE_FIG = 0
PLOT = 1

res_dir = "."


def EuclideanDistance(a,b):    
    d = np.sqrt(np.sum( (a-b)**2 ) )
    return(d)


def Zscore(x):
    z = (x-np.mean(x))/np.std(x)
    return(z)


if(len(sys.argv) != 4):
	print("Error!")
	print("Insert: pdb_code label_a label_b")
	exit()

pdb_name = sys.argv[1]
chain_a = sys.argv[2]
chain_b = sys.argv[3]



print("Loading complex %s"%(pdb_name))


try:
    surf_name1 =  "files/%s_%s_surf.csv"%(pdb_name, chain_a)
    surf_name2 =  "files/%s_%s_surf.csv"%(pdb_name, chain_b)
      
    
    surf_1_ = pd.read_csv(surf_name1)
    lag = len(surf_1_["x"])
    surf_1 = np.zeros((lag, 6))
    surf_1[:,:] = surf_1_[["x", "y", "z", "Nx", "Ny", "Nz"]]
    
    
    surf_2_ = pd.read_csv(surf_name2)
    lag = len(surf_2_["x"])
    surf_2 = np.zeros((lag, 6))
    surf_2[:,:] = surf_2_[["x", "y", "z", "Nx", "Ny", "Nz"]]
    
except:
    print("Cannot load %s not processed yet!"%(pdb_name))
    exit()
    


try:
    zern_1 = np.loadtxt("files/All_points_zernike_invariant_%s_%s_verso_up.dat"%(pdb_name, chain_a))
    zern_2 = np.loadtxt("files/All_points_zernike_invariant_%s_%s_verso_down.dat"%(pdb_name, chain_b))
except:
    print("Invariants fot complex %s are not presen!"%(pdb_name))
    exit()
    
l1 = np.shape(zern_1)[1]
l2 = np.shape(zern_2)[1]
    
    
d_all_theo = -2*np.dot(np.transpose(zern_1[1:,:]), zern_2[1:,:])
d1_square = np.sum( zern_1[1:,:]*zern_1[1:,:], axis=0)
d2_square = np.sum( zern_2[1:,:]*zern_2[1:,:], axis=0)

for i in range(l1):
    d_all_theo[i,:] += d2_square
for i in range(l2):
    d_all_theo[:,i] += d1_square
d_all_theo = np.sqrt(d_all_theo)    
   

color_1 = np.min(d_all_theo, axis=1)
color_2 = np.min(d_all_theo, axis=0)
    
    
len_1 = np.shape(surf_1)[0]
len_2 = np.shape(surf_2)[0]
    
    
# finding real patches
dd = 1. ## AA
patch_a, patch_b, mask_l1, mask_l2 = SF.ContactPoints_NewVersion(surf_1[:,:3], surf_2[:,:3], dd)
mask1 = np.zeros(len_1)
mask1[mask_l1.astype(int)] = 1

mask2 = np.zeros(np.shape(surf_2)[0])
mask2[mask_l2.astype(int)] = 1
    

np.savetxt("surf_newrot_%s_%s.txt"%(pdb_name, chain_a),np.column_stack([surf_1[:,:3], color_1, mask1]))
np.savetxt("surf_newrot_%s_%s.txt"%(pdb_name, chain_b),np.column_stack([surf_2[:,:3], color_2, mask2]))
    




