#!/usr/bin/env python
# coding: utf-8


import os, sys
import numpy as np
import scipy as sp

try:
    import matplotlib.pyplot as mpl
    PLOT_LAB = 1
except:
    print("Matplotlib is not present! Install it to see plots!")
    PLOT_LAB = 0

    
import pandas as pd
from scipy.spatial import distance_matrix

try:
    
    from mayavi import mlab
except:
    print("Mayavi library not present! Intall it to visulalize 3D plots!")
    exit()
 
# path to where libraries are..   
sys.path.append("bin/")


import ZernikeFunc as ZF
import SurfaceFunc as SF

import time 



Npixel = 25   # the plane edge in pixels...
Rs     = 6    # the radius of the sphere that includes the patch..
ZOrder = 20   # the Zernike expansion order..


if(len(sys.argv) != 4):
	print("Error! Invalid number of arguments!")
	print("Please, insert:")
	print("python ComputeZernike2DFromPatch.py   pdb_name  point_id  name_res")
	print("---- * ----")
	print("pdb_name, the file containing the surface to analyze;")
	print("point_id, the index of the point of the surface to be used as pin for creating the patch [0, #surface points -1];")
	print("name_res, the label to recognize the pacth (e.g. pdb code);")


pdb_name = sys.argv[1]
point_id = int(sys.argv[2])
name_res = sys.argv[3]

# loading patch files..
try:
        prot_ = pd.read_csv(pdb_name)
except:
        print("Protein %s not present!"%(name_res))
        exit()

    


lag = len(prot_["x"])
prot = np.zeros((lag, 6))
prot[:,:] = prot_[["x", "y", "z", "Nx", "Ny", "Nz"]]


#check!!
npoint = np.shape(prot)[0]
if( point_id < 0 or  point_id > npoint-1):
	print("Error! Invalid point index!")
	exit()
       

# creating surface objects
surf = SF.Surface(prot[:,:], patch_num = 0, r0 = Rs, theta_max = 45)


# extracting patch..
patch, mask = surf.BuildPatch(point_pos= point_id, Dmin=.5)
surf_obj.real_br = mask


# rotating patch in the xy plane..
rot_patch, rot_ag_patch_nv = surf.PatchReorientNew(patch, 1)
    

# finding cone origin...                
z = surf.FindOrigin(rot_patch)
            
            
           
# creating plane.. 
plane, weigths, dist_plane, thetas = surf.CreatePlane(patch=rot_patch, z_c=z , Np=Npixel)
new_plane = plane.copy()
new_plane = surf.FillTheGap_everywhere(plane_=plane)

## enlarging plane..
new_plane_ =  surf.EnlargePixels(new_plane)


try:
    zernike_env.img  = new_plane_
except:
    zernike_env = ZF.Zernike2d(new_plane_)

br_recon, br_coeff = zernike_env.ZernikeReconstruction(order=ZOrder, PLOT=0)
#br_coeff = zernike_env.ZernikeDecomposition(order=ZOrder)


coeff_inv = np.absolute(coeff)



if(PLOT_LAB):
    # plotting res..
    fig, ax = mpl.subplots(1,3, dpi = 150)
    ax[0].imshow(plane)
    ax[1].imshow(new_plane_)
    ax[2].imshow(recon.real)
    ax[0].axis('off')
    ax[1].axis('off')
    ax[2].axis('off')
    ax[0].set_title("original")
    ax[1].set_title("processed")
    ax[2].set_title("reconstructed")




    mpl.figure(dpi = 150)
    mpl.plot(coeff_inv)
    mpl.xlabel("index")
    mpl.ylabel("Invariant")
    mpl.show()



### saving invariants
np.savetxt("Zernike_inv_protein_%s_patch_%d.txt"%(name_res, point_id), coeff_inv)
