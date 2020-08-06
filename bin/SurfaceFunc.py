import os, sys
import numpy as np

try:
    import pandas as pd
    PD_LAB = 1
except:
    print("pandas not present")
    PD_LAB = 0
    #exit()

import matplotlib.pyplot as mpl


from scipy.spatial import distance_matrix


try:
    from mayavi import mlab
    MLAB_LAB = 1
except:
    print("mayavi not present")
    MLAB_LAB = 0
    #exit()

DEB_FIND_ORIENT = 0



    


def RotateMatrix(c, s, ax):

    A = np.zeros((3,3))
    if(ax == 0):
        A[0,0] = 1
        A[1,1] = c
        A[1,2] = -s
        A[2,1] = s
        A[2,2] = c
    if(ax == 1):
        A[1,1] = 1
        A[0,0] = c
        A[0,2] = -s
        A[2,0] = s
        A[2,2] = c
    if(ax == 2):
        A[2,2] = 1
        A[0,0] = c
        A[1,0] = s
        A[0,1] = -s
        A[1,1] = c
    return(A)
        

    

def RotatePatch(patch, mean_normal_v,axis_z_or,  pin):
    
    XY = True
    YZ = True
    XZ = True

    _ESP_ = 1e-10

    patch_trans = patch[:,:3] - pin
    
    Npoints  =np.shape(patch_trans)[0]
    
    ### defining rotating vectors...
    r_z = np.array([0,0, axis_z_or])
    r_vn = mean_normal_v.copy() #np.mean(normal_v, axis=0) #nter_atom_pos[1,:]

    
    P = np.abs(1 - r_z.dot(r_vn)/np.sqrt((r_z.dot(r_z))*(r_vn.dot(r_vn))))

    r_vn1 = r_vn.copy()
    
    while P > _ESP_:

        #r_vn1 = r_vn.copy()

        r1 = np.array([r_vn1[0],r_vn1[1]])
        r2 = np.array([r_z[0],r_z[1]])
        r1/=np.sqrt(r1.dot(r1))
        #r2/=np.sqrt(r2.dot(r2))
        
        cos_theta = r1[1] / np.sqrt(r1.dot(r1)) #(r1[0]*r2[0] + r1[1]*r2[1])/r1.dot(r1)
        sin_theta = r1[0] / np.sqrt(r1.dot(r1)) #(r1[0]*r2[1] - r1[1]*r2[0])/r1.dot(r1)
        
        
        R = RotateMatrix(cos_theta, sin_theta, 2)
        if(XY):
            for i in range(Npoints):
                patch_trans[i,:] = np.dot(R, patch_trans[i,:]) 
            r_vn1 = np.dot(R, r_vn1)
        
        #print("xy")
        #tmp = r_vn1
        #print("nter",tmp/np.sqrt(tmp.dot(tmp)))
        
       
        ## y-z plane
           
        r1 = np.array([r_vn1[1],r_vn1[2]])
        r2 = np.array([r_z[1],r_z[2]])
        r1/=np.sqrt(r1.dot(r1))
    
        #cos_theta = r1.dot(r2)/r1.dot(r1)
        #sin_theta = (-r1[0]*r2[1] + r1[1]*r2[0])/r1.dot(r1)
        if(axis_z_or > 0 ):
            cos_theta = r1[1] #(r1[0]*r2[0] + r1[1]*r2[1])/r1.dot(r1)
            sin_theta = r1[0] #(r1[0]*r2[1] - r1[1]*r2[0])/r1.dot(r1)
        else:
            cos_theta = -r1[1] #(r1[0]*r2[0] + r1[1]*r2[1])/r1.dot(r1)
            sin_theta = -r1[0] #(r1[0]*r2[1] - r1[1]*r2[0])/r1.dot(r1)
      
            
        R = RotateMatrix(cos_theta, sin_theta, 0)
        if(YZ):
            for i in range(Npoints):
                patch_trans[i,:] = np.dot(R, patch_trans[i,:]) 
            r_vn1 = np.dot(R, r_vn1)
      

        P = np.abs(1 - r_z.dot(r_vn1)/np.sqrt((r_z.dot(r_z))*(r_vn1.dot(r_vn1))))
    return(r_vn1, patch_trans)
        


def myflip(m, axis):
    if not hasattr(m, 'ndim'):
        m = asarray(m)
    indexer = [slice(None)] * m.ndim
    try:
        indexer[axis] = slice(None, None, -1)
    except IndexError:
        raise ValueError("axis=%i is invalid for the %i-dimensional input array"
                         % (axis, m.ndim))
    return m[tuple(indexer)]

def IsolateSurfaces(surface, minD = 1.):
    '''
    This function groups points nearer than minD..
    '''

    # squaring distance to avoid sqrt..
    minD2 = minD**2

    # computing number of surface points..
    l, tmp = np.shape(surface)

    # initialinzing label vectors..
    surf_label = np.ones(l)
    surf_tmp = np.ones(l)

    # starting from label = 2
    lab = 2

    #surf_label[0] = lab
    #surf_tmp[0] = lab


    #pos_s = np.where(surf_tmp == lab)

    # computing number of points without label..
    Nleft = np.sum(surf_tmp != 0)

    # starting iterating over different surfaces..
    while(Nleft > 0):

        count = 1

        pos__ = np.where(surf_tmp != 0)

        # seeding: first unlabeled point takes lab label..
        surf_label[pos__[0][0]] = lab
        surf_tmp[pos__[0][0]] = lab



        # iterating to find points belonging to the same surface...
        while(count > 0):

            #print(count)
            count = 0
            pos_s = np.where(surf_tmp == lab)
            # creating mask for points still to be processed..
            mask = np.logical_and(surf_tmp > 0, surf_tmp != lab)

            for i in pos_s[0]:
                #print("pos",i)

                # computin distances between points and the i-esime point...
                d = (surface[:,0] - surface[i,0])**2 + (surface[:,1] - surface[i,1])**2 + (surface[:,2] - surface[i,2])**2

                #print(np.sort(d)[:10])

                # finding near points...
                mm = d < minD2

                # removing processed point from system..
                surf_tmp[i] = 0

                #print("00", np.sum(surf_tmp == 0),"out of", len(surf_tmp), "lab", lab)

                # creating mask for  point still to be processed..
                mmm = np.logical_and(mm,mask)

                N = np.sum(mmm)
                #print("N",N)
                if(N >0):
                    surf_label[mmm] = lab
                    surf_tmp[mmm] = lab
                    count += 1
                mask = np.logical_and(surf_tmp > 0, surf_tmp != lab)

        # creating a new label..
        lab += 1
        # looking for how many points still to be processed...
        Nleft = np.sum(surf_tmp != 0)
        #print("left",Nleft)
    return(surf_label)


def FindBorder(new_plane_ab):
    '''
    This function finds the border of a figure in the plane...
    '''
    
    a,b = np.shape(new_plane_ab)
    p_h = np.ones((a,b))
    p_v = np.ones((a,b))

    index = np.arange(0,a)
    #print(index)

    for i in range(a):

        # horizontal
        tmp = new_plane_ab[i,:]
        lr = tmp != 0
        if(len(index[lr]) > 0 ):
            pos_1 = index[lr][0]
            pos_2 = index[lr][-1]

            p_h[i,:pos_1] = 0
            p_h[i,(pos_2+1):] = 0
        else:
            p_h[i,:] = 0


        # vertical
        tmp = new_plane_ab[:,i]
        lr = tmp != 0
        if(len(index[lr]) > 0 ):
            pos_1 = index[lr][0]
            pos_2 = index[lr][-1]

            p_v[:pos_1,i] = 0
            p_v[(pos_2+1):,i] = 0
        else:
            p_v[:,i] = 0


    return(p_v*p_h)



def ContactPoints(list_1, list_2, thresh):
    '''
    This function finds the groups of point of list1 and list2 that have a distance lesser that thresh 
    from at least one point of the other list.
    '''
    
    thresh2 = thresh**2
    contact_1 = [0,0,0,0]
    contact_2 = [0,0,0]
    L1 = np.shape(list_1)[0]
    L2 = np.shape(list_2)[0]

    mmm = np.zeros(4)

    for i in range(L1):
        if(i%1000 == 0):
            sys.stderr.write("\rProgress %d / %d"%(i, L1))
        d2 = (list_1[i,0] - list_2[:,0])**2 + (list_1[i,1] - list_2[:,1])**2 + (list_1[i,2] - list_2[:,2])**2
        mask = d2 < thresh2

        if(np.sum(mask) > 0):

            mmm[:3] = list_1[i,:3]
            mmm[3] =  np.min(d2[mask])

            contact_1 = np.row_stack([contact_1, mmm])
            contact_2 = np.row_stack([contact_2, list_2[mask,:3]])

    try:
        #print("ok")
        contact_2 = np.unique(contact_2,axis=0)

    except:
        aaa = contact_2.tolist()
        output = [0,0,0]
        for x in aaa:
            if x not in output:
                output = np.row_stack([output, x])
        contact_2 = output.copy() #np.matrix(output)

    L1 = np.shape(contact_2)[0]
    mmm = []
    for i in range(L1):
        if(i%1000 == 0):
            sys.stderr.write("\rProgress %d / %d"%(i, L1))
        d2 = (list_1[:,0] - contact_2[i,0])**2 + (list_1[:,1] - contact_2[i,1])**2 + (list_1[:,2] - contact_2[i,2])**2
        mmm.append(np.min(d2))
    contact_2 = np.column_stack([contact_2, np.array(mmm)])

    return(contact_1[1:,:], contact_2[1:,:])

def ContactPoints_NewVersion(list_1, list_2, thresh):
    '''
    This function finds the groups of point of list1 and list2 that have a distance lesser that thresh 
    from at least one point of the other list.
    '''
    
    thresh2 = thresh**2
    contact_1 = [0,0,0,0]
    contact_2 = [0,0,0]
    L1 = np.shape(list_1)[0]
    L2 = np.shape(list_2)[0]

    list_index_1 = []
    list_index_2 = []

    indexes_l1 = np.arange(L1)
    indexes_l2 = np.arange(L2)
    
    
    mmm = np.zeros(4)

    for i in range(L1):
        if(i%1000 == 0):
            sys.stderr.write("\rProgress %d / %d"%(i, L1))
        d2 = (list_1[i,0] - list_2[:,0])**2 + (list_1[i,1] - list_2[:,1])**2 + (list_1[i,2] - list_2[:,2])**2
        mask = d2 < thresh2

        if(np.sum(mask) > 0):

            mmm[:3] = list_1[i,:3]
            mmm[3] =  np.min(d2[mask])

            list_index_1 = np.concatenate([list_index_1, [i]])
            list_index_2 = np.concatenate([list_index_2, indexes_l2[mask]])
            
            contact_1 = np.row_stack([contact_1, mmm])
            contact_2 = np.row_stack([contact_2, list_2[mask,:3]])

    list_index_2 = np.unique(list_index_2)
    try:
        #print("ok")
        contact_2 = np.unique(contact_2,axis=0)
  
    except:
        aaa = contact_2.tolist()
        output = [0,0,0]
        for x in aaa:
            if x not in output:
                output = np.row_stack([output, x])
        contact_2 = output.copy() #np.matrix(output)

    L1 = np.shape(contact_2)[0]
    mmm = []
    for i in range(L1):
        if(i%1000 == 0):
            sys.stderr.write("\rProgress %d / %d"%(i, L1))
        d2 = (list_1[:,0] - contact_2[i,0])**2 + (list_1[:,1] - contact_2[i,1])**2 + (list_1[:,2] - contact_2[i,2])**2
        mmm.append(np.min(d2))
    contact_2 = np.column_stack([contact_2, np.array(mmm)])

    return(contact_1[1:,:], contact_2[1:,:], list_index_1, list_index_2)


def FindRealBR(surface1, surface2, _threshold_):
    '''
    Obsolete function...
    '''
    
    dist_mat_AbAg = distance_matrix(surface1, surface2, p=2, threshold=1000000)
    point_1, point_2 = np.where( dist_mat_AbAg < _threshold_)

    dist_mat_AbAg[:, :] = 0
    dist_mat_AbAg[point_1, point_2] = 1.

    AB = np.sum(dist_mat_AbAg, axis = 1)
    AG = np.sum(dist_mat_AbAg, axis = 0)
    mask_ag = AG > 0
    mask_ab = AB > 0
    return(mask_ab,mask_ag)

def EigRotation(points, eigvec):

    rot = np.transpose(np.dot(np.transpose(eigvec), np.transpose(points)))
    return(rot)


def Plot3DPointsAndVectors(x,y,z, u,v,w, color = []):

    ll = len(x)
    mask = np.random.choice(np.arange(ll),100, replace = False)
    if(MLAB_LAB):
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 600))
        mlab.clf()
        pp = mlab.points3d(x,y,z, scale_factor=0.2)
        pv = mlab.quiver3d(x[mask],y[mask],z[mask], u[mask],v[mask],w[mask])
        mm = mlab.points3d(0,0,0, scale_factor=1)

        pp.glyph.scale_mode = 'scale_by_vector'
        if(len(color) != 0):
            pp.mlab_source.dataset.point_data.scalars = color

        mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
        mlab.show()
    else:
        print("Library mahavi not present!")



def Plot3DPoints(x,y,z, color,size):

    if(MLAB_LAB):
        mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(800, 600))
        mlab.clf()
        pp = mlab.points3d(x,y,z, scale_factor=size)
        #mm = mlab.points3d(0,0,0, scale_factor=size*10)

        pp.glyph.scale_mode = 'scale_by_vector'
        pp.mlab_source.dataset.point_data.scalars = color

        mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))
        mlab.show()
    else:
        print("Library mahavi not present!")


def BuildCone(zmax, Ndisk):

    dz = zmax/float(Ndisk)
    z = 0

    n = 100
    res = [0,0,0]
    rad = np.linspace(0, 2*np.pi, n)
    for i in range(Ndisk):
        z += dz
        x = z*np.cos(rad)
        y = z*np.sin(rad)

        res = np.row_stack([res, np.column_stack([x,y,np.ones(n)*z])])
    return(res)

def ConcatenateFigPlots(list_):

    l = len(list_)
    res = list_[0]

    n, tmp = np.shape(list_[0])


    col_list = np.linspace(-100, 100, l)
    col = np.ones(n)*col_list[0]

    if(l>1):
        for i in range(1,l):
            res = np.row_stack([res, list_[i]])
            n, tmp = np.shape(list_[i])
            col = np.concatenate([col, np.ones(n)*col_list[i]])
    return(res, col)



def FixBridgeRealBS(patch_ab, patch_ag, Dpp):
    '''
    This function isolates the different groups of points in two given sets (patches) according to a cutoff distance Dpp.
    It associates each group to the closest group of the other set and returns a list of matched patches. 
    '''
    
    ## processing  patches to remove islands..
    #ab
    index_ab_bd_ = IsolateSurfaces(patch_ab, Dpp)

    val_ab, counts_ab = np.unique(index_ab_bd_,return_counts=True)
    print("val",val_ab,"counts", counts_ab)
    
    #ag
    index_ag_bd_ = IsolateSurfaces(patch_ag, Dpp)

    val_ag, counts_ag = np.unique(index_ag_bd_,return_counts=True)
    print("val", val_ag, "counts",counts_ag)
    
    
    # creating matrix of number of points in each group and label returned by IsolateSurface func.
    # Columns are ordered from the biggest to the smallest patch.
    tmp = np.row_stack([counts_ab, val_ab])
    
    #s_c_ab = np.flip(tmp[:, np.argsort(tmp[0, :])], axis=1)
    s_c_ab = myflip(tmp[:, np.argsort(tmp[0, :])], axis=1)
  
    tmp = np.row_stack([counts_ag, val_ag])
    
    s_c_ag = myflip(tmp[:, np.argsort(tmp[0, :])], axis=1)
    
    # finding center of mass of each group...
    cm_ab = []
    cm_ag = []
    for l in s_c_ab[1,:]:
        #print(l)
        tmp = patch_ab[index_ab_bd_ == l]
        cm_ab.append(np.mean(tmp[:,:3], axis = 0))
    for l in s_c_ag[1,:]:
        tmp = patch_ag[index_ag_bd_ == l]
        cm_ag.append(np.mean(tmp[:,:3], axis = 0))
        
    
    
    l_ag = np.shape(s_c_ag)[1]
    l_ab = np.shape(s_c_ab)[1]
    
    # computing distance matrix between the centers of the groups intra sets...
    D = np.zeros((l_ab, l_ag))
    for i in range(l_ab):
        for j in range(l_ag):
            D[i,j] = np.sum((np.array(cm_ab[i]) -  np.array(cm_ag[j]))**2)
        
    # associating groups according to the minimal distance...
    index_ab = []
    index_ag = []
    lmin = np.min([l_ag, l_ab])
    for i in range(lmin):
        if(lmin == l_ab):
            x = np.where(D[i,:] == np.min(D[i,:]))[0][0]
            index_ab.append(i)
            index_ag.append(x)
        else:
            x = np.where(D[:,i] == np.min(D[:,i]))[0][0]
            index_ab.append(x)
            index_ag.append(i)
        
        
    patch_ab_list = []
    patch_ag_list = []

    lab_ab = s_c_ab[1,index_ab]
    lab_ag = s_c_ag[1,index_ag]
    
    ## defining list of matched patches..
    for i in range(lmin):
        patch_ab_list.append(patch_ab[index_ab_bd_ == lab_ab[i]])
        patch_ag_list.append(patch_ag[index_ag_bd_ == lab_ag[i]])
        
    
    return(patch_ab_list, patch_ag_list)




class Surface:
    '''
    This class comprises a set of functions to analyse surface and patches.
    '''



    def __init__(self, surface, patch_num = 5, r0 = 11, theta_max = 45, real_br = []):

        if(type(surface) == str):

            self.surface = self.ReadSurface(surface)
        else:
            self.surface = surface

        #Nx, Ny, Nz  = np.shape(self.surface)

        #self.box_shape = np.array([Nx,Ny, Nz])



        self.patch_num = patch_num  ## number of patches to create
        self.r0 = r0                ## radius of the sphere to build the patch
        self.theta_max = theta_max  ## maximum degree of the cone



        self.radius_of_cilinder = 2
        self.threshold_on_layers = 5

        if(len(real_br) != 0):
            self.real_br = real_br # saving the mask     #self.surface[real_br,:]
        else:
            self.real_br = []

        #self.distance_matrix = distance_matrix(self.surface, self.surface, p=2, threshold=1000000)


    def EnlargePixels(self, plane, THRES = 300):


        nx, ny = np.shape(plane)
        p = plane.copy()

        while(nx < 400):
            tmp = np.zeros((nx*2, ny*2))
            tmp[::2,::2] = p
            tmp[1::2,::2] = p
            tmp[::2,1::2] = p
            tmp[1::2,1::2] = p

            nx, ny = np.shape(tmp)
            p = tmp.copy()
        return(p)

    def BuildPatch(self, point_pos, Dmin):

        d2 = (self.surface[:,0] - self.surface[point_pos,0])**2 + (self.surface[:,1] - self.surface[point_pos,1])**2 + (self.surface[:,2] - self.surface[point_pos,2])**2

        #mask = self.distance_matrix[point_pos, :] <= self.r0
        mask = d2 <= self.r0**2
        patch_points = self.surface[mask,:]



        ## processing patch to remove islands..
        index_ = IsolateSurfaces(patch_points, Dmin)

        val, counts = np.unique (index_,return_counts=True)
        #print(val, counts)
        pos_ = np.where(counts == np.max(counts))[0][0]
        lab_ = val[pos_]

        mmm = np.ones(len(mask))*-1000.
        mmm[mask] = index_
        new_mask = mmm == lab_

        patch_points__ = self.surface[new_mask,:]

        return(patch_points__, new_mask)

    def FindPatchOrientation(self,rot_protein, patch_mask):

        rot_a = rot_protein.copy()
        rot_p = rot_a[patch_mask,:]

        # finding center of rotated patch..
        cm = np.mean(rot_p[:,:3], axis=0)

        R =  self.radius_of_cilinder  ## AA   the cilinder radius
        THRES = self.threshold_on_layers  ## AA the threshold for different layers

        # translating surface to patch center..
        rot_a[:,:3] -= cm

        # finding points inside a cylinder of radius R having center in the origin and heigth along the z axis...
        D = np.sqrt(rot_a[:,0]**2 + rot_a[:,1]**2)
        # finding points of the surface inside the cylinder but not belonging to the patch...
        mask = np.logical_and(D <= R, np.logical_not(patch_mask))

        Z_ = rot_a[mask,2]


        count, zz = np.histogram(Z_, bins=20)
        zz = zz[:-1]

        if(DEB_FIND_ORIENT):
                # showing complex...
            mpl.figure()
            mpl.plot(zz, count)
            mpl.axvline(0)
            mpl.show()
            res, c = ConcatenateFigPlots(list_=[rot_a[:,:3], rot_a[patch_mask,:3]]) #, rot_a[mask,:]])
            Plot3DPoints(x=res[:,0], y=res[:,1], z=res[:,2], color=c, size=0.2)


        zz = zz[count > 0]


        check = zz > 0

        if(np.sum(check) == 0 ):
            orientation = 1.
        elif(np.sum(np.logical_not(check)) == 0 ):
            orientation = -1.
        else:

            positive = 0

            ## if points are a mix of positive and negative the first must be negative..
            negative = 1

            start = zz[0]
            for i in range(1,len(zz)):

                if(start < 0 and zz[i] >0 ):
                    positive += 1.
                    start = zz[i]

                if(start < 0):
                    if(zz[i] - start > THRES):
                        negative += 1
                    start = zz[i]
                elif(start > 0):
                    if(zz[i] - start > THRES):
                        positive += 1
                    start = zz[i]

            ## checking countings
            if(positive%2 == 0 and negative%2 != 0):
                orientation = 1
            elif(positive%2 != 0 and negative%2 == 0):
                orientation = -1.
            else:
                print("Attention! The code does not account for this case!\n Returened up orientation..")
                orientation = 2.

        return(orientation, zz)


    def CreatePlane(self, patch, z_c,Np = 20):

        _, lc = np.shape(patch)
        
        rot_p = patch.copy()

        # computing geometrical center..
        cm = np.mean(rot_p[:,:3], axis=0)

        # shifting patch to have the cone origin in [0,0,0]...
        rot_p[:,2] -= z_c

        #computing distances between points and the origin..
        weigths = np.sqrt(rot_p[:,0]**2 +  rot_p[:,1]**2 + rot_p[:,2]**2)

        # computing angles..
        thetas = np.arctan2(rot_p[:,1], rot_p[:,0])

        #computing distances in plane..
        dist_plane = np.sqrt(rot_p[:,0]**2 +  rot_p[:,1]**2)

        #computing the circle radius as the maximum distant point..
        R = np.max(dist_plane)*1.01

        # creating plane matrix..
        if(lc == 3):
            plane = np.zeros((Np,Np))
        else:
            plane = np.zeros((Np,Np), dtype = np.complex)
        #adapting points to pixels..
        rot_p[:,0] += R
        rot_p[:,1] -= R

        pos_plane = rot_p[:,:2] #np.abs(rot_p[:,:2])


        dR = 2.*R/Np
        rr_x = 0
        rr_y = 0
        for i in range(Np):
            rr_y = 0
            for j in range(Np):
                mask_x = np.logical_and(pos_plane[:,0]> rr_x, pos_plane[:,0]<= rr_x+dR)
                mask_y = np.logical_and(pos_plane[:,1]< -rr_y, pos_plane[:,1]>= -(rr_y+dR))
                mask = np.logical_and(mask_x, mask_y)
                if(len(weigths[mask]) >0):
                    w = np.mean(weigths[mask])
                    if(lc == 3):
                        plane[j,i] = w
                    else:
                        w_el = np.mean(patch[mask,3])
                        #print(w_el)
                        plane[j,i] = w + 1j*w_el 
                rr_y += dR
            rr_x += dR



        return(plane, weigths, dist_plane, thetas)

    def FindOrigin(self,rotated_pacth, CHECK = 0):
        '''
        This function finds the origin of a cone that inglobates the pacth with a maximum angle of 45 degrees.
        Input:
        - patch points (matrix)
        - CHECK, if 1 the cone is plotted.

        Output:
        - the origin (z-axis) of the cone

        TODO: generalize to a chosen degree...
        '''

        # copying patch matrix..ndO
        rot = rotated_pacth.copy()

        # computing distances of points from the origin (geometrical center) in the xy plane..
        dist_in_plane = np.sqrt(rot[:,0]**2 + rot[:,1]**2)

        # finding point with maximum distance..
        max_dist = np.max(dist_in_plane)

        pos_max = np.where(dist_in_plane == max_dist)[0]
        if(len(pos_max) >1):
            pos_max = pos_max[0]


        # translating patch to put the centre of the cone in the origin
        d = -max_dist + rot[pos_max, 2]
        #print("d",d)
        #print(np.shape(rot[:,2]))
        rot[:,2] -= d
        # looking for points outside the cone: their plane distance must be bigger than their z component..
        mask = dist_in_plane > np.abs(rot[:,2])


        # shifting the cone origin untill all points are inside the cone..
        while(np.sum(mask) != 0):

            ## finding the maximum distance only among points outside the cone
            new_d = np.max(dist_in_plane[mask])
            pos_max_new = np.where(dist_in_plane == new_d)[0]

            if(len(pos_max_new) >1):
                pos_max_new = pos_max_new[0]

            # shifting the patch..
            d = -new_d + rot[pos_max_new, 2]
            rot[:,2] -= d

            # looking for outer point outside..
            mask = dist_in_plane > np.abs(rot[:,2])

        # plotting cone + pacth..
        if(CHECK):
            cone = BuildCone(30, 50)
            all_ = np.row_stack([cone, rot])
            col = np.concatenate([np.ones(len(cone[:,0]))*-10,np.ones(len(rot[:,0]))*10])
            Plot3DPoints(all_[:,0],all_[:,1], all_[:,2], col, 0.4)

        # find new center of the pacth, corresponding to find the overal z shift
        new_cm = np.mean(rot[:,:3], axis=0)

        return(-new_cm[2])

    def FillTheInnerGap(self, plane_):
        '''
        This function fills the inner gaps (pixel with zero value) in the unit circle of a NxN plane.
        It replaces the zero pixel with the mean of the nearby pixels, only for those pixels that have non zero near pixels.
        Input:
        - plane (square matrix)
        
        Output:
        - Filled plane
        '''
      
        # copying plane..
        plane = plane_.copy()

        plane_bin = np.copy(plane)
        plane_bin[plane_bin != 0] = 1.

        # finding dimensions..
        xl, yl = np.shape(plane)

        # defining radius..
        r = int((xl-1)/2.)

        # defining radial..
        x, y = np.meshgrid(np.arange(-r,r+1), np.arange(-r,r+1))
        #Rmat = x**2 + np.flip(y, axis=0)**2
        Rmat = x**2 + myflip(y, axis=0)**2
        r2 = r**2

        # difining mask..
        tmp = np.zeros((3,3))
        tmp[1,1] = 1
        x_, y_ = np.where(tmp == 0)
        x_ -= 1
        y_ -= 1


        list_x, list_y = np.where(plane == 0)

        count = 0

        l_x = len(list_x)
        l_x_old = 2*l_x


        while(l_x < l_x_old):
            #print("old",l_x_old, "new", l_x)

            for i in range(l_x):
                x = list_x[i]
                y = list_y[i]

                if(Rmat[x,y] < r2):
                    XX = x_ + x
                    YY = y_ + y
                    mask_x = np.logical_and(XX>= 0, XX < xl)
                    mask_y = np.logical_and(YY>= 0, YY < yl)
                    mask = np.logical_and(mask_x, mask_y)
                    XX = XX[mask]
                    YY = YY[mask]
                    tmp = plane_bin[XX , YY]
                    count = np.sum(tmp)
                    if(count >= 6):
                        tmp = plane[XX , YY]
                        l__ = np.sum(tmp>0)
                        tmp = np.sum(tmp)/l__
                        plane[x,y] = tmp

            list_x, list_y = np.where(plane == 0)
            l_x_old = l_x
            l_x = len(list_x)
            plane_bin = np.copy(plane)
            plane_bin[plane_bin != 0] = 1.


        return(plane)

    def FillTheGap_everywhere(self, plane_):
        '''
        This function fills the gaps (pixel with zero value) in the unit circle of a NxN plane.
        It replaces the zero pixel with the mean of the nearby pixels.
        Input:
        - plane (square matrix)
        
        Output:
        - Filled plane
        '''
        plane = plane_.copy()

        xl, yl = np.shape(plane)

        r = int((xl-1)/2.)

        x, y = np.meshgrid(np.arange(-r,r+1), np.arange(-r,r+1))
        Rmat = x**2 + myflip(y, axis=0)**2
        r2 = r**2

        tmp = np.zeros((3,3))
        x_, y_ = np.where(tmp == 0)
        x_ -= 1
        y_ -= 1


        list_x, list_y = np.where(plane == 0)

        count = 0

        while(len(list_x) != 0 and count < 50):

            list_x, list_y = np.where(plane == 0)

            for i in range(len(list_x)):
                x = list_x[i]
                y = list_y[i]

                if(Rmat[x,y] < r2):

                    if((plane[x+1,y] !=0  and plane[x-1,y] != 0) or (plane[x,y+1] !=0  and plane[x,y-1] != 0) ):
                        tmp = plane[x_ + x , y_ + y]

                        l__ = np.sum(tmp!=0)
                        tmp = np.sum(tmp)/l__
                        plane[x,y] = tmp

            list_x, list_y = np.where(plane == 0)

            for i in range(len(list_x)):
                x = list_x[i]
                y = list_y[i]

                if(Rmat[x,y] < r2):

                    if((plane[x+1,y+1] !=0  and plane[x-1,y-1] != 0) or (plane[x-1,y+1] !=0  and plane[x+1,y-1] != 0) ):
                        tmp = plane[x_ + x , y_ + y]

                        l__ = np.sum(tmp!=0)
                        tmp = np.sum(tmp)/l__
                        plane[x,y] = tmp
            count += 1

        ### final cicle
        list_x, list_y = np.where(np.logical_and(plane == 0 , Rmat<r2))

        count = 0

        while(len(list_x) != 0 and count < 50):

            list_x, list_y = np.where(plane == 0)

            for i in range(len(list_x)):
                x = list_x[i]
                y = list_y[i]

                if(Rmat[x,y] < r2):

                    tmp = plane[x_ + x , y_ + y]

                    l__ = np.sum(tmp!=0)
                    if(l__ == 0):
                        l__ = 1
                    tmp = np.sum(tmp)/l__
                    plane[x,y] = tmp
            count += 1

        plane[Rmat>r2] = 0

        return(plane)



    def PatchReorient(self, patch_points):
        '''
        This function computes the eigenvalues of the covariance matrix of the given set of points and 
        rotates the point in order to have the xy plane aligned to the two maximum covariance axes.
        Input:
        - point matrix
        
        Output:
        - rotated point matrix
        - eigenvectors 3x3 matrix (can be used as rotation matrix)
        '''
      
        cm_patch =  np.mean(patch_points[:,:3], axis=0)
        patch_points[:,:3] -= cm_patch
      
        cov_mat = np.cov(np.transpose(patch_points[:,:3]))
        eig, eigvec = np.linalg.eig(cov_mat)
        idx = eig.argsort()[::-1]
        eig = eig[idx]
        eigvec = eigvec[:,idx]

        rot_p = EigRotation(patch_points[:,:3], eigvec)

        return(rot_p, eig, eigvec)


    def PatchReorientNew(self,patch_points, verso):
    
        ll = np.shape(patch_points)[0]
    
        mean_v = np.mean(patch_points[:,3:6], axis=0)
        pin = np.mean(patch_points[:,:3], axis=0)
    
        res, c11 = ConcatenateFigPlots(list_=[patch_points[:,:3], patch_points[:,:3] + patch_points[:,3:6]])

    
        phi, rot_pacth_all  = RotatePatch(res[:,:3], mean_v, verso,pin)
    
        rot_patch = rot_pacth_all[:ll,:3]
        rot_normal_vec = rot_pacth_all[ll:,:3] - rot_patch
    
        return(rot_patch, rot_normal_vec)
