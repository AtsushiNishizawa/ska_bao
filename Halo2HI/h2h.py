import numpy as np
import astropy.io.fits as pf
import h5py as h
import sys, os
import pylab as plt
from sklearn.neighbors import KDTree


def getcirc(x0,y0, r):
    t = np.linspace(0,np.pi*2.,100)
    x, y = x0+np.cos(t)*r, y0+np.sin(t)*r
    return x, y


#--parameters ---------------------------------------
Lbox      = 75000. # kpc/h
#fnum_list = [60, 68, 73, 85, 103, 111, 127, 45, 49, 54, 68, 85, 135]
fnum_list = [68]
nf        = 32
zlist     = np.loadtxt("/work/atsushi_data/illustris/illustris-3/zlist.txt")
zlist2    = {}
for i in range(zlist.shape[0]):
    key = "{0:03d}".format(int(zlist[i,0]))
    zlist2[key] = zlist[i,1]
#----------------------------------------------------

input_dir = "/work/atsushi_data/illustris/illustris-3/snapshots/"

#
# read x,y,z and vx,vy,vz from files for dark matter and gas
# 

for inum in fnum_list:

    fnum = "{0:03d}".format(inum)
    z   = zlist2[fnum]

    #"""
    for i in range(nf):
        fl = input_dir+"snap_"+fnum+".%i.hdf5"%(i)
        print ("    read ", fl)
        f = h.File(fl, "r")
        gas = f["PartType0"]
        if ( i == 0 ):
            x_gas = gas["Coordinates"             ].value
            nha   = gas["NeutralHydrogenAbundance"].value
            mass  = gas["Masses"                  ].value
        else:
            x_gas = np.append(x_gas, gas["Coordinates"             ].value, axis=0)
            nha   = np.append(nha  , gas["NeutralHydrogenAbundance"].value, axis=0)
            mass  = np.append(mass , gas["Masses"                  ].value, axis=0)
        f.close()
        del(gas)
    #"""


    for i in range(2):
        fgp = input_dir+"groups_%03i.%i.hdf5"%(inum, i)
        print ("    read ", fgp)
        f = h.File(fgp, "r")
        gp = f["Group"]
        if ( i == 0 ):
            xgp200  = gp["GroupPos"       ].value
            Mgp200c = gp["Group_M_Crit200"].value
            Rgp200c = gp["Group_R_Crit200"].value
        else:
            xgp200  = np.append(xgp200 , gp["GroupPos"       ].value, axis=0)
            Mgp200c = np.append(Mgp200c, gp["Group_M_Crit200"].value, axis=0)
            Rgp200c = np.append(Rgp200c, gp["Group_R_Crit200"].value, axis=0)
        f.close()
        del(gp)



    Ngp = Rgp200c.shape[0]
    print (Ngp)
    
    #objects near boundary
    #
    m = (Rgp200c<xgp200[:,0])*(xgp200[:,0]<Lbox-Rgp200c)*\
        (Rgp200c<xgp200[:,1])*(xgp200[:,1]<Lbox-Rgp200c)*\
        (Rgp200c<xgp200[:,2])*(xgp200[:,2]<Lbox-Rgp200c)
    print (m[~m].shape)
    
    #sanity check plots
    #
    if ( False ):
        plt.hist(Rgp200c, bins=100)
        plt.show()

        ic = 0
        for i in range(Ngp):
            if ( Rgp200c[i] > 10 ):
                plt.subplot(111,aspect=1.)
                x, y = getcirc(xgp200[i,0], xgp200[i,1], Rgp200c[i])
                plt.plot(x,y,ls="-", color="C%i"%(ic))
                plt.xlim(0,Lbox)
                plt.xlim(1,Lbox)
                ic += 1
                if ( ic == 10 ):
                    plt.show()
                    ic = 0
            

    #kd-tree finder
    #
    Ng = 30
    tree = KDTree(x_gas, leaf_size=300, metric="euclidean")
    ind = tree.query_radius(xgp200[:Ng,:], r=Rgp200c[:Ng])
    print (ind.shape)
    
    HI = mass*nha
    del(mass)
    del(nha)
    mHI = np.empty(Ng)
    for i in range(Ng):
        idx = ind[i]
        mHI[i] = np.sum(HI[idx])
    

    if ( True ):
        plt.plot(Mgp200c[:Ng], mHI, "C1o")
        plt.show()
