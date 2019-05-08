import numpy as np
import pyfits as pf
import h5py as h
import sys, os
import matplotlib.pyplot as plt
import utils as ut
from scipy.fftpack import fftn, fftfreq

#--parameters ---------------------------------------
Lbox      = 75000. # kpc/h
fnum_list = [60, 68, 73, 85, 103, 111, 127]
nf        = 32
zlist     = np.loadtxt("zlist.txt")
zlist2    = {}
print zlist.shape
print zlist
for i in range(zlist.shape[0]):
    key = "{0:03d}".format(int(zlist[i,0]))
    zlist2[key] = zlist[i,1]
    print key, zlist[i,1]
#----------------------------------------------------

#
# read x,y,z and vx,vy,vz from files for dark matter and gas
# 

for inum in fnum_list:

    fnum = "{0:03d}".format(inum)
    fits = "illustris/illustris-3/cutouts/full_"+fnum+".fits"

    if ( os.path.exists(fits) ):
        print fits, "already exists"
    else:
        for i in range(nf):
            fl = "illustris/illustris-3/snapshots/snap_"+fnum+".%i.hdf5"%(i)
            print "read ", fl
            f = h.File(fl, "r")
            gas = f["PartType0"]
            dmp = f["PartType1"]
            z   = zlist2[fnum]
            if ( i == 0 ):
                x_gas = gas["Coordinates"].value
                v_gas = gas["Velocities" ].value
                nha   = gas["NeutralHydrogenAbundance"].value
                x_dmp = dmp["Coordinates"].value
                v_dmp = dmp["Velocities" ].value
            else:
                x_gas = np.append(x_gas, gas["Coordinates"].value, axis=0)
                v_gas = np.append(v_gas, gas["Velocities" ].value, axis=0)
                nha   = np.append(nha  , gas["NeutralHydrogenAbundance"].value, axis=0)
                x_dmp = np.append(x_dmp, dmp["Coordinates"].value, axis=0)
                v_dmp = np.append(v_dmp, dmp["Velocities" ].value, axis=0)
            f.close()
            del(gas)
            del(dmp)
        ut.writefits(x_gas, x_dmp, v_gas, v_dmp, nha, z, fits)
