import h5py as H
import numpy as np
import os, sys
import glob
import matplotlib as mpl
mpl.use("Agg")
import pylab as plt
from sklearn.neighbors import KDTree

snapshot = 2
dir = "/str2/ando/gad/osaka/data_box/L205_N9_wmap0/"

f = dir+"groups_%03i/sub_%03i.*.hdf5"%(snapshot, snapshot)
sub_list = glob.glob(f)
f = dir+"snapshot_%03i/snapshot_%03i.*.hdf5"%(snapshot, snapshot)
snap_list = glob.glob(f)

#print (sub_list)
#print (snap_list)


Fwrite = True
GM_C200 = None
for sub in sub_list:
    print ("reading...", sub)
    d = H.File(sub, "r")
    #----#
    if ( Fwrite ):
        for k in d.keys():
            print (k)
            for item in d[k]:
                print ("  ", item)
    Fwrite = False
    #----#
    if ( GM_C200 is None ):
        GM_C200 = np.array(d["Group/Group_M_Crit200"])
        GR_C200 = np.array(d["Group/Group_R_Crit200"])
        GPos    = np.array(d["Group/GroupPos"])
    else:
        GM_C200 = np.append(GM_C200, np.array(d["Group/Group_M_Crit200"]))
        GR_C200 = np.append(GR_C200, np.array(d["Group/Group_R_Crit200"]))
        GPos    = np.append(GPos   , np.array(d["Group/GroupPos"]), axis=0)
    

mass_comp_limit = 10.**1.25
print ("mass_comp_limit=", mass_comp_limit)

mask = np.isfinite(GM_C200)&(GM_C200>mass_comp_limit)

GM_C200 = GM_C200[mask]
GR_C200 = GR_C200[mask]
GPos    = GPos   [mask,:]

print (GM_C200.shape)
print (GM_C200.shape)
print (GPos.shape)

if ( False ):
    plt.hist(np.log10(GM_C200), bins=100)
    plt.yscale("log")
    plt.show()




# read snapshot particles
#
ic = 0
GasPos = None
Fwrite = True
for snap in snap_list:
    d = H.File(snap, "r")
    print ("reading...", snap)
    #ic += 1
    #if ( ic > 2 ):
    #    continue

    #----#
    if ( Fwrite ):
        for k in d.keys():
            print (k)
            for item in d[k]:
                print ("  ", item)
    Fwrite = False
    #----#
    if ( GasPos is None ):
        GasPos  = np.array(d["PartType0/Coordinates"])
        HIfrac  = np.array(d["PartType0/HI"]         )
        GasMass = np.array(d["PartType0/Masses"]     )
    else:
        GasPos  = np.append(GasPos,  np.array(d["PartType0/Coordinates"]), axis=0)
        HIfrac  = np.append(HIfrac,  np.array(d["PartType0/HI"]         )        )
        GasMass = np.append(GasMass, np.array(d["PartType0/Masses"]     )        )

HIMass = HIfrac*GasMass

#
GasPos = GasPos[::100,:]
HIMass = HIMass[::100]
#

print (GasPos.shape)
print (GPos.shape)
Ngrp = GPos.shape[0]

tree_T = KDTree(GasPos, leaf_size=100, metric="euclidean")
ind, dist = tree_T.query_radius(GPos, r=GR_C200, return_distance=True, count_only=False)

HIMass_indiv = np.zeros(Ngrp)
fig = plt.figure()
for i in range(ind.shape[0]):
    if ( ind[i].shape[0] >1 ):
        #print (ind[i])
        HIMass_indiv[i] = np.sum(HIMass[ind[i]])
        plt.plot(GasPos[ind[i],0], GasPos[ind[i],1], "C0,")
fig.savefig("halo2HI.png")
plt.close(fig)


fig = plt.figure()
plt.plot(GM_C200, HIMass_indiv, "C0o")
fig.savefig("halo2HI_2.png")
