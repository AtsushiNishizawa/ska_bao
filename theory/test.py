#!/usr/bin/env python

import numpy as np
import os, sys
import matplotlib as mpl
import pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from classy import Class
#from mpi4py import MPI
import subprocess as sp
    

if __name__ == "__main__":

    #
    # MPI initialization
    # 
    #comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()
    #size = comm.Get_size()

    #
    # get command line argument
    #
    if ( len(sys.argv)==1 ):
        print "Usage: "
        print "  $ python pk_Legendre.py [zout]"
        print "  or"
        print "  $ mpirun -np 3 python pk_Legendre.py [zout]"
        print "  for RegPTcorr full calculation (for the first time for each zout)"
        sys.exit()
    zred = float(sys.argv[1])

    #
    # set WMAP-9 cosmological paarameters
    #
    omM0 = 0.2726
    omL0 = 0.7274
    omB0 = 0.0456
    sig8 = 0.8090
    ns   = 0.9630
    h    = 0.7040

    
    #
    # derived cosmological parameters and approximation of growth rate
    #
    #Hz = np.sqrt(omM0*(1+zred)**3+omL0)
    #omMz = omM0*(1+zred)**3./(Hz)**2
    #ff = omMz**0.55

    #
    # first calculate linear P(k) and full non-linear P(k) with fitting formulae from class
    #

    #"""
    cosmo = Class()
    zlist = "0, %f"%(zred)
    cosmo.set({'output':'mPk', 
               'Omega_cdm': omM0-omB0, 
               'Omega_k':1.-omM0-omL0, 
               'h':h,
               'sigma8': sig8,
               'n_s': ns,
               'non linear':'',
               'P_k_max_h/Mpc':10., 
               'z_pk':zlist,'output_verbose':1})
    cosmo.compute()
    print "----"   
    klin = 10**np.linspace(-4,0.0,100)
    pk0 = np.array([cosmo.pk(kk,0.  ) for kk in klin])
    pkz = np.array([cosmo.pk(kk,zred) for kk in klin])
        
    # save P(k) lin to be used for TNS model
    np.savetxt("pklin_z0_wmap9.dat", zip(klin, pk0))
    
    print "----" 
    #"""

    zlist_nl = str("0, %.1f"%(zred))
    print zlist_nl
    cosmo = Class()
    print omM0-omB0, 1.-omM0-omL0, h, sig8, ns, zlist_nl
    cosmo.set({'output':'mPk',
               'Omega_cdm': omM0-omB0, 
               'Omega_k':1.-omM0-omL0, 
               'h':h,
               'sigma8': sig8,
               'n_s': ns,
               'non linear':'HALOFIT',
               'P_k_max_h/Mpc':1.,
               'z_pk':zlist_nl,
               'output_verbose':1})
    cosmo.compute()
    print "----" 
    klin = 10**np.linspace(-4,0.0,100)
    pknl  = np.array([cosmo.pk(kk,zred) for kk in klin])
    print "----" 
    if ( True ):
        fig = plt.figure()
        plt.plot(klin,pk0,"C0--", label=r"$P^{lin}(k,0)$")
        plt.plot(klin,pkz,"C1--", label=r"$P^{lin}(k,z)$")
        plt.plot(klin,pknl,"C2-", label=r"$P^{NL}(k,z)$")
        plt.xscale("log")
        plt.yscale("log")
        plt.legend()
        plt.show()
        fig.savefig("temp.png")
        sys.exit()
        
