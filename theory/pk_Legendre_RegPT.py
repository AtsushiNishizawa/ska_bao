#!/usr/bin/env python
"""
P(k,mu) modeled by using Taruya Nishimichi Saito 2010 is decomposed in the Legendre modes.

      pk_kaiser = b**2 * (pk11 + 2.d0*beta*mu**2*pk12 + beta**2*mu**4*pk22)

      pk_A_term = b**3 * ( beta * mu**2*pk_A11 + beta**2*(mu**2*pk_A12 
                + mu**4*pk_A22)  + beta**3*(mu**4*pk_A23 + mu**6*pk_A33) )

      pk_B_term = b**4 * ( mu**2 * (beta**2*pk_B12 + beta**3*pk_B13 
               + beta**4*pk_B14) + mu**4 * (beta**2*pk_B22 + beta**3*pk_B23 
               + beta**4*pk_B24) + mu**6 * (beta**3*pk_B33 + beta**4*pk_B34)
               + mu**8 * beta**4*pk_B44 )
"""
#
# code created by A.J. Nishizawa Apr. 13, 2018
#

import numpy as np
import os, sys
import matplotlib as mpl
import pylab as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from classy import Class
from mpi4py import MPI
import subprocess as sp
    
        
def get_kmax_C(klin, pkz, th=0.7):
    k_max = 0.1 #initial guess
    six_pi2 = 6.*np.pi**2.
    kmin = klin[0]
    fpk = interp1d(np.log(klin), pkz, kind="linear")
    C = 0.
    while ( C <= th ):
        C = quad(lambda lnq: np.exp(lnq)*fpk(lnq), np.log(kmin), np.log(k_max))[0]*k_max**2./six_pi2
        k_max *= 1.01
    return k_max/1.01


def get_sigv_lin(klin, pklin):
    fpk = interp1d(np.log(klin), pklin, kind="linear")
    lnq0 = np.log(klin[0])
    lnq1 = np.log(klin[-1])
    six_pi2 = 6.*np.pi**2.
    sigv_lin = quad(lambda lnq: fpk(lnq)*np.exp(lnq), lnq0, lnq1)[0]/six_pi2
    return np.sqrt(sigv_lin)


if __name__ == "__main__":

    #
    # MPI initialization
    # 
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

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
    zred = int(sys.argv[1])

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
    # fitting option (monopole/quadrupole/hexadecapole)?
    #
    fit_ps0 = True
    fit_ps2 = True
    fit_ps4 = False
    fit_verbose = True
    fit_FoG = False #True
    fit_ff  = False #True
    #model_FoG = "None"
    #model_FoG = "Lorenz"
    model_FoG = "Gauss"
    Np0 = 1#100 # number of atempts to find best fit values
    Np1 = 100 # number of atempts to find best fit values
    Np2 = 100 # number of atempts to find best fit values
    if ( fit_FoG is False ):
        Np1 = 1
    if ( fit_ff is False ):
        Np2 = 1
    
    # if we use k-dependent bias, read k-dept bias data here
    
    # set prior (in log space)
    b0_min = np.log10(0.1)
    b0_max = np.log10(5.0)

    

    #
    # derived cosmological parameters and approximation of growth rate
    #
    Hz = np.sqrt(omM0*(1+zred)**3+omL0)
    omMz = omM0*(1+zred)**3./(Hz)**2
    f0 = omMz**0.55

    #
    # data file to be fitted
    #
    #dir="./"
    dir="/work/ando/im/dir25/count/power/result6/"
    fpleg = dir+"pleg_ill_512_redshift_z%i.txt"%(zred)

    #
    # first read linear P(k) and full non-linear P(k) with fitting formulae from class
    #
    klin, pkz  = np.loadtxt("PkTree/class_wmap9_z%i_pk.dat"%(zred+1), unpack=True)
    knl , pknl = np.loadtxt("PkTree/class_wmap9_z%i_pk_nl.dat"%(zred+1), unpack=True)

    #
    # read RegPT delta-delta, delta-v, v-v power spectrum
    #
    dd = np.loadtxt("PkRegPT/pkRegPT_dd.dat")
    dv = np.loadtxt("PkRegPT/pkRegPT_dv.dat")
    vv = np.loadtxt("PkRegPT/pkRegPT_vv.dat")
    kRP = dd[:,0]
    pkRegPT_dd = dd[:,zred*4-1]
    pkRegPT_dv = dv[:,zred*4-1]
    pkRegPT_vv = vv[:,zred*4-1]
        

    #
    # RegPTcorr to get A and B terms
    #

    fAI  = "PkRegPT/pkRegPT_Aterm_I_z%.1f.dat"%(zred)
    fAII = "PkRegPT/pkRegPT_Aterm_II_z%.1f.dat"%(zred)
    fA   = "PkRegPT/pkRegPT_Aterm_z%.1f.dat"%(zred)
    fB   = "PkRegPT/pkRegPT_Bterm_z%.1f.dat"%(zred)
    

    #
    # fitting by TNS+empirical power spectra
    #

    if ( rank == 0 ):


        #
        # read simulation data from file
        #
        (kdata, pk0_data, pk2_data, pk4_data, pk0_err, pk2_err, pk4_err) = np.loadtxt(fpleg, unpack=True)

        

        #
        # read RegPT data from files
        #
        (k_A, pkA11, pkA12, pkA22, pkA23, pkA33) = np.loadtxt(fA, unpack=True)
        (k_B, pkB12, pkB13, pkB14, pkB22, pkB23, pkB24, pkB33, pkB34, pkB44) = np.genfromtxt(fB, unpack=True)
        k_0, pkg0 = knl, pknl
        k_1, pkg1 = klin, pkz
        k_K, pk11, pk12, pk22 = kRP, pkRegPT_dd, pkRegPT_dv, pkRegPT_vv
        
        k_fit_max = get_kmax_C(klin, pkz, th=0.7)
        print "kmax=", k_fit_max
        
        sigv_lin = get_sigv_lin(klin, pkz)
        print "sigv_lin=", sigv_lin
        sigv_min = 0.#np.log10(sigv_lin*1e-4)
        sigv_max = sigv_lin*5 #np.log10(sigv_lin*1e4)

        
        Ninte = 100
        mu = np.linspace(-1,1,Ninte)
        dmu = mu[1]-mu[0]
        L0 = mu*0.+1.
        L2 = 0.5*(3*mu**2-1.)
        L4 = 0.125*(35*mu**4-30*mu**2+3.)

        #kpmin = np.max(np.array((k_A.min(), k_B.min(), k_K.min())))*1.01
        #kpmax = np.min(np.array((k_A.max(), k_B.max(), k_K.max())))*0.99
        kpmin = 0.01
        kpmax = 1.5
        Nk = 100
        kp = 10.**np.linspace(np.log10(kpmin), np.log10(kpmax), Nk)

        print kpmin, kpmax, k_K.min(), k_K.max()


        fpk11 = interp1d(k_K, pk11, kind="linear")
        fpk12 = interp1d(k_K, pk12, kind="linear")
        fpk22 = interp1d(k_K, pk22, kind="linear")
        fpkg0 = interp1d(k_0, pkg0, kind="linear")
        fpkg1 = interp1d(k_1, pkg1, kind="linear")

        fpkA11 = interp1d(k_A, pkA11, kind="linear")
        fpkA12 = interp1d(k_A, pkA12, kind="linear")
        fpkA22 = interp1d(k_A, pkA22, kind="linear")
        fpkA23 = interp1d(k_A, pkA23, kind="linear")
        fpkA33 = interp1d(k_A, pkA33, kind="linear")

        fpkB12 = interp1d(k_B, pkB12, kind="linear")
        fpkB13 = interp1d(k_B, pkB13, kind="linear")
        fpkB14 = interp1d(k_B, pkB14, kind="linear")
        fpkB22 = interp1d(k_B, pkB22, kind="linear")
        fpkB23 = interp1d(k_B, pkB23, kind="linear")
        fpkB24 = interp1d(k_B, pkB24, kind="linear")
        fpkB33 = interp1d(k_B, pkB33, kind="linear")
        fpkB34 = interp1d(k_B, pkB34, kind="linear")
        fpkB44 = interp1d(k_B, pkB44, kind="linear")

        pk11 = fpk11(kp)
        pk12 = fpk12(kp)
        pk22 = fpk22(kp)
        pkg0 = fpkg0(kp)
        pkg1 = fpkg1(kp)

        pkA11=fpkA11(kp)
        pkA12=fpkA12(kp)
        pkA22=fpkA22(kp)
        pkA23=fpkA23(kp)
        pkA33=fpkA33(kp)

        pkB12=fpkB12(kp)
        pkB13=fpkB13(kp)
        pkB14=fpkB14(kp)
        pkB22=fpkB22(kp)
        pkB23=fpkB23(kp)
        pkB24=fpkB24(kp)
        pkB33=fpkB33(kp)
        pkB34=fpkB34(kp)
        pkB44=fpkB44(kp)

        chisq_min   = 1e30
        chisqL0_min = 1e30
        chisqL1_min = 1e30
        if ( fit_verbose ):
            fout   = open("chisq_fit.dat", "w")
            foutL0 = open("chisqL0_fit.dat", "w")
            foutL1 = open("chisqL1_fit.dat", "w")
        
        for ip0 in range(Np0):
            #b0 = 10.**np.random.uniform(b0_min, b0_max, 1)[0]
            b0 = 1.0
    
            for ip1 in range(Np1):                
                if ( fit_FoG ):
                    #sigv = 10.**np.random.uniform(sigv_min, sigv_max, 1)[0]
                    sigv = np.random.uniform(sigv_min, sigv_max, 1)[0]
                else:
                    sigv = sigv_lin


                for ip2 in range(Np2):
                    if ( fit_ff ):
                        ff = np.random.uniform(0,1,1)[0]
                    else:
                        ff = f0
                    beta = ff/b0

                

                    if ( model_FoG == "Gauss" ):
                        DFoG = np.exp(-(kp[:,None]*mu[None,:]*ff*sigv)**2.)
                    elif ( model_FoG == "Lorenz" ):
                        DFoG = 1./(1.+(kp[:,None]*mu[None,:]*ff*sigv)**2.)
                    else:
                        DFoG = 1.

                

                    pk_kaiser   = b0**2*(pk11[:,None]+2.0*beta*mu[None,:]**2*pk12[:,None] + beta**2*mu[None,:]**4*pk22[:,None])
                    pk_kaiserL0 = b0**2*pkg0[:,None]*(1+beta*mu[None,:]**2)**2
                    pk_kaiserL1 = b0**2*pkg1[:,None]*(1+beta*mu[None,:]**2)**2
    
                    pk_Aterm = b0**3*(beta   * mu[None,:]**2*pkA11[:,None]+
                                      beta**2*(mu[None,:]**2*pkA12[:,None]+mu[None,:]**4*pkA22[:,None])+
                                      beta**3*(mu[None,:]**4*pkA23[:,None]+mu[None,:]**6*pkA33[:,None]))

                    pk_Bterm = b0**4*(mu[None,:]**2*(beta**2*pkB12[:,None] + beta**3*pkB13[:,None] + beta**4*pkB14[:,None]) +
                                      mu[None,:]**4*(beta**2*pkB22[:,None] + beta**3*pkB23[:,None] + beta**4*pkB24[:,None]) +
                                      mu[None,:]**6*(beta**3*pkB33[:,None] + beta**4*pkB34[:,None])+
                                      mu[None,:]**8* beta**4*pkB44[:,None] )

                    psL0  = np.dot(L0, (DFoG*(pk_kaiser  +pk_Aterm+pk_Bterm)).T)*dmu *1./2.
                    psL0L0= np.dot(L0, (DFoG*(pk_kaiserL0+pk_Aterm+pk_Bterm)).T)*dmu *1./2.
                    psL0L1= np.dot(L0, (DFoG*(pk_kaiserL1                  )).T)*dmu *1./2.

                    psL2  = np.dot(L2, (DFoG*(pk_kaiser  +pk_Aterm+pk_Bterm)).T)*dmu *5./2.
                    psL2L0= np.dot(L2, (DFoG*(pk_kaiserL0+pk_Aterm+pk_Bterm)).T)*dmu *5./2.
                    psL2L1= np.dot(L2, (DFoG*(pk_kaiserL1                  )).T)*dmu *5./2.
                    
                    psL4  = np.dot(L4, (DFoG*(pk_kaiser  +pk_Aterm+pk_Bterm)).T)*dmu *9./2.
                    psL4L0= np.dot(L4, (DFoG*(pk_kaiserL0+pk_Aterm+pk_Bterm)).T)*dmu *9./2.
                    psL4L1= np.dot(L4, (DFoG*(pk_kaiserL1                  )).T)*dmu *9./2.


                    fpsL0  = interp1d(kp, psL0  , kind="linear")
                    fpsL0L0= interp1d(kp, psL0L0, kind="linear")
                    fpsL0L1= interp1d(kp, psL0L1, kind="linear")
                    fpsL2  = interp1d(kp, psL2  , kind="linear")
                    fpsL2L0= interp1d(kp, psL2L0, kind="linear")
                    fpsL2L1= interp1d(kp, psL2L1, kind="linear")
                    fpsL4  = interp1d(kp, psL4  , kind="linear")
                    fpsL4L0= interp1d(kp, psL4L0, kind="linear")
                    fpsL4L1= interp1d(kp, psL4L1, kind="linear")

                    mask = kdata <= k_fit_max

                    chisq = 0.
                    if ( fit_ps0 ):
                        chisq += np.sum( (pk0_data[mask]-fpsL0(kdata[mask]))**2 ) / np.sum( (pk0_err[mask])**2 )
                    if ( fit_ps2 ):
                        chisq += np.sum( (pk2_data[mask]-fpsL2(kdata[mask]))**2 ) / np.sum( (pk2_err[mask])**2 )
                    if ( fit_ps4 ):
                        chisq += np.sum( (pk4_data[mask]-fpsL4(kdata[mask]))**2 ) / np.sum( (pk4_err[mask])**2 )
                    if ( chisq < chisq_min ):
                        chisq_min = chisq
                        b0_best = b0
                        sigv_best = sigv
                        ff_best = ff
                    if ( fit_verbose ):
                        fout.write("%f %f %f\n"%(b0, sigv, chisq))

                    chisqL0 = 0.
                    if ( fit_ps0 ):
                        chisqL0 += np.sum( (pk0_data[mask]-fpsL0L0(kdata[mask]))**2 ) / np.sum( (pk0_err[mask])**2 )
                    if ( fit_ps2 ):
                        chisqL0 += np.sum( (pk2_data[mask]-fpsL2L0(kdata[mask]))**2 ) / np.sum( (pk2_err[mask])**2 )
                    if ( fit_ps4 ):
                        chisqL0 += np.sum( (pk4_data[mask]-fpsL4L0(kdata[mask]))**2 ) / np.sum( (pk4_err[mask])**2 )
                    if ( chisqL0 < chisqL0_min ):
                        chisqL0_min = chisqL0
                        b0_bestL0 = b0
                        sigv_bestL0 = sigv
                        ff_bestL0 = ff
                    if ( fit_verbose ):
                        foutL0.write("%f %f %f\n"%(b0, sigv, chisq))

                    chisqL1 = 0.
                    if ( fit_ps0 ):
                        chisqL1 += np.sum( (pk0_data[mask]-fpsL0L1(kdata[mask]))**2 ) / np.sum( (pk0_err[mask])**2 )
                    if ( fit_ps2 ):
                        chisqL1 += np.sum( (pk2_data[mask]-fpsL2L1(kdata[mask]))**2 ) / np.sum( (pk2_err[mask])**2 )
                    if ( fit_ps4 ):
                        chisqL1 += np.sum( (pk4_data[mask]-fpsL4L1(kdata[mask]))**2 ) / np.sum( (pk4_err[mask])**2 )
                    if ( chisqL1 < chisqL1_min ):
                        chisqL1_min = chisqL1
                        b0_bestL1 = b0
                        sigv_bestL1 = sigv
                        ff_bestL1 = ff
                    if ( fit_verbose ):
                        foutL1.write("%f %f %f\n"%(b0, sigv, chisq))


                    if ( Np0*Np1*Np2 > 784 ):
                        n = Np0*Np1*Np2/784
                    else:
                        n = 1
                    if ( ((ip0*Np1+ip1)*Np2+ip2)%(n) ==  0 ):
                        sys.stderr.write("%d %d %d %f %f %f %f %f %f %f\n"%(ip0, ip1, ip2, chisq, chisq_min, b0_best, sigv_best, sigv_lin, b0, sigv))
                    
#----- end of loop------#

        if ( fit_verbose ):
            fout.close()
            foutL0.close()
            foutL1.close()

        # for best-fit curves
        b0 = b0_best
        sigv = sigv_best
        ff = ff_best
        print "best fit bias and sigv parameters =%f , %f , %f (chisq.= %f)"%(b0_best, sigv_best, ff_best, chisq_min)
        b0L0 = b0_bestL0
        sigvL0 = sigv_bestL0
        ffL0 = ff_bestL0
        print "best fit bias and sigv parameters =%f , %f , %f (chisq.= %f)"%(b0_bestL0, sigv_bestL0, ff_bestL0, chisqL0_min)
        b0L1 = b0_bestL1
        sigvL1 = sigv_bestL1
        ffL1 = ff_bestL1
        print "best fit bias and sigv parameters =%f , %f , %f (chisq.= %f)"%(b0_bestL1, sigv_bestL1, ff_bestL1, chisqL1_min)
        print "linear growth rate at z=%i is %.3f"%(zred, f0)

        if ( model_FoG == "Gauss" ):
            DFoG   = np.exp(-(kp[:,None]*mu[None,:]*ff  *sigv  )**2.)
            DFoGL0 = np.exp(-(kp[:,None]*mu[None,:]*ffL0*sigvL0)**2.)
            DFoGL1 = np.exp(-(kp[:,None]*mu[None,:]*ffL1*sigvL1)**2.)
        elif ( model_FoG == "Lorenz" ):
            DFoG   = 1./(1.+(kp[:,None]*mu[None,:]*ff  *sigv  )**2.)
            DFoGL0 = 1./(1.+(kp[:,None]*mu[None,:]*ffL0*sigvL0)**2.)
            DFoGL1 = 1./(1.+(kp[:,None]*mu[None,:]*ffL1*sigvL1)**2.)
        else:
            DFoG   = 1.
            DFoGL0 = 1.
            DFoGL1 = 1.

        beta   = ff  /b0
        betaL0 = ffL0/b0L0
        betaL1 = ffL1/b0L1
    
        pk_kaiser   = b0  **2*(pk11[:,None]+2.0*beta*mu[None,:]**2*pk12[:,None] + beta**2*mu[None,:]**4*pk22[:,None])
        pk_kaiserL0 = b0L0**2*pkg0[:,None]*(1+betaL0*mu[None,:]**2)**2
        pk_kaiserL1 = b0L1**2*pkg1[:,None]*(1+betaL1*mu[None,:]**2)**2
    
        pk_Aterm = b0**3*(beta   * mu[None,:]**2*pkA11[:,None]+
                          beta**2*(mu[None,:]**2*pkA12[:,None]+mu[None,:]**4*pkA22[:,None])+
                          beta**3*(mu[None,:]**4*pkA23[:,None]+mu[None,:]**6*pkA33[:,None]))
        pk_AtermL0 = b0L0**3*(betaL0   * mu[None,:]**2*pkA11[:,None]+
                              betaL0**2*(mu[None,:]**2*pkA12[:,None]+mu[None,:]**4*pkA22[:,None])+
                              betaL0**3*(mu[None,:]**4*pkA23[:,None]+mu[None,:]**6*pkA33[:,None]))

        pk_Bterm = b0**4*(mu[None,:]**2*(beta**2*pkB12[:,None] + beta**3*pkB13[:,None] + beta**4*pkB14[:,None]) +
                          mu[None,:]**4*(beta**2*pkB22[:,None] + beta**3*pkB23[:,None] + beta**4*pkB24[:,None]) +
                          mu[None,:]**6*(beta**3*pkB33[:,None] + beta**4*pkB34[:,None])+
                          mu[None,:]**8* beta**4*pkB44[:,None] )
        pk_BtermL0 = b0L0**4*(mu[None,:]**2*(betaL1**2*pkB12[:,None] + betaL1**3*pkB13[:,None] + betaL1**4*pkB14[:,None]) +
                              mu[None,:]**4*(betaL1**2*pkB22[:,None] + betaL1**3*pkB23[:,None] + betaL1**4*pkB24[:,None]) +
                              mu[None,:]**6*(betaL1**3*pkB33[:,None] + betaL1**4*pkB34[:,None])+
                              mu[None,:]**8* betaL1**4*pkB44[:,None] )
        
        psL0  = np.dot(L0, (DFoG*(pk_kaiser    +pk_Aterm  +pk_Bterm  )).T)*dmu *1./2.
        psL2  = np.dot(L2, (DFoG*(pk_kaiser    +pk_Aterm  +pk_Bterm  )).T)*dmu *5./2.
        psL4  = np.dot(L4, (DFoG*(pk_kaiser    +pk_Aterm  +pk_Bterm  )).T)*dmu *9./2.
        psL0L0= np.dot(L0, (DFoGL0*(pk_kaiserL0+pk_AtermL0+pk_BtermL0)).T)*dmu *1./2.
        psL2L0= np.dot(L2, (DFoGL0*(pk_kaiserL0+pk_AtermL0+pk_BtermL0)).T)*dmu *5./2.
        psL4L0= np.dot(L4, (DFoGL0*(pk_kaiserL0+pk_AtermL0+pk_BtermL0)).T)*dmu *9./2.
        psL0L1= np.dot(L0, (DFoGL1*(pk_kaiserL1                      )).T)*dmu *1./2.
        psL2L1= np.dot(L2, (DFoGL1*(pk_kaiserL1                      )).T)*dmu *5./2.
        psL4L1= np.dot(L4, (DFoGL1*(pk_kaiserL1                      )).T)*dmu *9./2.
            

        if ( True ):
            fig = plt.figure()
            plt.plot(kp, kp*psL0  , "C0-", label=r"$P^{(s),{\rm TNS}}_0(k)$"   )
            plt.plot(kp, kp*psL0L0, "C0--", label=r"$P^{(s),{\rm nl}}_0(k)$")
            plt.plot(kp, kp*psL0L1, "C0:", label=r"$P^{(s)^{\rm lin}}_0(k)$")
            plt.errorbar(kdata, kdata*pk0_data, kdata*pk0_err, color="C0", marker="o", ls="none")
            plt.plot(kp, kp*psL2  , "C1-", label=r"$P^{(s)}_2(k)$")
            plt.plot(kp, kp*psL2L0, "C1--")
            plt.plot(kp, kp*psL2L1, "C1:" )
            plt.errorbar(kdata, kdata*pk2_data, kdata*pk2_err, color="C1", marker="s", ls="none")
            plt.plot(kp, kp*psL4  , "C2-", label=r"$P^{(s)}_4(k)$")
            plt.plot(kp, kp*psL4L0, "C2--")
            plt.plot(kp, kp*psL4L1, "C2:" )
            plt.errorbar(kdata, kdata*pk4_data, kdata*pk4_err, color="C2", marker="d", ls="none")
            plt.legend(loc="best")
            plt.xlabel(r"$k [h/{\rm Mpc}]$")
            plt.ylabel(r"$k P^{(s)}_{\ell}(k) [{\rm Mpc}^2/h^2]$")
            #plt.yscale("log")
            plt.title(r"$b0=%.2f, \sigma_v=%.2e$"%(b0_best, sigv_best))
            plt.axvline(x=k_fit_max, color="k", ls="--")
            plt.xlim(kpmin, 1.0)
            plt.ylim(ymin=-50)
            #plt.show()
            fig.savefig("pk_Legendre_bestfit_z%.1f.eps"%(zred), bbox_inches="tight", pad_inches=0.1)


        np.savetxt("test_psk_z%i.dat"%(zred), np.array([kp, psL0L1, pkg0, pkg1]).T)
