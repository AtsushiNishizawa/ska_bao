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
    # fitting option (monopole/quadrupole/hexadecapole)?
    #
    fit_ps0 = True
    fit_ps2 = True
    fit_ps4 = False
    fit_verbose = True
    fit_FoG = True
    #model_FoG = "None"
    #model_FoG = "Lorenz"
    model_FoG = "Gauss"
    Np0 = 300 # number of atempts to find best fit values
    Np1 = 300 # number of atempts to find best fit values
    if ( fit_FoG is False ):
        Np1 = 1
    
    # if we use k-dependent bias, read k-dept bias data here
    
    # set prior (in log space)
    b0_min = np.log10(0.1)
    b0_max = np.log10(5.0)

    

    #
    # derived cosmological parameters and approximation of growth rate
    #
    Hz = np.sqrt(omM0*(1+zred)**3+omL0)
    omMz = omM0*(1+zred)**3./(Hz)**2
    ff = omMz**0.55

    #
    # data file to be fitted
    #
    fpleg = "/work/ando/im/dir25/count/power/result6/pleg_ill_512_redshift_z"+str(zred)+".txt"

    #
    # first calculate linear P(k) and full non-linear P(k) with fitting formulae from class
    #

    if ( rank == 0 ):
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
        zlist_nl = "0, %.1f"%(zred)
        cosmo = Class()
        cosmo.set({'output':'mPk', 
                'Omega_cdm': omM0-omB0, 
                'Omega_k':1.-omM0-omL0, 
                'h':h,
                'sigma8': sig8,
                'n_s': ns,
                'non linear':'halofit',
                'P_k_max_h/Mpc':10.,
                'z_pk':zlist_nl,'output_verbose':1})
        cosmo.compute()
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

    sys.exit()        
    
    #
    #
    # RegPTcorr to get A and B terms
    #
    #

    fAI  = "pkRegPT_Aterm_I_z%.1f.dat"%(zred)
    fAII = "pkRegPT_Aterm_II_z%.1f.dat"%(zred)
    fA   = "pkRegPT_Aterm_z%.1f.dat"%(zred)
    fB   = "pkRegPT_Bterm_z%.1f.dat"%(zred)


    if ( rank==0 ):
        if ( os.path.exists(fAI) ):
            print "file %s exists : skipped"%(fAI)
        else:
            com = "RegPTcorr/A_term_I.exe "
            com += "%f "%(zred)
            com += "1 "
            com += "100 "
            com += "pklin_z0_wmap9.dat "
            com += "%f "%(h)
            com += "%f "%(2.723)
            com += "%f "%(ns)
            com += "%f "%(sig8)
            com += "%f "%(omM0)
            com += "%f "%(omB0)
            com += "%f "%(-1.0)
            com += fAI
            command = com.split(" ")
            print rank, com
            try:
                sp.check_output(command)
                print "A_term_I.exe done"
            except:
                print "error on A_term_I.exe"

    elif ( (size==1) or ((size>1) and (rank==1)) ):
        if ( os.path.exists(fAII) ):
            print "file %s exists : skipped"%(fAII)
        else:
            com = "RegPTcorr/A_term_II.exe "
            com += "%f "%(zred)
            com += "1 "
            com += "100 "
            com += "pklin_z0_wmap9.dat "
            com += "%f "%(h)
            com += "%f "%(2.723)
            com += "%f "%(ns)
            com += "%f "%(sig8)
            com += "%f "%(omM0)
            com += "%f "%(omB0)
            com += "%f "%(-1.0)
            com += fAII
            command = com.split(" ")
            print rank, com
            try:
                sp.check_output(command)
                print "A_term_II.exe done"
            except:
                print "error on A_term_II.exe"

    elif ( (size==1) or ((size>1) and (rank==2)) ):
        if ( os.path.exists(fB) ):
            print "file %s exists : skipped"%(fB)
        else:
            com = "RegPTcorr/B_term.exe "
            com += "%f "%(zred)
            com += "1 "
            com += "100 "
            com += "pklin_z0_wmap9.dat "
            com += "%f "%(h)
            com += "%f "%(2.723)
            com += "%f "%(ns)
            com += "%f "%(sig8)
            com += "%f "%(omM0)
            com += "%f "%(omB0)
            com += "%f "%(-1.0)
            com += fB
            command = com.split(" ")
            print rank, com
            try:
                sp.check_output(command)
                print "B_term.exe done"
            except:
                print "error on B_term.exe"

    comm.Barrier()

    if ( rank == 0 ):
        if ( os.path.exists(fA) ):
            print "file %s exists : skipped"%(fA)
        else:
            com = "RegPTcorr/sum_pkRegPT_Aterm.exe "
            com += fAI +" "
            com += fAII +" "
            com += fA
            command = com.split(" ")
            try:
                sp.check_output(command)
                print "sum_pkRegPT_Aterm.exe done"
            except:
                print "error on sum_pkRegPT_Aterm.exe"


                
    comm.Barrier()
    

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
        k_K, pkg = klin, pknl

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


        #fpk11 = interp1d(k_K, pk11, kind="linear")
        #fpk12 = interp1d(k_K, pk12, kind="linear")
        #fpk22 = interp1d(k_K, pk22, kind="linear")
        fpkg = interp1d(klin, pkg, kind="linear")

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

        #pk11 = fpk11(kp)
        #pk12 = fpk12(kp)
        #pk22 = fpk22(kp)
        pkg = fpkg(kp)

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

        chisq_min = 1e30
        if ( fit_verbose ):
            fout = open("chisq_fit.dat", "w")
        
        for ip0 in range(Np0):
            b0 = 10.**np.random.uniform(b0_min, b0_max, 1)[0]
            beta = ff/b0 # = f/b
    
            for ip1 in range(Np1):
                
                if ( fit_FoG ):
                    #sigv = 10.**np.random.uniform(sigv_min, sigv_max, 1)[0]
                    sigv = np.random.uniform(sigv_min, sigv_max, 1)[0]
                else:
                    sigv = sigv_lin

                if ( model_FoG == "Gauss" ):
                    DFoG = np.exp(-(kp[:,None]*mu[None,:]*ff*sigv)**2.)
                elif ( model_FoG == "Lorenz" ):
                    DFoG = 1./(1.+(kp[:,None]*mu[None,:]*ff*sigv)**2.)
                else:
                    DFoG = 1.

                #pk_kaiser = b0**2*(pk11[:,None]+2.0*beta*mu[None,:]**2*pk12[:,None] + beta**2*mu[None,:]**4*pk22[:,None])
                pk_kaiser= b0**2*pkg[:,None]*(1+beta*mu[None,:]**2)**2
    
                pk_Aterm = b0**3*(beta   * mu[None,:]**2*pkA11[:,None]+
                                  beta**2*(mu[None,:]**2*pkA12[:,None]+mu[None,:]**4*pkA22[:,None])+
                                  beta**3*(mu[None,:]**4*pkA23[:,None]+mu[None,:]**6*pkA33[:,None]))

                pk_Bterm = b0**4*(mu[None,:]**2*(beta**2*pkB12[:,None] + beta**3*pkB13[:,None] + beta**4*pkB14[:,None]) +
                                  mu[None,:]**4*(beta**2*pkB22[:,None] + beta**3*pkB23[:,None] + beta**4*pkB24[:,None]) +
                                  mu[None,:]**6*(beta**3*pkB33[:,None] + beta**4*pkB34[:,None])+
                                  mu[None,:]**8* beta**4*pkB44[:,None] )

                
        
                psL0 = np.dot(L0, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu
                psL2 = np.dot(L2, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu
                psL4 = np.dot(L4, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu

                fpsL0 = interp1d(kp, psL0, kind="linear")
                fpsL2 = interp1d(kp, psL2, kind="linear")
                fpsL4 = interp1d(kp, psL4, kind="linear")

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
                if ( fit_verbose ):
                    fout.write("%f %f %f\n"%(b0, sigv, chisq))

                if ( Np0*Np1 > 784 ):
                    n = Np0*Np1/784
                else:
                    n = 1
                if ( (ip0*Np1+ip1)%(n) ==  0 ):
                    sys.stderr.write("%d %d %f %f %f %f %f %f %f\n"%(ip0, ip1, chisq, chisq_min, b0_best, sigv_best, sigv_lin, b0, sigv))
                    


        if ( fit_verbose ):
            fout.close()

        # for best-fit curves
        b0 = b0_best
        sigv = sigv_best
        print "best fit bias and sigv parameters =%f , %f (chisq.= %f)="%(b0_best, sigv_best, chisq_min)

        if ( model_FoG == "Gauss" ):
            DFoG = np.exp(-(kp[:,None]*mu[None,:]*ff*sigv)**2.)
        elif ( model_FoG == "Lorenz" ):
            DFoG = 1./(1.+(kp[:,None]*mu[None,:]*ff*sigv)**2.)
        else:
            DFoG = 1.

        beta = ff/b0
    
        pk_kaiser= b0**2*pkg[:,None]*(1+beta*mu[None,:]**2)**2
    
        pk_Aterm = b0**3*(beta   * mu[None,:]**2*pkA11[:,None]+
                          beta**2*(mu[None,:]**2*pkA12[:,None]+mu[None,:]**4*pkA22[:,None])+
                          beta**3*(mu[None,:]**4*pkA23[:,None]+mu[None,:]**6*pkA33[:,None]))

        pk_Bterm = b0**4*(mu[None,:]**2*(beta**2*pkB12[:,None] + beta**3*pkB13[:,None] + beta**4*pkB14[:,None]) +
                          mu[None,:]**4*(beta**2*pkB22[:,None] + beta**3*pkB23[:,None] + beta**4*pkB24[:,None]) +
                          mu[None,:]**6*(beta**3*pkB33[:,None] + beta**4*pkB34[:,None])+
                          mu[None,:]**8* beta**4*pkB44[:,None] )
        
        psL0 = np.dot(L0, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu
        psL2 = np.dot(L2, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu
        psL4 = np.dot(L4, (DFoG*(pk_kaiser+pk_Aterm+pk_Bterm)).T)*dmu

            

        if ( True ):
            fig = plt.figure()
            plt.plot(kp, kp*psL0, "C0-", label=r"$P^{(s)}_0(k)$")
            plt.errorbar(kdata, kdata*pk0_data, kdata*pk0_err, color="C0", marker="o", ls="none")
            plt.plot(kp, kp*psL2, "C1--", label=r"$P^{(s)}_2(k)$")
            plt.errorbar(kdata, kdata*pk2_data, kdata*pk2_err, color="C1", marker="s", ls="none")
            plt.plot(kp, kp*psL4, "C2:", label=r"$P^{(s)}_4(k)$")
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
