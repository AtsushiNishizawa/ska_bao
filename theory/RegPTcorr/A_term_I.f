c ******************************************************* c
c                                                         c
c      Corrections to Scoccimarro (2004) formula          c
c      in redshift-space power spectrum: A term           c
c                                                         c  
c      Explicit formula suited for RegPT treatment        c  
c                                                         c  
c            Time-stamp: <2016-11-15 23:17:44 ataruya>    c
c ******************************************************* c
c
c     Note--. 
c
c     All the parameter 'ikmax' in this program must be 
c     the same. 
c
c     Exponential cutoff is taken into account 
c     (see function integ_fp). 
c
      program A_term_I
c
      implicit none
c
      integer  ik, ikmax, ik_max, ik_max1, ik_max2
      integer  itype, ibox, isigmav
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax)
      real*8  akcorr(ikmax)
      real*8  pk0corr(ikmax), pk2corr(ikmax), pk4corr(ikmax)
      real*8  pk_A1(ikmax), pk_A2(ikmax), pk_A3(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8  pi, Lbox, zred, ss, sigmav, sigma_v, growth
      real*8  h_input, Tcmb_input, n_s_input, sigma8_input
      real*8  Omega_m_input, Omega_b_input, w_de_input
      character infile*128, arg*128, outfile*128
      common /pk_data/ ak, pk, ik_max
      common /pkcorr/ akcorr, pk0corr, pk2corr, pk4corr, 
     &     pk_A1, pk_A2, pk_A3, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
      common /sigma_vz/ sigma_v
      pi = 4.d0 * atan(1.d0)
c     -------------------------------------------------
c
      write(6,*)
      write(6,*) '*** Corrections to Scoccimarro (2004) formula ***'
      write(6,*) '*** (A term) in redshift-space power spectrum ***'
      write(6,*)
      write(6,*)
      write(6,*) 'Output redshift ?'
!      read(5,*) zred
      call getarg(1, arg)
      read(arg,*) zred
      write(6,*) 'including finite-volume effect? :type (0 or 1, L_box)'
      write(6,*) ' (e.g., (0, 300) for boxsize L_box=300Mpc/h ) '
      write(6,*) ' (      (1, ???) for ignoring finite-volume effect) '
!     read(5,*) ibox, Lbox
      call getarg(2,arg)
      read(arg,*) ibox
      call getarg(3,arg)
      read(arg,*) Lbox
c
c     /////// Load (linear) matter power spectrum ///////
c
      call getarg(4,infile)
      call load_matterpower_data(infile)
c
      write(6,*) ' loading (linear) matter power spectrum, done '
c
c     //////// Set cosmological parameters ////////
c
      call getarg(5, arg); read(arg,*)h_input
      call getarg(6, arg); read(arg,*)Tcmb_input
      call getarg(7, arg); read(arg,*)n_s_input
      call getarg(8, arg); read(arg,*)sigma8_input
      call getarg(9, arg); read(arg,*)Omega_m_input
      call getarg(10,arg); read(arg,*)Omega_b_input
      call getarg(11,arg); read(arg,*)w_de_input
      call set_cosmological_params(
     +     h_input,
     +     Tcmb_input,
     +     n_s_input,
     +     sigma8_input,
     +     Omega_m_input,
     +     Omega_b_input,
     +     w_de_input)
c
c     /////// Sigma8 normalization ///////
c
      call normalization_trapez
c
      write(6,*) ' normalization done '
c
c     //////// Truncation of low-k, high-k  in linear P(k) ////////
c
      call calc_sigmav(ss)
      sigmav = ss
      write(6,'(A,1p1e14.6)') 'sigma_v=',
     &     ss * growth(zred)
!      write(6,*)
!      write(6,*) 'Use this value ? y[0], n[1]'
!      read(5,*) isigmav
!      if(isigmav.eq.1) then
!         write(6,*) 'Input sigma_v'
!         read(5,*) sigmav
!         sigmav = sigmav / growth(zred)
!      endif
c
      sigma_v = sigmav * growth(zred)
c
      call calc_running_sigmav2(growth(zred))
c
      call truncation_k(ibox, Lbox)
c
c     /////// Correction to Scoccimarro (2004) formula ///////
c
      call calc_correction(zred)
c
      write(6,*) ' 1-loop P(k) done '
c
c     /////// Summing up all contributions ///////
c
      call getarg(12,outfile)
      call calc_pkred(zred,sigmav, outfile)
c
      write(6,*) ' summing up all contributions done '
c
c     /////// Save output data ///////
c
cc      open(10,file='corr_pkred2.dat',status='unknown')
cc      open(11,file='pkcorr2_red_RegPT_I.dat',status='unknown')
c
cc      do ik =1, ik_max
cc         write(10,'(1p8e18.10)') ak(ik), pk(ik), 
cc     &        pk0EH(ik), pk0corr(ik), pk2EH(ik), 
cc     &        pk2corr(ik), pk4EH(ik), pk4corr(ik) 
cc         write(11,'(1p4e18.10)') ak(ik), pk_A1(ik), pk_A2(ik), 
cc     &        pk_A3(ik)
cc      enddo
c
cc      write(6,*)
cc      write(6,*) ' Output file: corr_pkred2.dat, pkcorr2_red.dat'
c
cc      close(10)
cc      close(11)
c
      end
c
c ******************************************************* c
c
      subroutine load_matterpower_data(infile)
c
c ******************************************************* c
c
c     input file is assumed to be the matter power spectrum data 
c     created by CAMB code. 
c
      implicit none
c
      integer ik, ikk, ikmax, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax), pi
      common /pk_data/ ak, pk, ik_max
      character(len=128) infile
      pi = 4.d0 * datan(1.d0)
c     --------------------------------------------------------
c
!      write(6,*) ' Type file name of matter P(k) data'
!      read(5,*) infile
      write(6,*)"read pk from", infile
      
      open(9, file=trim(infile)) !, status='unknown')
c
      do ik=1, ikmax
         read(9,*,END=10) ak(ik), pk(ik)
      enddo
c
 10   continue
      ik_max = ik - 1
c
      write(6,*) 'ak(1)=', ak(1)
      write(6,*) 'ak(ik_max)=', ak(ik_max)
      write(6,*) 'ik_max=', ik_max
      write(6,*)
c
      close(9)
c
      end
c
c ******************************************************* c
c
      subroutine set_cosmological_params(
     +     h_input,
     +     Tcmb_input,
     +     n_s_input,
     +     sigma8_input,
     +     Omega_m_input,
     +     Omega_b_input,
     +     w_de_input)
c
c ******************************************************* c
c
c     Set cosmological parameters for linear growth factor,
c     growth-rate parameter, and sigma8 normalization 
c
c     Note.--
c     assume a flat cosmological model (Omega_v = 1-Omega_m)
c
      implicit none
      integer iparams
      real*8 h, Tcmb, n_s, sigma8, Omega_m, Omega_b, Omega_v, w_de
      real*8 omega0
      real*8  h_input, Tcmb_input, n_s_input, sigma8_input
      real*8  Omega_m_input, Omega_b_input, w_de_input
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
      common /no_wiggle_param/ Omega_b, omega0, h, Tcmb, n_s
c     ---------------------------------------------------
c
!      write(6,*)
!      write(6,*) '*** Set cosmological parameters ***' 
!      write(6,*)
c
c     Default cosmological parameters (you can change them appropriately)
c
!      h = 0.6731d0
!      Tcmb = 2.7255d0
!      n_s = 0.9655d0
!      sigma8 = 0.8207d0
!      Omega_m = 0.3132d0
!      Omega_b = 0.049044d0
!      Omega_v = 1.d0 - Omega_m
!      w_de = -1.d0
c
! 5    write(6,*) '[1] h        =',h
!      write(6,*) '[2] T_cmb[K] =',Tcmb
!      write(6,*) '[3] n_s      =',n_s
!      write(6,*) '[4] sigma_8  =',sigma8
!      write(6,*) '[5] Omega_m  =',Omega_m
!      write(6,*) '[6] Omega_b  =',Omega_b
!      write(6,*) '[7] w_de     =',w_de
!      write(6,*) 
!      write(6,*) 'Note--. spatial curvature is forced to be flat, '
!      write(6,*) '                   i.e., Omega_v = 1 - Omega_m  '
!      write(6,*)
!      write(6,*) 'change cosmological parameter? [1-7] or n[0]'
!      read(5,*)  iparams
!      if(iparams.eq.0) goto 8
!      if(iparams.eq.1) then
!         write(6,*) 'type h'
!         read(5,*) h
!      elseif(iparams.eq.2) then
!         write(6,*) 'type Tcmb'
!         read(5,*) Tcmb
!      elseif(iparams.eq.3) then
!         write(6,*) 'type n_s'
!         read(5,*) n_s
!      elseif(iparams.eq.4) then
!         write(6,*) 'type sigma8'
!         read(5,*) sigma8
!      elseif(iparams.eq.5) then
!         write(6,*) 'type Omega_m'
!         read(5,*) Omega_m
!      elseif(iparams.eq.6) then
!         write(6,*) 'type Omega_b'
!         read(5,*) Omega_b
!      elseif(iparams.eq.7) then
!         write(6,*) 'type w'
!         read(5,*) w_de
!      else 
!         stop
!      endif
!      goto 5
c
!     8    continue

      !---AJN modified
      h       = h_input
      Tcmb    = Tcmb_input
      n_s     = n_s_input
      sigma8  = sigma8_input
      Omega_m = Omega_m_input
      Omega_b = Omega_b_input
      w_de    = w_de_input
      omega_v = 1.-omega_m
      !---AJN modified
      
      omega0 = Omega_m
c
      end
c
c ******************************************************* c
c
      subroutine normalization_trapez
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax
      parameter(ikmax=3000)
      real*8  Omega_m, Omega_v, w_de, sigma8
      real*8  ak(ikmax), pk(ikmax)
      real*8  r_th, const, x
      real*8  W_TH, sigma_a, sigma_b, pi
      common /pk_data/ ak, pk, ik_max
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
      pi = 4.d0 * atan(1.d0)
      r_th = 8.d0
c     ---------------------------------------------------
c
      x = ak(1) * r_th
      if(x.lt.1.d-3) then
         W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
      else
         W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
      endif
      sigma_a = W_TH * W_TH * pk(1) * ak(1) * ak(1)
      sigma_a = sigma_a / (2.d0 * pi * pi)
c
      const = 0.d0 
      do ik=2, ik_max
         x = ak(ik) * r_th
         if(x.lt.1.d-3) then
            W_TH = 1.d0 - x*x / 10.d0 + x**4 / 280.d0 
         else
            W_TH = 3.d0 * (sin(x) - x * cos(x))/x/x/x
         endif
         sigma_b = W_TH * W_TH * pk(ik) * ak(ik) * ak(ik) 
         sigma_b = sigma_b / (2.d0 * pi * pi)
         const = const + 
     &        (sigma_a + sigma_b) * ( ak(ik) - ak(ik-1) )/ 2.d0
         sigma_a = sigma_b
      enddo
c
      do ik=1, ik_max
         pk(ik) = sigma8 * sigma8 / const * pk(ik) 
      enddo
c
      end
c
c
c ******************************************************* c
c
      subroutine calc_running_sigmav2(dd)
c
c ******************************************************* c
c
c     computing sigma_v2 with running UV cutoff, k_UV = k
c
      implicit none
c
      integer ik, ikk, ikmax, ik_max, ik_max_k, iq
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax), pi, dd
      real*8  k, sv2
      real*8  ak_k(ikmax), sigmav2_k(ikmax)
      common /pk_data/ ak, pk, ik_max
      common /running_sigmav/ ak_k, sigmav2_k, ik_max_k
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      ik_max_k = ik_max
c
      do ik=1, ik_max-1
c
         ak_k(ik) = ak(ik)
c
         if(ik.eq.1) then
            sigmav2_k(ik) = pk(1) * ak(1) / 2.d0 / (6.d0*pi**2)
         else
            sigmav2_k(ik) = pk(1) * ak(1) / 2.d0 / (6.d0*pi**2)
            do iq=2, ik
               sigmav2_k(ik) = sigmav2_k(ik) + 
     &              ( pk(iq) + pk(iq-1) ) * ( ak(iq) - ak(iq-1) )
     &              / 2.d0 / (6.d0*pi**2)
            enddo
         endif
c
         sigmav2_k(ik) = sigmav2_k(ik) * dd**2
c
      enddo
c
      end
c
c ******************************************************* c
c
      subroutine truncation_k(ibox, Lbox)
c
c ******************************************************* c
c
      implicit none
c
      integer ik, ikk, ikmax, ik_max, ibox
      parameter(ikmax=3000)
      real*8 ak(ikmax), pk(ikmax)
      real*8 akk(ikmax), pkk(ikmax)
      real*8 kmin, kmax, Lbox, pi
      common /pk_data/ ak, pk, ik_max
      common /kmin_kmax/ kmin, kmax
      pi = 4.d0 * datan(1.d0)
c     -----------------------------------------------------
c
      kmin = 5.d-4   ! default value
c
      if(ibox.eq.0) kmin = 2.d0 * pi / Lbox
c
      kmax = 10.d0    ! default value
c
      do ik=1, ik_max
         akk(ik) = ak(ik)
         pkk(ik) = pk(ik)
      enddo
c
      ikk = 1
      do ik=1, ik_max
         if(akk(ik).ge.kmin .and. akk(ik).le.kmax
     &        .and. mod(ik,2).eq.0 ) then
            ak(ikk) = akk(ik)
            pk(ikk) = pkk(ik)
            ikk = ikk + 1
         endif
      enddo
c
      ik_max = ikk -1
c
      write(6,*) 'ak(1)=', ak(1)
      write(6,*) 'ak(ik_max)=', ak(ik_max)
      write(6,*) 'ik_max=', ik_max
      write(6,*)
c
      end
c
c ******************************************************* c
c
      subroutine find_pk(kk, pklin)
c
c ******************************************************* c
c
      implicit none
      integer ik_max, ikmax
      integer j, jmin, jmax
      parameter(ikmax=3000)
      real*8 ak(ikmax), pk(ikmax), kk, s, ds, pklin, kmin, kmax
      common /pk_data/ ak, pk, ik_max
      common /kmin_kmax/ kmin, kmax      
c     -------------------------------------------
c
      if(kk.ge.kmin .and. kk.le.kmax) then
c
         call hunt(ak, ik_max, kk, j)
c
         jmin = j - 2
         jmax = j + 2
         if(jmin.lt.1) jmin = 1
         if(jmax.ge.ik_max) jmax = ik_max
c     
         call polint(ak(jmin),pk(jmin),jmax-jmin+1,kk,s,ds)
         pklin = s
      else
         pklin = 0.d0
      endif
c      
      end
c
c ******************************************************* c
c
      subroutine find_running_sigmav2(kk, running_v2)
c
c ******************************************************* c
c
c     finding sigma_v2 with running UV cutoff, k_UV = k/2
c
      implicit none
      integer ik_max_k, ikmax
      integer j, jmin, jmax
      parameter(ikmax=3000)
      real*8 ak_k(ikmax), sigmav2_k(ikmax), k, kk, s, ds, running_v2
      common /running_sigmav/ ak_k, sigmav2_k, ik_max_k
c     -------------------------------------------
c
      do j=1, ik_max_k
         if(ak_k(j).ge.(kk/2.d0)) goto 1
      enddo
c
      if(kk/2.d0 .lt. ak_k(1)) s = sigmav2_k(1) 
      if(kk/2.d0 .gt. ak_k(ik_max_k)) s = sigmav2_k(ik_max_k)
      goto 2
c
cc      !!! for some reasons, hunt does not work 
cc      call hunt(ak_k, ik_max_k, 0.5d0*kk, j) 
c
 1    jmin = j - 2
      jmax = j + 2
      if(jmin.lt.1) jmin = 1
      if(jmax.ge.ik_max_k) jmax = ik_max_k
c
      call polint(ak_k(jmin),sigmav2_k(jmin),jmax-jmin+1,0.5d0*kk,s,ds)
 2    running_v2 = s
c
      end
c
c ******************************************************* c
c
      subroutine calc_correction(zred)
c
c ******************************************************* c
c
      implicit none
      integer ik, ik_max, ikmax, isub, ip
      integer ix, ixmax
      parameter(ikmax=3000)
      parameter(ixmax=300)
      real*8  pi
      real*8  ak(ikmax), pk(ikmax)
      real*8  pk_A111(ikmax), pk_A121(ikmax), pk_A221(ikmax)
      real*8  pk_A212(ikmax), pk_A222(ikmax), pk_A322(ikmax)
      real*8  pk_tA111(ikmax), pk_tA121(ikmax), pk_tA221(ikmax)
      real*8  pk_tA212(ikmax), pk_tA222(ikmax), pk_tA322(ikmax)
      real*8  kmin, kmax, xmin, xmax, mumin, mumax
      real*8  k, ww(ixmax), xx(ixmax)
      real*8  int_A(6), int_tA(6), zred, zz
      common /redshift/ zz 
      common /pk_data/ ak, pk, ik_max
      common /wave_number/  k, xmin, xmax
      common /corr_pk/ pk_A111, pk_A121, pk_A221, pk_A212, pk_A222, 
     &     pk_A322, pk_tA111, pk_tA121, pk_tA221, pk_tA212,
     &     pk_tA222, pk_tA322
      pi = 4.d0 * datan(1.d0)
c     ---------------------------------------------------
c
      zz = zred
      kmin = ak(1)
      kmax = ak(ik_max) 
c
      do 10 ik=1, ik_max
c
         k = ak(ik)
         pk_A111(ik)=0.d0 
         pk_A121(ik)=0.d0
         pk_A221(ik)=0.d0
         pk_A212(ik)=0.d0
         pk_A222(ik)=0.d0
         pk_A322(ik)=0.d0
         pk_tA111(ik)=0.d0 
         pk_tA121(ik)=0.d0
         pk_tA221(ik)=0.d0
         pk_tA212(ik)=0.d0
         pk_tA222(ik)=0.d0
         pk_tA322(ik)=0.d0
c
         xmin = kmin / k
         xmax = kmax / k
c
c     ////// Gauss-Legendre integration //////  c
c
cc         if(k.lt.0.2) isub =200 
cc         if(k.ge.0.2) isub =0 
         isub = 0
c
         call gauleg(log(xmin),log(xmax),xx,ww,ixmax-isub)
c
         do ix=1, ixmax-isub
            xx(ix)= dexp(xx(ix))
            call integ_fp(xx(ix), int_A, int_tA)
            pk_A111(ik) = pk_A111(ik)+ww(ix)*int_A(1)
            pk_A121(ik) = pk_A121(ik)+ww(ix)*int_A(2)
            pk_A221(ik) = pk_A221(ik)+ww(ix)*int_A(3)
            pk_A212(ik) = pk_A212(ik)+ww(ix)*int_A(4)
            pk_A222(ik) = pk_A222(ik)+ww(ix)*int_A(5)
            pk_A322(ik) = pk_A322(ik)+ww(ix)*int_A(6)
            pk_tA111(ik) = pk_tA111(ik)+ww(ix)*int_tA(1)
            pk_tA121(ik) = pk_tA121(ik)+ww(ix)*int_tA(2)
            pk_tA221(ik) = pk_tA221(ik)+ww(ix)*int_tA(3)
            pk_tA212(ik) = pk_tA212(ik)+ww(ix)*int_tA(4)
            pk_tA222(ik) = pk_tA222(ik)+ww(ix)*int_tA(5)
            pk_tA322(ik) = pk_tA322(ik)+ww(ix)*int_tA(6)
         enddo
c
         pk_A111(ik) = 2.d0 * pk_A111(ik) * k**3 / (2.*pi)**2
         pk_A121(ik) = 2.d0 * pk_A121(ik) * k**3 / (2.*pi)**2
         pk_A221(ik) = 2.d0 * pk_A221(ik) * k**3 / (2.*pi)**2
         pk_A212(ik) = 2.d0 * pk_A212(ik) * k**3 / (2.*pi)**2
         pk_A222(ik) = 2.d0 * pk_A222(ik) * k**3 / (2.*pi)**2
         pk_A322(ik) = 2.d0 * pk_A322(ik) * k**3 / (2.*pi)**2
         pk_tA111(ik) = 2.d0 * pk_tA111(ik) * k**3 / (2.*pi)**2
         pk_tA121(ik) = 2.d0 * pk_tA121(ik) * k**3 / (2.*pi)**2
         pk_tA221(ik) = 2.d0 * pk_tA221(ik) * k**3 / (2.*pi)**2
         pk_tA212(ik) = 2.d0 * pk_tA212(ik) * k**3 / (2.*pi)**2
         pk_tA222(ik) = 2.d0 * pk_tA222(ik) * k**3 / (2.*pi)**2
         pk_tA322(ik) = 2.d0 * pk_tA322(ik) * k**3 / (2.*pi)**2
c
c
      write(6,'(i3,1p13e18.10)') ik,k, pk_A111(ik), pk_A121(ik),
     &        pk_A221(ik), pk_A212(ik), pk_A222(ik), pk_A322(ik), 
     &        pk_tA111(ik), pk_tA121(ik), pk_tA221(ik), pk_tA212(ik),
     &        pk_tA222(ik), pk_tA322(ik)
c
 10   continue
c
      end
c
c ******************************************************* c
c
      subroutine integ_fp(x, int_A, int_tA)
c
c ******************************************************* c
c
c     int_A(1):  A111,     int_tA(1):  tA111
c     int_A(2):  A121,     int_tA(2):  tA121
c     int_A(3):  A221,     int_tA(3):  tA221
c     int_A(4):  A212,     int_tA(4):  tA212
c     int_A(5):  A222,     int_tA(5):  tA222
c     int_A(6):  A322,     int_tA(6):  tA322
c
      implicit none
      integer ip, imu, imu_max
      parameter(imu_max=10)
      real*8  int_A(6), int_tA(6), xmin, xmax, mumin, mumax
      real*8  integ_fp_A(6), integ_fp_tA(6)
      real*8  k, x, wmu(imu_max), mmu(imu_max)
      real*8  sigma_v, exp_fact, p, kp
      common /wave_number/  k, xmin, xmax
      common /sigma_vz/ sigma_v
c     ------------------------------------------  c
c
      mumin = max(-1.0, (1.+x**2-xmax**2)/2./x)
      mumax = min( 1.0, (1.+x**2-xmin**2)/2./x)
c
      if(x.ge.0.5d0) mumax= 0.5d0/x
c
      call gauleg(mumin, mumax, mmu, wmu, imu_max)
c
      do ip=1, 6
         int_A(ip) = 0.d0
         int_tA(ip) = 0.d0
      enddo
c
      do imu=1, imu_max
c
         p = k * x
         exp_fact =  k*k*(1.d0 + x*x - x*mmu(imu)) * sigma_v**2
ccc         if(exp_fact.ge.100) goto 10
c
         call fp_A(x, mmu(imu), integ_fp_A, integ_fp_tA)
c
         do ip=1, 6
            int_A(ip) = int_A(ip) + wmu(imu) * integ_fp_A(ip)
            int_tA(ip) = int_tA(ip) + wmu(imu) * integ_fp_tA(ip)
         enddo
c
      enddo
c
 10   continue
c
      end
c
c ******************************************************* c
c
      subroutine  fp_A(x, mu, integ_fp_A, integ_fp_tA)
c
c ******************************************************* c
c
c     integ_fp_A(1):  A111,    integ_fp_tA(1):  tA111
c     integ_fp_A(2):  A121,    integ_fp_tA(2):  tA121
c     integ_fp_A(3):  A221,    integ_fp_tA(3):  tA221
c     integ_fp_A(4):  A212,    integ_fp_tA(4):  tA212
c     integ_fp_A(5):  A222,    integ_fp_tA(5):  tA222
c     integ_fp_A(6):  A322,    integ_fp_tA(6):  tA322
c
      implicit none
      real*8 mu, k, x
      real*8 xmax, xmin, integ_fp_A(6), integ_fp_tA(6)
      real*8 p, kp, zred
      real*8 b211A, b221A, b212A, b222A, b211tA, b221tA, b212tA, b222tA
      common /redshift/ zred
      common /wave_number/  k, xmin, xmax
c     --------------------------------------- c
c
      p = k * x
      kp = k * dsqrt(1.d0 + x*x - 2.d0*x*mu)
c
      integ_fp_A(1) = x * mu 
      integ_fp_A(2) = - x*x * (3.*x*mu-2.) * (mu*mu-1.) 
     &     / (1.+x*x-2.*mu*x) /2.
      integ_fp_A(3) = x*( 2.*mu + x*(2.-6.*mu*mu) + 
     &     x*x*mu*(-3.+5.*mu*mu) ) / (1.+x*x-2.*mu*x) / 2.
      integ_fp_A(4) = integ_fp_A(1)
      integ_fp_A(5) = integ_fp_A(2)
      integ_fp_A(6) = integ_fp_A(3)
c
      integ_fp_tA(1) = - x*x * (x*mu-1.) / (1.+x*x-2.*x*mu) 
      integ_fp_tA(2) = x*x * (3.*x*mu-1.) * (mu*mu-1.)
     &     / (1.+x*x-2.*x*mu) /2.
      integ_fp_tA(3) = x*x * (-1. + 3.*x*mu + 3.*mu*mu - 5.*x*mu**3 ) 
     &     / (1.+x*x-2.*x*mu) /2.
      integ_fp_tA(4) = integ_fp_tA(1)
      integ_fp_tA(5) = integ_fp_tA(2)
      integ_fp_tA(6) = integ_fp_tA(3)
c
      call Bispec(zred, k, p, kp, b211A, b221A, b212A, b222A, 
     &     b211tA, b221tA, b212tA, b222tA)
c
      integ_fp_A(1) = integ_fp_A(1) * b211A * x
      integ_fp_A(2) = integ_fp_A(2) * b221A * x
      integ_fp_A(3) = integ_fp_A(3) * b221A * x
      integ_fp_A(4) = integ_fp_A(4) * b212A * x
      integ_fp_A(5) = integ_fp_A(5) * b222A * x
      integ_fp_A(6) = integ_fp_A(6) * b222A * x
c
      integ_fp_tA(1) = integ_fp_tA(1) * b211tA * x
      integ_fp_tA(2) = integ_fp_tA(2) * b221tA * x
      integ_fp_tA(3) = integ_fp_tA(3) * b221tA * x
      integ_fp_tA(4) = integ_fp_tA(4) * b212tA * x 
      integ_fp_tA(5) = integ_fp_tA(5) * b222tA * x
      integ_fp_tA(6) = integ_fp_tA(6) * b222tA * x
c
      end
c
c ******************************************************* c
c
      subroutine calc_pkred(zred, sigmav, outfile)
c
c ******************************************************* c
c
c     Summing up all contributions to redshift P(k) in PT
c     and calculating monopole, quadrupole and hexadecapole
c     spectra
c
      implicit none
      integer ik, ikmax, ik_max, ik_max1, ik_max2
      parameter(ikmax=3000)
      real*8  zred, growth, f, dd2, ff
      real*8  Omega_m, Omega_v, w_de, sigma8
      real*8  Omega_b, h, Tcmb, n_s
      real*8  ak(ikmax), pk(ikmax), pk_EH(ikmax)
      real*8  akcorr(ikmax)
      real*8  pk_A111(ikmax), pk_A121(ikmax), pk_A221(ikmax)
      real*8  pk_A212(ikmax), pk_A222(ikmax), pk_A322(ikmax)
      real*8  pk_tA111(ikmax), pk_tA121(ikmax), pk_tA221(ikmax)
      real*8  pk_tA212(ikmax), pk_tA222(ikmax), pk_tA322(ikmax)
c
      real*8  pk0corr(ikmax), pk2corr(ikmax), pk4corr(ikmax)
      real*8  pk_A1(ikmax), pk_A2(ikmax), pk_A3(ikmax)
      real*8  akEH(ikmax), pk0EH(ikmax), pk2EH(ikmax), pk4EH(ikmax)
      real*8  fact, alpha, sigmav, omega0
      character outfile*128
      common /alpha_param/ alpha
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
      common /no_wiggle_param/ Omega_b, omega0, h, Tcmb, n_s
      common /pk_data/ ak, pk, ik_max
      common /corr_pk/ pk_A111, pk_A121, pk_A221, pk_A212, pk_A222, 
     &     pk_A322, pk_tA111, pk_tA121, pk_tA221, pk_tA212,
     &     pk_tA222, pk_tA322
      common /pkcorr/ akcorr, pk0corr, pk2corr, pk4corr, 
     &     pk_A1, pk_A2, pk_A3, ik_max1
      common /pkEH/ akEH, pk0EH, pk2EH, pk4EH, ik_max2
c     -------------------------------------------------------- c
      ik_max1 = ik_max
      ik_max2 = ik_max
c
c     growth factor and its logarithmic derivative
c     (assuming flat universe)
      dd2 = growth(zred)**2
      ff = f(zred)
c
      write(6,'(A,1p3e18.10)') 'd, f, sigmav=', 
     &     dd2, ff, sigmav*dsqrt(dd2)
c
c     ////// corrections to power spectrum  ////// c
c
!      open(11,file='pkRegPT_Aterm_I.dat',status='unknown')
      open(11,file=trim(outfile),status='unknown')
c
      do ik=1, ik_max1
c
         alpha = (ak(ik)*ff*sigmav)**2 * dd2
c
         akcorr(ik) = ak(ik)
c
         pk_A1(ik) = ff * (pk_A111(ik) + pk_tA111(ik) ) +
     &        ff*ff * (pk_A121(ik) + pk_tA121(ik) ) 
         pk_A2(ik) = ff*ff * ( pk_A221(ik) + pk_A212(ik) + 
     &        pk_tA221(ik) + pk_tA212(ik) ) + 
     &        ff*ff*ff * ( pk_A222(ik) + pk_tA222(ik) ) 
         pk_A3(ik) = ff*ff*ff * ( pk_A322(ik) + pk_tA322(ik) )
c
         pk0corr(ik) = fact(1,0) * pk_A1(ik)
     &        + fact(2,0) * pk_A2(ik) + fact(3,0) * pk_A3(ik)
c
         pk2corr(ik) = fact(1,2) * pk_A1(ik)
     &        + fact(2,2) * pk_A2(ik) + fact(3,2) * pk_A3(ik)
c
         pk4corr(ik) = fact(1,4) * pk_A1(ik)
     &        + fact(2,4) * pk_A2(ik) + fact(3,4) * pk_A3(ik)
c
         write(11,'(1p6e18.10)') akcorr(ik), 
     &        pk_A111(ik) + pk_tA111(ik), ! mu^2 * f
     &        pk_A121(ik) + pk_tA121(ik), ! mu^2 * f*f
     &        pk_A221(ik) + pk_A212(ik) + pk_tA221(ik) + pk_tA212(ik), ! mu^4* f*f*f
     &        pk_A222(ik) + pk_tA222(ik), ! mu^4 * f*f*f
     &        pk_A322(ik) + pk_tA322(ik) ! mu^6 * f*f*f
c
      enddo
c
      close(11)
      write(6,*) ' Output file: ', trim(outfile)
c
c     ////// no-wiggle linear power spectrum  ////// c
c
      call no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
      ff = dabs(ff)
c
      do ik=1, ik_max2
         akEH(ik) = ak(ik)
         pk(ik) = pk(ik) * dd2 
         pk0EH(ik) = (1.d0+2.d0/3.d0*ff+1.d0/5.d0*ff*ff)* pk_EH(ik)*dd2
         pk2EH(ik) = ( 4.d0/3.d0*ff + 4.d0/7.d0*ff*ff )* pk_EH(ik)*dd2
         pk4EH(ik) = 8.d0/35.d0*ff*ff * pk_EH(ik)*dd2
      enddo
c
      end
c
c ******************************************************* c
c
      subroutine calc_sigmav(sigmav)
c
c ******************************************************* c
c
      implicit none
      integer ikmax, ik, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk(ikmax), sigmav, pi
      common /pk_data/ ak, pk, ik_max
c     -------------------------------------------------
      pi = 4.d0 * datan(1.d0)
c
      sigmav = 0.d0
c
      do ik=1, ik_max-1
         sigmav = sigmav + (pk(ik+1) + pk(ik)) * 
     &        (ak(ik+1)-ak(ik)) /2.d0
      enddo
c
      sigmav = dsqrt(sigmav /(2*pi*pi) /3.d0)
c
      end
c
c ************************************************ c
c
      function fact(n, l)
c
c ************************************************ c
c
c     (2l+1)/2 * integ dmu  mu^(2n) * exp(-alpha*mu^2) * P_l(mu)
c
      implicit none
      integer n, l
      real*8 fact, nn, alpha, gamhalf, gammp
      common /alpha_param/ alpha
c     ---------------------------- c
      nn = dble(n)
c
      if(alpha.gt.0.05) then
c
         if(l.eq.0) then
            fact = gamhalf(n) * gammp(0.5+nn,alpha)
            fact = fact / alpha**(nn+0.5) / 4.d0
         elseif(l.eq.2) then
            fact = alpha * gamhalf(n) * gammp(0.5+nn,alpha) 
     &           - 3.d0 * gamhalf(n+1) * gammp(1.5+nn,alpha) 
            fact = fact / alpha**(nn+1.5) * (-5.d0/8.d0)
         elseif(l.eq.4) then
            fact = 12.*gamhalf(n)*gammp(0.5+nn,alpha)/alpha**(n+0.5)  
     &           -120.*gamhalf(n+1)*gammp(1.5+nn,alpha)/alpha**(n+1.5)
     &           +140.*gamhalf(n+2)*gammp(2.5+nn,alpha)/alpha**(n+2.5)
            fact = fact * 9./128.
         endif
c
      else
c
         if(l.eq.0) then
            fact = 1./(2.+4.*nn) - alpha/(6.+4.*nn) 
     &           + alpha**2/(20.+8.*nn) 
         elseif(l.eq.2) then
            fact = nn/(3.+8.*nn+4.*nn**2) 
     &           - (nn+1.)*alpha/(15.+16.*nn+4.*nn**2) 
     &           + (nn+2.)*alpha**2/(70.+48.*nn+8.*nn**2)
            fact = fact * 5.d0
         elseif(l.eq.4) then
            fact = dble(n*(n-1))/dble(15+46*n+36*n**2+8*n**3) 
     &           - dble(n*(n+1))/dble(105+142*n+60*n**2+8*n**3)*alpha
     &           + dble((n+1)*(n+2))/dble(315+286*n+84*n**2+8*n**3)
     &           *alpha**2/2.d0
            fact = fact * 18.d0
         endif
c         
      endif
c
      fact = fact * (1.d0 + (-1.d0)**(2.*nn)) 
c
      end
c
c ************************************************ c
c
      function gamhalf(n)
c
c ************************************************ c
c
c     Gamma(n+1/2) up to n=6
c
      integer n
      real*8 gamhalf, pi
      pi = 4.d0 *datan(1.d0)
c
      if(n.eq.0) gamhalf = 1.d0
      if(n.eq.1) gamhalf = 0.5d0
      if(n.eq.2) gamhalf = 0.75d0
      if(n.eq.3) gamhalf = 1.875d0
      if(n.eq.4) gamhalf = 6.5625d0
      if(n.eq.5) gamhalf = 29.53125d0
      if(n.eq.6) gamhalf = 162.421875d0
      if(n.eq.7) gamhalf = 1055.7421875d0
      if(n.eq.8) gamhalf = 7918.06640625d0
c
      gamhalf = gamhalf * dsqrt(pi)
c
      end
c
c ************************************************ c
c
      subroutine no_wiggle_Pk(sigma8, ik_max, ak, pk_EH)
c
c ************************************************ c
c
      implicit none
      integer ik, ikmax, ik_max
      parameter(ikmax=3000)
      real*8  ak(ikmax), pk_EH(ikmax)
      real*8  sigma8
      real*8  const, r_th, ks, kmin, kmax, s1, s2
      real*8  Pk_lin_EH
      real*8  var
      external var, midpnt, midinf
      common /R_tophat/ r_th
c     ----------------------------------
c
c ///// normalization by sigma8 ///// c
c
      r_th = 8.d0
      ks = 1.d0 / r_th
      kmin = 0.00001 * ks
      kmax = 1.d3
c
      call qromo(var, kmin, ks, s1, midpnt)
      call qromo(var, ks, kmax, s2, midinf)
c
      const = sigma8**2 / (s1 + s2)
c
      do ik = 1, ik_max
         pk_EH(ik) = const * pk_lin_EH(ak(ik))
      enddo
c
      end
c ************************************************ c
c
      function var(k) 
c
c ************************************************ c
c
      implicit none
      real*8  var, k, x, w_th, Pk_lin_EH
      real*8  pi, r_th
      common /R_tophat/ r_th
c     -----------------------------------
c
      pi = 4.d0 * datan(1.d0)
c
      x = k * r_th
      w_th = 3.d0 * (sin(x)-x*cos(x)) / x**3 
      var = k**2 * w_th**2 * Pk_lin_EH(k) / (2.d0*pi**2)
c
      end
c
c ************************************************ c
c
      function Pk_lin_EH(k)
c
c ************************************************ c
c
c     compute un-normalized linear P(k) 
c     based on eq.(29) of Eisenstein & Hu (1998)
c     (no-wiggle approximation)  
c
      implicit none
      real*8 k, ss, alpha_gam, theta_cmb
      real*8 gamma_eff, q, L0, C0, T_EH, Pk_lin_EH
      real*8 omegab, omega0, h, Tcmb, n_s
      common /no_wiggle_param/ omegab, omega0, h, Tcmb, n_s
c     -----------------------------------------
c
c ///// fitting formula for no-wiggle P(k) (Eq.[29] of EH98)
c
      ss = 44.5 * h * dlog( 9.83 / (omega0*h*h) ) / 
     &     dsqrt( 1.d0 + 10.d0 * (omegab*h*h)**0.75 )
      alpha_gam = 1.d0 
     &     - 0.328 * dlog( 431. * omega0*h*h ) * omegab/omega0
     &     + 0.38 * dlog( 22.3 * omega0*h*h ) * (omegab/omega0)**2
      theta_cmb = Tcmb / 2.70 
      gamma_eff = omega0 * h * 
     &     ( alpha_gam + (1.d0 - alpha_gam) / (1.d0 + (0.43*k*ss)**4) )
c
      q = k * theta_cmb**2 / gamma_eff
      L0 = dlog( 2.d0 * dexp(1.d0) + 1.8 * q ) 
      C0 = 14.2 + 731.d0 / ( 1.d0 + 62.5 * q )
c
      T_EH = L0 / (L0 + C0*q*q )
c
      Pk_lin_EH = k ** n_s * T_EH**2 
c
      end
c
c ******************************************************* c
c
      function growth(zred)
c
c ******************************************************* c
c
c     Linear growth factor for flat cosmology
c
      implicit none
      real*8 growth, zred, Omega_m, Omega_v, w_de, sigma8, a, b, c
      real*8 zred1, zi, zz, gz, g0
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
c     --------------------------------------------
c
      a = -1.d0 / (3.d0 * w_de)
      b = (w_de - 1.d0)/ (2.d0 * w_de)
      c = 1.d0 - 5.d0 / (6.d0 * w_de)
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3.d0*w_de) 
c
      call HYGFX(a,b,c,zi, g0)  
      call HYGFX(a,b,c,zz, gz)  
c
      growth = (gz/g0) / zred1
c
      end
c
c ******************************************************* c
c
      function f(zred)
c
c ******************************************************* c
c
c     d lnD_+(z) / d ln a,   as function of redshift 
c
      implicit none
c
      real*8  f, zred, zred1, Omega_m, Omega_v, w_de, sigma8
      real*8  a, b, c, zi, zz, g1, g2
      common /cosmological_param/ Omega_m, Omega_v, w_de, sigma8
c     ---------------------------------------------------
c
      zred1 = 1.d0 + zred
      zi = - Omega_v / Omega_m
      zz = zi * zred1**(3.d0*w_de) 
c
      a = 1.d0 - 1.d0 / (3.d0 * w_de)
      b = 1.5d0 - 1.d0 / (2.d0 * w_de)
      c = 2.d0 - 5.d0 / (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g1)  
c
      a = - 1.d0 / (3.d0 * w_de)
      b = (-1.d0 + w_de) / (2.d0 * w_de)
      c = 1.d0 - 5.d0 / (6.d0 * w_de)
      call HYGFX(a,b,c,zz, g2)  
c
      f = 1.d0 + 3.d0 * (-1.d0 + w_de) / (6.d0*w_de -5.d0) 
     &     * zz * g1 / g2
c
      end
c
c ******************************************************* c
c
        SUBROUTINE HYGFX(A,B,C,X,HF)
c
c ******************************************************* c
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMAX for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMAX(C,GC)
           CALL GAMMAX(C-A-B,GCAB)
           CALL GAMMAX(C-A,GCA)
           CALL GAMMAX(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMAX(C,G1)
           CALL GAMMAX(1.0D0+A/2.0-B,G2)
           CALL GAMMAX(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(A+M,GAM)
              CALL GAMMAX(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMAX(A,GA)
              CALL GAMMAX(B,GB)
              CALL GAMMAX(C,GC)
              CALL GAMMAX(C-A,GCA)
              CALL GAMMAX(C-B,GCB)
              CALL GAMMAX(C-A-B,GCAB)
              CALL GAMMAX(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE GAMMAX(X,GA)
c
c ******************************************************* c
C
C       ==================================================
C       Purpose: Compute gamma function â(x)
C       Input :  x  --- Argument of â(x)
C                       ( x is not equal to 0,-1,-2,úúú)
C       Output:  GA --- â(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END
c
c ******************************************************* c
c
        SUBROUTINE PSI(X,PS)
c
c ******************************************************* c
C
C       ======================================
C       Purpose: Compute Psi function
C       Input :  x  --- Argument of psi(x)
C       Output:  PS --- psi(x)
C       ======================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XA=DABS(X)
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        S=0.0D0
        IF (X.EQ.INT(X).AND.X.LE.0.0) THEN
           PS=1.0D+300
           RETURN
        ELSE IF (XA.EQ.INT(XA)) THEN
           N=XA
           DO 10 K=1 ,N-1
10            S=S+1.0D0/K
           PS=-EL+S
        ELSE IF (XA+.5.EQ.INT(XA+.5)) THEN
           N=XA-.5
           DO 20 K=1,N
20            S=S+1.0/(2.0D0*K-1.0D0)
           PS=-EL+2.0D0*S-1.386294361119891D0
        ELSE
           IF (XA.LT.10.0) THEN
              N=10-INT(XA)
              DO 30 K=0,N-1
30               S=S+1.0D0/(XA+K)
              XA=XA+N
           ENDIF
           X2=1.0D0/(XA*XA)
           A1=-.8333333333333D-01
           A2=.83333333333333333D-02
           A3=-.39682539682539683D-02
           A4=.41666666666666667D-02
           A5=-.75757575757575758D-02
           A6=.21092796092796093D-01
           A7=-.83333333333333333D-01
           A8=.4432598039215686D0
           PS=DLOG(XA)-.5D0/XA+X2*(((((((A8*X2+A7)*X2+
     &        A6)*X2+A5)*X2+A4)*X2+A3)*X2+A2)*X2+A1)
           PS=PS-S
        ENDIF
        IF (X.LT.0.0) PS=PS-PI*DCOS(PI*X)/DSIN(PI*X)-1.0D0/X
        RETURN
        END
c
c ************************************************************
c
      SUBROUTINE qromo(func,a,b,ss,choose)
c
c ************************************************************
c
      implicit none
      INTEGER JMAX,JMAXP,K,KM
      REAL*8 a,b,func,ss,EPS
      EXTERNAL func,choose
      PARAMETER (EPS=1.d-6, JMAX=14, JMAXP=JMAX+1, K=5, KM=K-1)
CU    USES polint
      INTEGER j
      REAL*8 dss,h(JMAXP),s(JMAXP)
      h(1)=1.
      do 11 j=1,JMAX
        call choose(func,a,b,s(j),j)
        if (j.ge.K) then
          call polint(h(j-KM),s(j-KM),K,0.d0,ss,dss)
          if (abs(dss).le.EPS*abs(ss)) return
        endif
        s(j+1)=s(j)
        h(j+1)=h(j)/9.
11    continue
      stop 'too many steps in qromo'
      END
c
c ************************************************************
c
      SUBROUTINE polint(xa,ya,n,x,y,dy)
c
c ************************************************************
c
      implicit none
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)stop 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C
c ************************************************ c
c
      SUBROUTINE hunt(xx,n,x,jlo)
c
c ************************************************ c
c
      implicit none
      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
c
C ************************************************************
C
      SUBROUTINE midpnt(func,a,b,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 ddel,del,sum,tnm,x
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
C
C ************************************************************
C
      SUBROUTINE midinf(funk,aa,bb,s,n)
C
C ************************************************************
C
      implicit none
      INTEGER n
      REAL*8 aa,bb,s,funk
      EXTERNAL funk
      INTEGER it,j
      REAL*8 a,b,ddel,del,sum,tnm,func,x
      func(x)=funk(1./x)/x**2
      b=1./aa
      a=1./bb
      if (n.eq.1) then
        s=(b-a)*func(0.5*(a+b))
      else
        it=3**(n-2)
        tnm=it
        del=(b-a)/(3.*tnm)
        ddel=del+del
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+ddel
          sum=sum+func(x)
          x=x+del
11      continue
        s=(s+(b-a)*sum/tnm)/3.
      endif
      return
      END
c
c ************************************************************
c
      SUBROUTINE gauleg(x1,x2,x,w,n)
c
c ************************************************************
c
      INTEGER n
      REAL*8 x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)
      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
        z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software z!0(0.
c ************************************************************
c
      FUNCTION gammp(a,x)
c
c ************************************************************
c
      REAL*8 a,gammp,x
CU    USES gcf,gser
      REAL*8 gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)stop 'bad arguments in gammp'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammp=gamser
      else
        call gcf(gammcf,a,x,gln)
        gammp=1.-gammcf
      endif
      return
      END
c
c
c ************************************************************
c
      SUBROUTINE gser(gamser,a,x,gln)
c
c ************************************************************
c
      INTEGER ITMAX
      REAL*8 a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.d-7)
CU    USES gammln
      INTEGER n
      REAL*8 ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)stop 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      stop 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END
c
c ************************************************************

      SUBROUTINE gcf(gammcf,a,x,gln)
c
c ************************************************************
c
      INTEGER ITMAX
      REAL*8 a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.d-7,FPMIN=1.d-30)
CU    USES gammln
      INTEGER i
      REAL*8 an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      stop 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END
c
c ************************************************************
c
      FUNCTION gammln(xx)
c
c ************************************************************
c
      REAL*8 gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c
