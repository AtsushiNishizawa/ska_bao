c ******************************************************* c
c
      subroutine Bispec(zred, k, p, kp, b211A, b221A, b212A, b222A, 
     &     b211tA, b221tA, b212tA, b222tA)
c                                                         c
c            Time-stamp: <2012-04-18 16:08:22 ataruya>    c
c ******************************************************* c
c
c     bispectrum calculation with RegPT
c
c      b211A = Bispec(2, 1, 1, p, kp, k, zred)
c      b221A = Bispec(2, 2, 1, p, kp, k, zred)
c      b212A = Bispec(2, 1, 2, p, kp, k, zred)
c      b222A = Bispec(2, 2, 2, p, kp, k, zred)
c      b211tA = Bispec(2, 1, 1, kp, p, k, zred)
c      b221tA = Bispec(2, 2, 1, kp, p, k, zred)
c      b212tA = Bispec(2, 1, 2, kp, p, k, zred)
c      b222tA = Bispec(2, 2, 2, kp, p, k, zred)
c
      implicit none
c
      real*8  zred, k, p, kp
      real*8  b211A, b221A, b212A, b222A
      real*8  b211tA, b221tA, b212tA, b222tA
      real*8  sigmav2k, sigmav2p, sigmav2kp, exp_fact
      real*8  dd2, growth, pklin_k, pklin_p, pklin_kp
      real*8  G1_1_k, G1_2_k, G1_1_p, G1_2_p, G1_1_kp, G1_2_kp
      real*8  G2_1_kp_k, G2_2_kp_k, G2_1_p_k, G2_2_p_k, G2_1_p_kp
      real*8  G2_2_p_kp
c     -----------------------------------------   c
c
c     !  Note that this sigma_v2 already includes redshift dependence
      call find_running_sigmav2(k, sigmav2k)
      call find_running_sigmav2(p, sigmav2p)
      call find_running_sigmav2(kp, sigmav2kp)
c
      dd2 = growth(zred)**2
c
      exp_fact = (k*k*sigmav2k + p*p*sigmav2p + kp*kp*sigmav2kp) / 2.d0
c
      if(exp_fact.ge.100) then 
         b211A = 0.d0
         b221A = 0.d0
         b212A = 0.d0
         b222A = 0.d0
         b211tA = 0.d0
         b221tA = 0.d0
         b212tA = 0.d0
         b222tA = 0.d0
         goto 1
      endif
c
      call find_pk(k, pklin_k)
      call find_pk(p, pklin_p)
      call find_pk(kp, pklin_kp)
c
c     bispectrum at tree-level order
c
      call G1reg_1loop(k, p, kp, dd2, sigmav2k, sigmav2p, sigmav2kp,
     &     G1_1_k, G1_2_k, G1_1_p, G1_2_p, G1_1_kp, G1_2_kp)
      call G2reg_1loop(k, p, kp, dd2, sigmav2k, sigmav2p, sigmav2kp,
     &     G2_1_kp_k, G2_2_kp_k, G2_1_p_k, G2_2_p_k, G2_1_p_kp, 
     &     G2_2_p_kp)
c
      b211A = 2.d0 * ( G2_2_kp_k * G1_1_kp * G1_1_k * pklin_kp * pklin_k
     &     + G2_1_p_k * G1_2_p * G1_1_k * pklin_p * pklin_k
     &     + G2_1_p_kp * G1_2_p * G1_1_kp * pklin_p * pklin_kp )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b221A = 2.d0 * ( G2_2_kp_k * G1_2_kp * G1_1_k * pklin_kp * pklin_k
     &     + G2_2_p_k * G1_2_p * G1_1_k * pklin_p * pklin_k
     &     + G2_1_p_kp * G1_2_p * G1_2_kp * pklin_p * pklin_kp )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b212A = 2.d0 * ( G2_2_kp_k * G1_1_kp * G1_2_k * pklin_kp * pklin_k
     &     + G2_1_p_k * G1_2_p * G1_2_k * pklin_p * pklin_k
     &     + G2_2_p_kp * G1_2_p * G1_1_kp * pklin_p * pklin_kp )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b222A = 2.d0 * ( G2_2_kp_k * G1_2_kp * G1_2_k * pklin_kp * pklin_k
     &     + G2_2_p_k * G1_2_p * G1_2_k * pklin_p * pklin_k
     &     + G2_2_p_kp * G1_2_p * G1_2_kp * pklin_p * pklin_kp )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b211tA = 2.d0 * ( G2_2_p_k * G1_1_p * G1_1_k * pklin_p * pklin_k
     &     + G2_1_kp_k * G1_2_kp * G1_1_k * pklin_kp * pklin_k
     &     + G2_1_p_kp * G1_2_kp * G1_1_p * pklin_kp * pklin_p )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b221tA = 2.d0 * ( G2_2_p_k * G1_2_p * G1_1_k * pklin_p * pklin_k
     &     + G2_2_kp_k * G1_2_kp * G1_1_k * pklin_kp * pklin_k
     &     + G2_1_p_kp * G1_2_kp * G1_2_p * pklin_kp * pklin_p )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b212tA = 2.d0 * ( G2_2_p_k * G1_1_p * G1_2_k * pklin_p * pklin_k
     &     + G2_1_kp_k * G1_2_kp * G1_2_k * pklin_kp * pklin_k
     &     + G2_2_p_kp * G1_2_kp * G1_1_p * pklin_kp * pklin_p )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
      b222tA = 2.d0 * ( G2_2_p_k * G1_2_p * G1_2_k * pklin_p * pklin_k
     &     + G2_2_kp_k * G1_2_kp * G1_2_k * pklin_kp * pklin_k
     &     + G2_2_p_kp * G1_2_kp * G1_2_p * pklin_kp * pklin_p )
     &     * dexp(-exp_fact) * dd2**2
cc     &     * dd2**2
c
 1    continue
c
      end
c
c ******************************************************* c
c
      subroutine G1reg_1loop(k, p, kp, dd2, sigmav2k, sigmav2p, 
     &     sigmav2kp, G1_1_k, G1_2_k, G1_1_p, G1_2_p, G1_1_kp, G1_2_kp)
c
c ******************************************************* c
c
c  !  Note that sigma_v2 in the arguments includes redshift dependence
c
      implicit none
c
      real*8  k, p, kp, dd2, sigmav2k, sigmav2p, sigmav2kp
      real*8  G1_1_k, G1_2_k, G1_1_p, G1_2_p, G1_1_kp, G1_2_kp
c     -------------------------- 
c
      if(k.ge.1.d-2) then
         call one_loop_Gamma1(1, k, G1_1_k)
         call one_loop_Gamma1(2, k, G1_2_k)
         G1_1_k = 1.d0 + k*k*sigmav2k/2.d0 + dd2*G1_1_k
         G1_2_k = 1.d0 + k*k*sigmav2k/2.d0 + dd2*G1_2_k
      else
         G1_1_k = 1.d0
         G1_2_k = 1.d0
      endif
c
      if(p.ge.1.d-2 .and. k.ge.1.d-2) then
         call one_loop_Gamma1(1, p, G1_1_p)
         call one_loop_Gamma1(2, p, G1_2_p)
         G1_1_p = 1.d0 + p*p*sigmav2p/2.d0 + dd2*G1_1_p
         G1_2_p = 1.d0 + p*p*sigmav2p/2.d0 + dd2*G1_2_p
      else
         G1_1_p = 1.d0
         G1_2_p = 1.d0
      endif
c
      if(kp.ge.3.d-2 .and. k.ge.1.d-2) then
         call one_loop_Gamma1(1, kp, G1_1_kp)
         call one_loop_Gamma1(2, kp, G1_2_kp)
         G1_1_kp = 1.d0 + kp*kp*sigmav2kp/2.d0 + dd2*G1_1_kp
         G1_2_kp = 1.d0 + kp*kp*sigmav2kp/2.d0 + dd2*G1_2_kp
      else
         G1_1_kp = 1.d0
         G1_2_kp = 1.d0
      endif
c
      end
c
c ******************************************************* c
c
      subroutine G2reg_1loop(k, p, kp, dd2, sigmav2k, sigmav2p, 
     &     sigmav2kp, G2_1_kp_k, G2_2_kp_k, G2_1_p_k, G2_2_p_k, 
     &     G2_1_p_kp, G2_2_p_kp)
c
c ******************************************************* c
c
c  !  Note that sigma_v2 in the arguments includes redshift dependence
c
      implicit none
c
      real*8  k, p, kp, dd2, sigmav2k, sigmav2p, sigmav2kp
      real*8  G2_1_kp_k, G2_2_kp_k, G2_1_p_k, G2_2_p_k
      real*8  G2_1_p_kp, G2_2_p_kp, G2_tree
c     -------------------------- 
c
      if(max(kp, k, p).ge.1d-2 .and. k.ge.1.d-2 ) then 
         call one_loop_Gamma2d(kp, k, p, G2_1_kp_k)
         call one_loop_Gamma2v(kp, k, p, G2_2_kp_k)
         call one_loop_Gamma2d(p, k, kp, G2_1_p_k)
         call one_loop_Gamma2v(p, k, kp, G2_2_p_k)
         call one_loop_Gamma2d(p, kp, k, G2_1_p_kp)
         call one_loop_Gamma2v(p, kp, k, G2_2_p_kp)
c
         G2_1_kp_k = G2_tree(1, kp, k, p) * (1.d0 + p*p*sigmav2p/2.d0)
     &        + dd2*G2_1_kp_k
         G2_2_kp_k = G2_tree(2, kp, k, p) * (1.d0 + p*p*sigmav2p/2.d0)
     &        + dd2*G2_2_kp_k
         G2_1_p_k = G2_tree(1, p, k, kp) *(1.d0 + kp*kp*sigmav2kp/2.d0)
     &        + dd2*G2_1_p_k
         G2_2_p_k = G2_tree(2, p, k, kp) *(1.d0 + kp*kp*sigmav2kp/2.d0)
     &        + dd2*G2_2_p_k
         G2_1_p_kp = G2_tree(1, p, kp, k) *(1.d0 + k*k*sigmav2k/2.d0)
     &        + dd2*G2_1_p_kp
         G2_2_p_kp = G2_tree(2, p, kp, k) *(1.d0 + k*k*sigmav2k/2.d0)
     &        + dd2*G2_2_p_kp
c
      else
         G2_1_kp_k = G2_tree(1, kp, k, p)
         G2_2_kp_k = G2_tree(2, kp, k, p)
         G2_1_p_k = G2_tree(1, p, k, kp)
         G2_2_p_k = G2_tree(2, p, k, kp)
         G2_1_p_kp = G2_tree(1, p, kp, k)
         G2_2_p_kp = G2_tree(2, p, kp, k)
      endif
c     
      end
c
c ******************************************************* c
c
      function G2_tree(a, k1, k2, k3)
c
c ******************************************************* c
c
      implicit none
c
      integer a
      real*8  G2_tree, k1, k2, k3, k1dk2
c     --------------------------------------------------------
c
      k1dk2 = ( k3**2 - k1**2 - k2**2 ) / 2.d0
c
      if(a.eq.1) 
     &     G2_tree = 5.d0/7.d0 + 0.5d0*k1dk2*
     &     (1.d0/k1**2 + 1.d0/k2**2) + 2.d0/7.d0*(k1dk2/(k1*k2))**2
      if(a.eq.2) 
     &     G2_tree = 3.d0/7.d0 + 0.5d0*k1dk2*
     &     (1.d0/k1**2 + 1.d0/k2**2) + 4.d0/7.d0*(k1dk2/(k1*k2))**2
c
      end
c
c
c ******************************************************* c
c
      function LFunc(k,q)
c
c ******************************************************* c
c
      implicit none
      real*8 LFunc, k, q
c     -------------------------------------------------
c
      LFunc = dlog((k + q)**2/(k - q)**2)
c
      end
c
c ******************************************************* c
c
      function WFunc(k1, k2, k3, q)
c
c ******************************************************* c
c
      implicit none
      real*8 aa, bb
      real*8 WFunc, k1, k2, k3, q
c     -------------------------------------------------
c
      aa = -4.d0*k3**2*q**2 - 2.d0*(k1**2 - q**2)*(k2**2 - q**2) 
      bb = 4.d0*k3*q*dsqrt(k1**2*k2**2 + (-k1**2 - k2**2 + k3**2)*q**2 + 
     -     q**4)
c
      WFunc = dlog( (aa-bb)/(aa+bb) )
c
      end
c
c ******************************************************* c
c
      function betafunc(i,z)
c
c ******************************************************* c
c
c     Incomplete beta function of the type, B(z, i, 0)
c     (i=2, 4, 6)
c
c     this is indeed expressed in terms of elementary functions
c
      implicit none
      integer i
      real*8  z, a, betafunc
c     ----------------------------------
c
      betafunc = 0.d0
c
      if(z.lt.0.1d0) then
         a = dble(i)
         betafunc = z**i*(1.d0/a + z/(1.d0 + a) + z**2/(2.d0 + a) + 
     &        z**3/(3.d0 + a) + z**4/(4.d0 + a) + z**5/(5.d0 + a) + 
     &        z**6/(6.d0 + a))
      else
         if(i.eq.2) betafunc =
     &        (z**2*(-2.d0/z - (2.d0*dlog(1.d0 - z))/z**2))/2.d0
         if(i.eq.4) betafunc = 
     &        (z**4*((-2.d0*(6.d0 + 3.d0*z + 2.d0*z**2))/(3.d0*z**3) - 
     &        (4.d0*dlog(1.d0 - z))/z**4))/4.d0
         if(i.eq.6) betafunc =             
     &           (z**6*((-60.d0 - 30.d0*z - 20.d0*z**2 - 15.d0*z**3 - 
     &           12.d0*z**4)/(10.d0*z**5) - 
     &           (6.d0*dlog(1 - z))/z**6))/6.d0
      endif
c
      end
c
c ******************************************************* c
c
      function small_beta(k, q)
c
c ******************************************************* c
c
c     Regular function involving incomplete beta function: 
c
c     (k**2 - q**2) * Beta( 4*k*q/(k+q)**2, 4, 0) / (k + q)**2
c
      implicit none
      real*8 small_beta, k, q
      real*8 y, betafunc
c     ----------------------------------
c
      y = (k - q) / (k + q)
c
      if(dabs(y).lt.1.d-6) then
         small_beta = 0.d0
      else
         small_beta = y * betafunc(4, 1.d0 - y**2) 
      endif
c
cc      write(6,'(A,1p2e18.10)') 'small_beta=',small_beta, y
c
      end
c
c ******************************************************* c
c
      function big_beta(k1, k2, k3, q)
c
c ******************************************************* c
c
c     Regular function involving incomplete beta function:  
c
c     x * Beta(1-X**2, 2, 0) / sqrt( x + a**2 )
c
c     where a, x and X are defined as
c
c     a = k3 * q
c     x = (k1**2 - q**2) * (k2**2 - q**2)
c     X = {1 - sqrt(x + a**2)} / {1 + sqrt(x + a**2)} 
c

      implicit none
      real*8 big_beta, k1, k2, k3, q
      real*8 x, a, y, betafunc
c     ----------------------------------
c
      x = (k1**2-q**2) * (k2**2-q**2)
      a = k3 * q
      y = dsqrt(x+a*a)/a
c     
      if(dabs(y-1).lt.1.d-5) then
         big_beta = 0.d0
      elseif(dabs(y-1).lt.1.d-2) then
         big_beta = 2.d0*a*(-1.d0 + y)*(-1.d0 + dlog(4.d0) - 
     &        2.d0*dlog(dabs(-1.d0 + y))) + a*(-1.d0 + y)**2*
     &        (3.d0 - dlog(4.d0) + 2.d0*dlog(dabs(-1.d0 + y)))
cc         write(6,*) 'asymptotic expansion'
      else 
         big_beta = x * betafunc(2, 1-(1.d0-y)**2/(1.d0+y)**2) / (a*y)
cc         write(6,*) 'exact expression'
      endif
c
ccc      write(6,'(A,1p1e18.10)') 'big_beta=',big_beta
c
      end
