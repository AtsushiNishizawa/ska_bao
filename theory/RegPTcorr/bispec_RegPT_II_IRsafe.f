c ******************************************************* c
c
      subroutine Bispec_II_IRsafe(kk1, kk2, kk3, qq, a, b, c, 
     &     bk222_abc, bk321_abc)
c                                                         c
c            Time-stamp: <2016-11-30 11:39:15 ataruya>    c
c ******************************************************* c
c
c     integrand of bispectrum RegPT at one-loop order --- II ---
c
c     Note !!! 
c     Exponential damping is not yet multiplied in this subroutine.
c
      implicit none
c
      integer a, b, c, i
      real*8  kk1(3), kk2(3), kk3(3), qq(3)
      real*8  pp1(3), pp2(3), pp3(3), rr1(3), rr2(3), rr3(3)
      real*8  q, k1, k2, k3, p1, p2, p3, r1, r2, r3
      real*8  pk1, pk2, pk3
      real*8  pkq, pkp1, pkp2, pkp3, pkr1, pkr2, pkr3
      real*8  bk222_abc, bk321_abc
      real*8  F2_sym, F3_sym, pi
      pi = 4.d0 * datan(1.d0)
c     -----------------------------------------   c
c
      bk222_abc = 0.d0
      bk321_abc = 0.d0
c
      do i=1, 3
         pp1(i) = kk1(i) - qq(i)
         pp2(i) = kk2(i) - qq(i)
         pp3(i) = kk3(i) - qq(i)
         rr1(i) = kk1(i) + qq(i)
         rr2(i) = kk2(i) + qq(i)
         rr3(i) = kk3(i) + qq(i)
      enddo
      p1 = dsqrt(pp1(1)**2 + pp1(2)**2 + pp1(3)**2)
      p2 = dsqrt(pp2(1)**2 + pp2(2)**2 + pp2(3)**2)
      p3 = dsqrt(pp3(1)**2 + pp3(2)**2 + pp3(3)**2)
      r1 = dsqrt(rr1(1)**2 + rr1(2)**2 + rr1(3)**2)
      r2 = dsqrt(rr2(1)**2 + rr2(2)**2 + rr2(3)**2)
      r3 = dsqrt(rr3(1)**2 + rr3(2)**2 + rr3(3)**2)
c
      q  = dsqrt( qq(1)**2 + qq(2)**2  + qq(3)**2 )  
      k1 = dsqrt( kk1(1)**2 + kk1(2)**2 + kk1(3)**2 ) 
      k2 = dsqrt( kk2(1)**2 + kk2(2)**2 + kk2(3)**2 ) 
      k3 = dsqrt( kk3(1)**2 + kk3(2)**2 + kk3(3)**2 ) 
c
      call find_pk(q,  pkq )
      call find_pk(k1, pk1 )
      call find_pk(k2, pk2 )
      call find_pk(k3, pk3 )
      call find_pk(p1, pkp1)
      call find_pk(p2, pkp2)
      call find_pk(p3, pkp3)
      call find_pk(r1, pkr1)
      call find_pk(r2, pkr2)
      call find_pk(r3, pkr3)
c
c     /////// IR-safe integrand for bk222  ///////   c
c
      if(p1.gt.q .and. r2.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(a, pp1, qq) * F2_sym(b, rr2, -qq) *
     &        F2_sym(c, -rr2, -pp1) * pkp1 * pkq * pkr2
      endif
      if(r1.gt.q .and. p2.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(a, rr1, -qq) * F2_sym(b, pp2, qq) *
     &        F2_sym(c, -pp2, -rr1) * pkr1 * pkq * pkp2 
      endif
c
      if(p3.gt.q .and. r2.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(c, pp3, qq) * F2_sym(b, rr2, -qq) *
     &        F2_sym(a, -rr2, -pp3) * pkp3 * pkq * pkr2 
      endif
      if(r3.gt.q .and. p2.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(c, rr3, -qq) * F2_sym(b, pp2, qq) *
     &        F2_sym(a, -pp2, -rr3) * pkr3 * pkq * pkp2 
      endif
c
      if(p1.gt.q .and. r3.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(a, pp1, qq) * F2_sym(c, rr3, -qq) *
     &        F2_sym(b, -rr3, -pp1) * pkp1 * pkq * pkr3 
      endif
      if(r1.gt.q .and. p3.gt.q) then
         bk222_abc = bk222_abc + 
     &        F2_sym(a, rr1, -qq) * F2_sym(c, pp3, qq) *
     &        F2_sym(b, -pp3, -rr1) * pkr1 * pkq * pkp3 
      endif
c
      bk222_abc = bk222_abc / 2.d0
c
c     /////// IR-safe integrand for bk321  ///////   c
c
      if(p2.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(a, -kk3, -pp2, -qq) * F2_sym(b, pp2, qq) 
     &        * pkp2 * pkq * pk3 +
     &        F3_sym(c, -kk1, -pp2, -qq) * F2_sym(b, pp2, qq) 
     &        * pkp2 * pkq * pk1 )
      endif
      if(r2.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(a, -kk3, -rr2, qq) * F2_sym(b, rr2, -qq) 
     &        * pkr2 * pkq * pk3 +
     &        F3_sym(c, -kk1, -rr2, qq) * F2_sym(b, rr2, -qq) 
     &        * pkr2 * pkq * pk1 )
      endif
c
      if(p3.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(a, -kk2, -pp3, -qq) * F2_sym(c, pp3, qq) 
     &        * pkp3 * pkq * pk2 +
     &        F3_sym(b, -kk1, -pp3, -qq) * F2_sym(c, pp3, qq) 
     &        * pkp3 * pkq * pk1 ) 
      endif
      if(r3.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(a, -kk2, -rr3, qq) * F2_sym(c, rr3, -qq) 
     &        * pkr3 * pkq * pk2 +
     &        F3_sym(b, -kk1, -rr3, qq) * F2_sym(c, rr3, -qq) 
     &        * pkr3 * pkq * pk1 )
      endif
c
      if(p1.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(b, -kk3, -pp1, -qq) * F2_sym(a, pp1, qq) 
     &        * pkp1 * pkq * pk3 +
     &        F3_sym(c, -kk2, -pp1, -qq) * F2_sym(a, pp1, qq) 
     &        * pkp1 * pkq * pk2 )
      endif
      if(r1.gt.q) then
         bk321_abc = bk321_abc + (
     &        F3_sym(b, -kk3, -rr1, qq) * F2_sym(a, rr1, -qq) 
     &        * pkr1 * pkq * pk3 +
     &        F3_sym(c, -kk2, -rr1, qq) * F2_sym(a, rr1, -qq) 
     &        * pkr1 * pkq * pk2 )
      endif
c
      bk222_abc = 8.d0 * bk222_abc / (2.*pi)**3
      bk321_abc = 6.d0 * bk321_abc / (2.*pi)**3
c
      end
c
c ******************************************************* c
c
      function sigmaab(n, a, b)
c
c ******************************************************* c
c
      implicit none
      integer n, a, b
      real*8  sigmaab
c
      if(a.eq.1 .and. b.eq.1) then
         sigmaab = 2.d0*dble(n) + 1.d0
      elseif(a.eq.1 .and. b.eq.2) then
         sigmaab = 2.d0
      elseif(a.eq.2 .and. b.eq.1) then
         sigmaab = 3.d0
      elseif(a.eq.2 .and. b.eq.2) then
         sigmaab = 2.d0*dble(n)
      else
         sigmaab = 0.d0
      endif
c
      sigmaab = sigmaab / (2.d0*dble(n) + 3.d0) / (dble(n) - 1.d0)
c
      end
c
c ******************************************************* c
c
      function gam_matrix(a, b, c, p, q)
c
c ******************************************************* c
c
      implicit none
      integer a, b, c
      real*8  gam_matrix, p(3), q(3), pp, qq, pq
c
      pp =  p(1)**2 + p(2)**2 + p(3)**2 
      qq =  q(1)**2 + q(2)**2 + q(3)**2 
      pq =  p(1)*q(1) + p(2)*q(2) + p(3)*q(3)
c
      if(a.eq.1 .and. b.eq.1 .and. c.eq.2) then
         gam_matrix = 1.d0 + pq / qq 
      elseif(a.eq.1 .and. b.eq.2 .and. c.eq.1) then
         gam_matrix = 1.d0 + pq / pp 
      elseif(a.eq.2 .and. b.eq.2 .and. c.eq.2) then
         gam_matrix = pq * ( pp+qq+2.d0*pq ) /( pp*qq )
      else
         gam_matrix = 0.d0
      endif
c
      gam_matrix = gam_matrix / 2.d0
c
      end
c
c ******************************************************* c
c
      function F2_sym(a, p, q)
c
c ******************************************************* c
c
      implicit none
      integer a
      real*8 F2_sym, p(3), q(3), pp, qq, mu
c     -------------------------------------------
c
      pp = dsqrt( p(1)**2 + p(2)**2 + p(3)**2 )
      qq = dsqrt( q(1)**2 + q(2)**2 + q(3)**2 )
      mu = ( p(1)*q(1) + p(2)*q(2) + p(3)*q(3) )/( pp*qq )
c
      if(a.eq.1) F2_sym = 
     &     5.d0/7.d0 + mu/2.d0*(pp/qq + qq/pp) + 2.d0/7.d0*mu**2 
      if(a.eq.2) F2_sym = 
     &     3.d0/7.d0 + mu/2.d0*(pp/qq + qq/pp) + 4.d0/7.d0*mu**2 
c
      end
c
c ******************************************************* c
c
      function F3(a, p, q, r)
c
c ******************************************************* c
c
      implicit none
      integer a, i
      real*8 F3, F2_sym, p(3), q(3), r(3), qr(3), sigmaab, gam_matrix
c
      qr(1) = q(1) + r(1)
      qr(2) = q(2) + r(2)
      qr(3) = q(3) + r(3)
c
      if( (dabs(qr(1)).lt.1.d-30 .and. dabs(qr(2)).lt.1.d-30 .and. 
     &     dabs(qr(3)).lt.1.d-30)  .or. 
     &     (dabs(p(1)+qr(1)).lt.1.d-30 .and. dabs(p(2)+qr(2)).lt.1.d-30 
     &     .and. dabs(p(3)+qr(3)).lt.1.d-30) 
     &     ) then
         F3 = 0.d0
      else
         F3 = ( sigmaab(3,a,1) * gam_matrix(1,1,2,p,qr) 
     &        + sigmaab(3,a,2) * gam_matrix(2,2,2,p,qr) ) 
     &        * F2_sym(2,q,r) 
     &        + sigmaab(3,a,1) * gam_matrix(1,2,1,p,qr) 
     &        * F2_sym(1,q,r)
         F3 = 2.d0 * F3
      endif
c
ccc      write(6,*) 'F3=',F3, (p(i),i=1,3),(q(i),i=1,3),(r(i),i=1,3)
      end
c
c ******************************************************* c
c
      function F3_sym(a, p, q, r)
c
c ******************************************************* c
c
c     fully symmetric F3 :  p <--> q <--> r
c
      implicit none
      integer a
      real*8 F3_sym, F3, p(3), q(3), r(3)
c
      F3_sym = ( F3(a, p, q, r) + F3(a, r, p, q) + F3(a, q, r, p) )/3.d0
c
      end
c
