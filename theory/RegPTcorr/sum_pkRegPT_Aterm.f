c ******************************************************* c
c
      program sum_pkRegPT_Aterm
c
c ******************************************************* c      
c
      implicit none
      integer  ik, ikmax, ik_max, ik_max1, ik_max2
      parameter(ikmax=3000)
      real*8  ak, pkI_A11, pkII_A11 ! mu^2 * f
      real*8  pkI_A12, pkII_A12     ! mu^2 * f**2
      real*8  pkI_A22, pkII_A22     ! mu^4 * f**2
      real*8  pkI_A23, pkII_A23     ! mu^4 * f**3
      real*8  pkI_A33, pkII_A33     ! mu^6 * f**3
      character ifile1*128, ifile2*128, outfile*128
c
      call getarg(1, ifile1)
      call getarg(2, ifile2)
      call getarg(3, outfile)
      open(10, file=ifile1, status='unknown')
      open(11, file=ifile2, status='unknown')
      !open(12, file='pkRegPT_Aterm.dat', status='unknown')
      open(12, file=outfile, status='unknown')
      do ik=1, ikmax
         read(10,*,END=10) ak, pkI_A11, pkI_A12, pkI_A22, pkI_A23,
     &        pkI_A33
         read(11,*,END=10) ak, pkII_A11, pkII_A12, pkII_A22, pkII_A23,
     &        pkII_A33
         write(12,'(1p6e18.10)') ak, pkI_A11+pkII_A11, pkI_A12+pkII_A12,
     &        pkI_A22+pkII_A22, pkI_A23+pkII_A23, pkI_A33+pkII_A33
      enddo
 10   close(10)
      close(11)
      close(12)
      end
