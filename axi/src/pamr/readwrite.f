cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     w2f: write multigrid data to binary file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine w2f(V,NR,NX)

      implicit none

      integer NR,NX
      real*8 V(NR,NX)

      integer i, j
      integer unum, nbytes

      nbytes = 8 * NR * NX
c      write(*,*) nbytes

      unum=20

      open(unit=unum, file='V.bin', status='replace',
     &     form='unformatted', access='direct', recl=8)

      do j = 1, NX
      do i = 1, NR
        write(unum, rec=((j-1)*NR+i)) V(i,j)
      end do
      end do

      close(unum)

      return
      end subroutine w2f

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     rff: read multigrid data from binary file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rff(V,NR,NX)

      implicit none

      integer NR,NX
      real*8 V(NR,NX)

      integer i, j
      integer unum, nbytes

      nbytes = 8 * NR * NX
c      write(*,*) nbytes

      unum=20

      open(unit=unum, file='V.bin', status='old',
     &     form='unformatted', access='direct', recl=8)

      do j = 1, NX
      do i = 1, NR
        read(unum, rec=((j-1)*NR+i)) V(i,j)
      end do
      end do

      close(unum)

      return
      end subroutine rff
