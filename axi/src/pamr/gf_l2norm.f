cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gfl2norm: manually computes the L2-norm of an arbitrary function
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function gfl2norm(gf_np1,NR,NX,gw)

      implicit none

      include 'globals.inc'

      real*8  gf_np1(1:NR,1:NX)
      integer NR,NX
      integer gw(4)
      integer i, j
      integer gwc(4)
      real*8  norm
      integer numpts

      norm   = 0.0d0
      numpts = 0

      ! make a local copy of ghost widths
      gwc(1) = gw(1)
      gwc(2) = gw(2)
      gwc(3) = gw(3)
      gwc(4) = gw(4)

      ! does not consider ghost points
      do i=1+gwc(1), NR-gwc(2)
      do j=1+gwc(3), NX-gwc(4)
        norm   = norm + gf_np1(i,j)**2d0
        numpts = numpts+1
      end do
      end do

      gfl2norm = norm

      return
      end function gfl2norm
