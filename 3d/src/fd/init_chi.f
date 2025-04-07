      subroutine init_chi(x,y,z,Nx,Ny,Nz,amp,c,chiX0,chiY0,chiZ0,d,delta
     &,r0,wwX,wwY,wwZ,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 amp
      real*8 c
      real*8 chiX0
      real*8 chiY0
      real*8 chiZ0
      real*8 d
      real*8 delta
      real*8 r0
      real*8 wwX
      real*8 wwY
      real*8 wwZ
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = amp * exp(-0.1D1 * (sqrt((d * exp(c * x(i)) - 0.1D1 * d * exp
     #(-0.1D1 * c * x(i)) - 0.1D1 * chiX0) ** 2 / wwX ** 2 + (d * exp(c 
     #* y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j)) - 0.1D1 * chiY0) ** 2
     # / wwY ** 2 + (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(
     #k)) - 0.1D1 * chiZ0) ** 2 / wwZ ** 2) - 0.1D1 * r0) ** 2 / delta *
     #* 2)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
