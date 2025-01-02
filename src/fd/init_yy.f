      subroutine init_yy(x,y,z,Nx,Ny,Nz,c,d,myzero,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 c
      real*8 d
      real*8 myzero
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = myzero * x(i) * z(k) + d * exp(c * y(j)) - 0.1D1 * d * exp(-0
     #.1D1 * c * y(j))
      res(i,j,k)=qb
      end do
      end do
      end do
      END
