      subroutine init_dydy2(x,y,z,Nx,Ny,Nz,c,d,myzero,res)
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
      qb = -0.1D1 * (exp(0.2D1 * c * y(j)) - 0.1D1) * exp(c * y(j)) / d 
     #/ c / (exp(0.4D1 * c * y(j)) + 0.2D1 * exp(0.2D1 * c * y(j)) + 0.1
     #D1) * (d ** 2 * (exp(0.2D1 * c * y(j)) + 0.1D1) ** 2 * exp(-0.2D1 
     #* c * y(j))) ** (-0.1D1 / 0.2D1) + myzero * x(i) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
