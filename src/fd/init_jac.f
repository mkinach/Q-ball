      subroutine init_jac(x,y,z,Nx,Ny,Nz,c,d,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 c
      real*8 d
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = c ** 3 * d ** 3 * exp(-0.1D1 * c * (x(i) + y(j) + z(k))) * (e
     #xp(0.2D1 * c * x(i)) + 0.1D1) * (exp(0.2D1 * c * y(j)) + 0.1D1) * 
     #(exp(0.2D1 * c * z(k)) + 0.1D1)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
