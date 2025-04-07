      subroutine init_magx(n_Ay,n_Az,n_dydY,n_dzdZ,Nx,Ny,Nz,hy,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hy
      real*8 hz
      real*8 n_Ay(Nx,Ny,Nz)
      real*8 n_Az(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (n_Az(i, j - 2, k) - 0.8D1 * n_Az(i, 
     #j - 1, k) + 0.8D1 * n_Az(i, j + 1, k) - 0.1D1 * n_Az(i, j + 2, k))
     # / hy * n_dydY(i, j, k) - 0.8333333333333333D-1 * (n_Ay(i, j, k - 
     #2) - 0.8D1 * n_Ay(i, j, k - 1) + 0.8D1 * n_Ay(i, j, k + 1) - 0.1D1
     # * n_Ay(i, j, k + 2)) / hz * n_dzdZ(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
