      subroutine init_elecy(n_At,n_By,n_dydY,Nx,Ny,Nz,hy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hy
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_By(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=3, Ny-2, 1
      do k=1, Nz, 1
      qb = 0.8333333333333333D-1 * (n_At(i, j - 2, k) - 0.8D1 * n_At(i, 
     #j - 1, k) + 0.8D1 * n_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, k))
     # / hy * n_dydY(i, j, k) - 0.1D1 * n_By(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
