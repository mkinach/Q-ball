      subroutine init_elecz(n_At,n_Bz,n_dzdZ,Nx,Ny,Nz,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hz
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_Bz(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (n_At(i, j, k - 2) - 0.8D1 * n_At(i, 
     #j, k - 1) + 0.8D1 * n_At(i, j, k + 1) - 0.1D1 * n_At(i, j, k + 2))
     # / hz * n_dzdZ(i, j, k) - 0.1D1 * n_Bz(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
