      subroutine init_elecx(n_At,n_Bx,n_dxdX,Nx,Ny,Nz,hx,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hx
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_Bx(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = 0.8333333333333333D-1 * (n_At(i - 2, j, k) - 0.8D1 * n_At(i -
     # 1, j, k) + 0.8D1 * n_At(i + 1, j, k) - 0.1D1 * n_At(i + 2, j, k))
     # / hx * n_dxdX(i, j, k) - 0.1D1 * n_Bx(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
