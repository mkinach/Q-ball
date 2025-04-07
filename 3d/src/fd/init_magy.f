      subroutine init_magy(n_Ax,n_Az,n_dxdX,n_dzdZ,Nx,Ny,Nz,hx,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hx
      real*8 hz
      real*8 n_Ax(Nx,Ny,Nz)
      real*8 n_Az(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=1, Ny, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (n_Ax(i, j, k - 2) - 0.8D1 * n_Ax(i, 
     #j, k - 1) + 0.8D1 * n_Ax(i, j, k + 1) - 0.1D1 * n_Ax(i, j, k + 2))
     # / hz * n_dzdZ(i, j, k) - 0.8333333333333333D-1 * (n_Az(i - 2, j, 
     #k) - 0.8D1 * n_Az(i - 1, j, k) + 0.8D1 * n_Az(i + 1, j, k) - 0.1D1
     # * n_Az(i + 2, j, k)) / hx * n_dxdX(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
