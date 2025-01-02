      subroutine init_magz(n_Ax,n_Ay,n_dxdX,n_dydY,Nx,Ny,Nz,hx,hy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hx
      real*8 hy
      real*8 n_Ax(Nx,Ny,Nz)
      real*8 n_Ay(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=1, Nz, 1
      qb = 0.8333333333333333D-1 * (n_Ay(i - 2, j, k) - 0.8D1 * n_Ay(i -
     # 1, j, k) + 0.8D1 * n_Ay(i + 1, j, k) - 0.1D1 * n_Ay(i + 2, j, k))
     # / hx * n_dxdX(i, j, k) - 0.8333333333333333D-1 * (n_Ax(i, j - 2, 
     #k) - 0.8D1 * n_Ax(i, j - 1, k) + 0.8D1 * n_Ax(i, j + 1, k) - 0.1D1
     # * n_Ax(i, j + 2, k)) / hy * n_dydY(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
