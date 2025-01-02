      subroutine ire_ax(n_Ax,n_dxdX,n_dxdX2,n_dydY,n_dydY2,n_dzdZ,n_dzdZ
     &2,n_phi1,n_phi2,nm1_Ax,np1_Ax,Nx,Ny,Nz,e,ht,hx,hy,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 e
      real*8 ht
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 n_Ax(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dxdX2(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dydY2(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 n_dzdZ2(Nx,Ny,Nz)
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 nm1_Ax(Nx,Ny,Nz)
      real*8 np1_Ax(Nx,Ny,Nz)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      do j=2, Ny-1, 1
      do k=2, Nz-1, 1
      qb = 0.1D1 / ht ** 2 * (nm1_Ax(i, j, k) - 0.2D1 * n_Ax(i, j, k) + 
     #np1_Ax(i, j, k)) + (0.2D1 * n_Ax(i, j, k) * n_phi1(i, j, k) ** 2 +
     # 0.2D1 * n_Ax(i, j, k) * n_phi2(i, j, k) ** 2) * e ** 2 + (-0.1D1 
     #* n_phi1(i, j, k) * (-0.1D1 * n_phi2(i - 1, j, k) + n_phi2(i + 1, 
     #j, k)) / hx * n_dxdX(i, j, k) + n_phi2(i, j, k) * (-0.1D1 * n_phi1
     #(i - 1, j, k) + n_phi1(i + 1, j, k)) / hx * n_dxdX(i, j, k)) * e -
     # 0.1D1 / hx ** 2 * (n_Ax(i - 1, j, k) - 0.2D1 * n_Ax(i, j, k) + n_
     #Ax(i + 1, j, k)) * n_dxdX(i, j, k) ** 2 - 0.5000000000000000D0 * (
     #-0.1D1 * n_Ax(i - 1, j, k) + n_Ax(i + 1, j, k)) / hx * n_dxdX2(i, 
     #j, k) - 0.1D1 / hy ** 2 * (n_Ax(i, j - 1, k) - 0.2D1 * n_Ax(i, j, 
     #k) + n_Ax(i, j + 1, k)) * n_dydY(i, j, k) ** 2 - 0.500000000000000
     #0D0 * (-0.1D1 * n_Ax(i, j - 1, k) + n_Ax(i, j + 1, k)) / hy * n_dy
     #dY2(i, j, k) - 0.1D1 / hz ** 2 * (n_Ax(i, j, k - 1) - 0.2D1 * n_Ax
     #(i, j, k) + n_Ax(i, j, k + 1)) * n_dzdZ(i, j, k) ** 2 - 0.50000000
     #00000000D0 * (-0.1D1 * n_Ax(i, j, k - 1) + n_Ax(i, j, k + 1)) / hz
     # * n_dzdZ2(i, j, k)
      res = res + qb**2
      end do
      end do
      end do
      res = sqrt(res/(1*Nx*Ny*Nz))
      END
