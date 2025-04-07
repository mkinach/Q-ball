      subroutine ire_phi2_poly(n_At,n_Ax,n_Ay,n_Az,n_chi,n_dxdX,n_dxdX2,
     &n_dydY,n_dydY2,n_dzdZ,n_dzdZ2,n_phi1,n_phi2,n_pi1,nm1_phi2,np1_phi
     &2,Nx,Ny,Nz,cc,e,g,h,ht,hx,hy,hz,m,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 cc
      real*8 e
      real*8 g
      real*8 h
      real*8 ht
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 m
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_Ax(Nx,Ny,Nz)
      real*8 n_Ay(Nx,Ny,Nz)
      real*8 n_Az(Nx,Ny,Nz)
      real*8 n_chi(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dxdX2(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dydY2(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 n_dzdZ2(Nx,Ny,Nz)
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 n_pi1(Nx,Ny,Nz)
      real*8 nm1_phi2(Nx,Ny,Nz)
      real*8 np1_phi2(Nx,Ny,Nz)
      real*8 res
      real*8 qb
      res = 0.0D0
      do i=2, Nx-1, 1
      do j=2, Ny-1, 1
      do k=2, Nz-1, 1
      qb = 0.1D1 / ht ** 2 * (nm1_phi2(i, j, k) - 0.2D1 * n_phi2(i, j, k
     #) + np1_phi2(i, j, k)) - 0.1D1 * (n_phi2(i, j, k) * n_At(i, j, k) 
     #** 2 - 0.1D1 * n_phi2(i, j, k) * n_Ax(i, j, k) ** 2 - 0.1D1 * n_ph
     #i2(i, j, k) * n_Az(i, j, k) ** 2 - 0.1D1 * n_phi2(i, j, k) * n_Ay(
     #i, j, k) ** 2) * e ** 2 - 0.1D1 * (0.2D1 * n_pi1(i, j, k) * n_At(i
     #, j, k) - 0.1D1 * (-0.1D1 * n_phi1(i, j, k - 1) + n_phi1(i, j, k +
     # 1)) / hz * n_dzdZ(i, j, k) * n_Az(i, j, k) - 0.1D1 * (-0.1D1 * n_
     #phi1(i - 1, j, k) + n_phi1(i + 1, j, k)) / hx * n_dxdX(i, j, k) * 
     #n_Ax(i, j, k) - 0.1D1 * (-0.1D1 * n_phi1(i, j - 1, k) + n_phi1(i, 
     #j + 1, k)) / hy * n_dydY(i, j, k) * n_Ay(i, j, k)) * e - 0.1D1 / h
     #y ** 2 * (n_phi2(i, j - 1, k) - 0.2D1 * n_phi2(i, j, k) + n_phi2(i
     #, j + 1, k)) * n_dydY(i, j, k) ** 2 - 0.5000000000000000D0 * (-0.1
     #D1 * n_phi2(i, j - 1, k) + n_phi2(i, j + 1, k)) / hy * n_dydY2(i, 
     #j, k) - 0.1D1 / hz ** 2 * (n_phi2(i, j, k - 1) - 0.2D1 * n_phi2(i,
     # j, k) + n_phi2(i, j, k + 1)) * n_dzdZ(i, j, k) ** 2 - 0.500000000
     #0000000D0 * (-0.1D1 * n_phi2(i, j, k - 1) + n_phi2(i, j, k + 1)) /
     # hz * n_dzdZ2(i, j, k) - 0.1D1 / hx ** 2 * (n_phi2(i - 1, j, k) - 
     #0.2D1 * n_phi2(i, j, k) + n_phi2(i + 1, j, k)) * n_dxdX(i, j, k) *
     #* 2 - 0.5000000000000000D0 * (-0.1D1 * n_phi2(i - 1, j, k) + n_phi
     #2(i + 1, j, k)) / hx * n_dxdX2(i, j, k) + m ** 2 * n_phi2(i, j, k)
     # - 0.1D1 * g * (n_phi1(i, j, k) ** 2 + n_phi2(i, j, k) ** 2) * n_p
     #hi2(i, j, k) + h * (n_phi1(i, j, k) ** 2 + n_phi2(i, j, k) ** 2) *
     #* 2 * n_phi2(i, j, k) + cc * n_chi(i, j, k) ** 2 * n_phi2(i, j, k)
      res = res + qb**2
      end do
      end do
      end do
      res = sqrt(res/(1*Nx*Ny*Nz))
      END
