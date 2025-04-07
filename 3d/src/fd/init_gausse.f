      subroutine init_gausse(n_At,n_Bx,n_By,n_Bz,n_dxdX,n_dxdX2,n_dydY,n
     &_dydY2,n_dzdZ,n_dzdZ2,n_qden,Nx,Ny,Nz,e,hx,hy,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 e
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_Bx(Nx,Ny,Nz)
      real*8 n_By(Nx,Ny,Nz)
      real*8 n_Bz(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dxdX2(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dydY2(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 n_dzdZ2(Nx,Ny,Nz)
      real*8 n_qden(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (-0.1D1 * n_At(i - 2, j, k) + 0.16D2 
     #* n_At(i - 1, j, k) - 0.30D2 * n_At(i, j, k) + 0.16D2 * n_At(i + 1
     #, j, k) - 0.1D1 * n_At(i + 2, j, k)) / hx ** 2 * n_dxdX(i, j, k) *
     #* 2 + 0.8333333333333333D-1 * (n_At(i - 2, j, k) - 0.8D1 * n_At(i 
     #- 1, j, k) + 0.8D1 * n_At(i + 1, j, k) - 0.1D1 * n_At(i + 2, j, k)
     #) / hx * n_dxdX2(i, j, k) - 0.8333333333333333D-1 * (n_Bx(i - 2, j
     #, k) - 0.8D1 * n_Bx(i - 1, j, k) + 0.8D1 * n_Bx(i + 1, j, k) - 0.1
     #D1 * n_Bx(i + 2, j, k)) / hx * n_dxdX(i, j, k) + 0.833333333333333
     #3D-1 * (-0.1D1 * n_At(i, j - 2, k) + 0.16D2 * n_At(i, j - 1, k) - 
     #0.30D2 * n_At(i, j, k) + 0.16D2 * n_At(i, j + 1, k) - 0.1D1 * n_At
     #(i, j + 2, k)) / hy ** 2 * n_dydY(i, j, k) ** 2 + 0.83333333333333
     #33D-1 * (n_At(i, j - 2, k) - 0.8D1 * n_At(i, j - 1, k) + 0.8D1 * n
     #_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, k)) / hy * n_dydY2(i, j,
     # k) - 0.8333333333333333D-1 * (n_By(i, j - 2, k) - 0.8D1 * n_By(i,
     # j - 1, k) + 0.8D1 * n_By(i, j + 1, k) - 0.1D1 * n_By(i, j + 2, k)
     #) / hy * n_dydY(i, j, k) + 0.8333333333333333D-1 * (-0.1D1 * n_At(
     #i, j, k - 2) + 0.16D2 * n_At(i, j, k - 1) - 0.30D2 * n_At(i, j, k)
     # + 0.16D2 * n_At(i, j, k + 1) - 0.1D1 * n_At(i, j, k + 2)) / hz **
     # 2 * n_dzdZ(i, j, k) ** 2 + 0.8333333333333333D-1 * (n_At(i, j, k 
     #- 2) - 0.8D1 * n_At(i, j, k - 1) + 0.8D1 * n_At(i, j, k + 1) - 0.1
     #D1 * n_At(i, j, k + 2)) / hz * n_dzdZ2(i, j, k) - 0.83333333333333
     #33D-1 * (n_Bz(i, j, k - 2) - 0.8D1 * n_Bz(i, j, k - 1) + 0.8D1 * n
     #_Bz(i, j, k + 1) - 0.1D1 * n_Bz(i, j, k + 2)) / hz * n_dzdZ(i, j, 
     #k) - 0.1D1 * e * n_qden(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
