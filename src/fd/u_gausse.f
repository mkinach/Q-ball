      subroutine u_gausse(np1_At,np1_Bx,np1_By,np1_Bz,np1_dxdX,np1_dxdX2
     &,np1_dydY,np1_dydY2,np1_dzdZ,np1_dzdZ2,np1_qden,Nx,Ny,Nz,e,hx,hy,h
     &z,phys_bdy,res)
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
      real*8 np1_At(Nx,Ny,Nz)
      real*8 np1_Bx(Nx,Ny,Nz)
      real*8 np1_By(Nx,Ny,Nz)
      real*8 np1_Bz(Nx,Ny,Nz)
      real*8 np1_dxdX(Nx,Ny,Nz)
      real*8 np1_dxdX2(Nx,Ny,Nz)
      real*8 np1_dydY(Nx,Ny,Nz)
      real*8 np1_dydY2(Nx,Ny,Nz)
      real*8 np1_dzdZ(Nx,Ny,Nz)
      real*8 np1_dzdZ2(Nx,Ny,Nz)
      real*8 np1_qden(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (-0.1D1 * np1_At(i - 2, j, k) + 0.16D
     #2 * np1_At(i - 1, j, k) - 0.30D2 * np1_At(i, j, k) + 0.16D2 * np1_
     #At(i + 1, j, k) - 0.1D1 * np1_At(i + 2, j, k)) / hx ** 2 * np1_dxd
     #X(i, j, k) ** 2 + 0.8333333333333333D-1 * (np1_At(i - 2, j, k) - 0
     #.8D1 * np1_At(i - 1, j, k) + 0.8D1 * np1_At(i + 1, j, k) - 0.1D1 *
     # np1_At(i + 2, j, k)) / hx * np1_dxdX2(i, j, k) - 0.83333333333333
     #33D-1 * (np1_Bx(i - 2, j, k) - 0.8D1 * np1_Bx(i - 1, j, k) + 0.8D1
     # * np1_Bx(i + 1, j, k) - 0.1D1 * np1_Bx(i + 2, j, k)) / hx * np1_d
     #xdX(i, j, k) + 0.8333333333333333D-1 * (-0.1D1 * np1_At(i, j - 2, 
     #k) + 0.16D2 * np1_At(i, j - 1, k) - 0.30D2 * np1_At(i, j, k) + 0.1
     #6D2 * np1_At(i, j + 1, k) - 0.1D1 * np1_At(i, j + 2, k)) / hy ** 2
     # * np1_dydY(i, j, k) ** 2 + 0.8333333333333333D-1 * (np1_At(i, j -
     # 2, k) - 0.8D1 * np1_At(i, j - 1, k) + 0.8D1 * np1_At(i, j + 1, k)
     # - 0.1D1 * np1_At(i, j + 2, k)) / hy * np1_dydY2(i, j, k) - 0.8333
     #333333333333D-1 * (np1_By(i, j - 2, k) - 0.8D1 * np1_By(i, j - 1, 
     #k) + 0.8D1 * np1_By(i, j + 1, k) - 0.1D1 * np1_By(i, j + 2, k)) / 
     #hy * np1_dydY(i, j, k) + 0.8333333333333333D-1 * (-0.1D1 * np1_At(
     #i, j, k - 2) + 0.16D2 * np1_At(i, j, k - 1) - 0.30D2 * np1_At(i, j
     #, k) + 0.16D2 * np1_At(i, j, k + 1) - 0.1D1 * np1_At(i, j, k + 2))
     # / hz ** 2 * np1_dzdZ(i, j, k) ** 2 + 0.8333333333333333D-1 * (np1
     #_At(i, j, k - 2) - 0.8D1 * np1_At(i, j, k - 1) + 0.8D1 * np1_At(i,
     # j, k + 1) - 0.1D1 * np1_At(i, j, k + 2)) / hz * np1_dzdZ2(i, j, k
     #) - 0.8333333333333333D-1 * (np1_Bz(i, j, k - 2) - 0.8D1 * np1_Bz(
     #i, j, k - 1) + 0.8D1 * np1_Bz(i, j, k + 1) - 0.1D1 * np1_Bz(i, j, 
     #k + 2)) / hz * np1_dzdZ(i, j, k) - 0.1D1 * e * np1_qden(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
