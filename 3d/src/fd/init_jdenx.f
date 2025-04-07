      subroutine init_jdenx(n_At,n_Ax,n_Ay,n_Az,n_Bx,n_By,n_Bz,n_dxdX,n_
     &dydY,n_dzdZ,n_phi1,n_phi2,n_pi1,n_pi2,y,z,Nx,Ny,Nz,c,d,e,hx,hy,hz,
     &res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 c
      real*8 d
      real*8 e
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_Ax(Nx,Ny,Nz)
      real*8 n_Ay(Nx,Ny,Nz)
      real*8 n_Az(Nx,Ny,Nz)
      real*8 n_Bx(Nx,Ny,Nz)
      real*8 n_By(Nx,Ny,Nz)
      real*8 n_Bz(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 n_pi1(Nx,Ny,Nz)
      real*8 n_pi2(Nx,Ny,Nz)
      real*8 y(Ny)
      real*8 z(Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (0.8333333333333333D-1 * (d * exp(c *
     # y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) * (n_Ax(i, j, k - 2) 
     #- 0.8D1 * n_Ax(i, j, k - 1) + 0.8D1 * n_Ax(i, j, k + 1) - 0.1D1 * 
     #n_Ax(i, j, k + 2)) / hz * n_dzdZ(i, j, k) - 0.8333333333333333D-1 
     #* (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) * (n_Az
     #(i - 2, j, k) - 0.8D1 * n_Az(i - 1, j, k) + 0.8D1 * n_Az(i + 1, j,
     # k) - 0.1D1 * n_Az(i + 2, j, k)) / hx * n_dxdX(i, j, k) - 0.833333
     #3333333333D-1 * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * 
     #z(k))) * (n_Ax(i, j - 2, k) - 0.8D1 * n_Ax(i, j - 1, k) + 0.8D1 * 
     #n_Ax(i, j + 1, k) - 0.1D1 * n_Ax(i, j + 2, k)) / hy * n_dydY(i, j,
     # k) + 0.8333333333333333D-1 * (d * exp(c * z(k)) - 0.1D1 * d * exp
     #(-0.1D1 * c * z(k))) * (n_Ay(i - 2, j, k) - 0.8D1 * n_Ay(i - 1, j,
     # k) + 0.8D1 * n_Ay(i + 1, j, k) - 0.1D1 * n_Ay(i + 2, j, k)) / hx 
     #* n_dxdX(i, j, k)) * (n_At(i - 2, j, k) - 0.8D1 * n_At(i - 1, j, k
     #) + 0.8D1 * n_At(i + 1, j, k) - 0.1D1 * n_At(i + 2, j, k)) / hx * 
     #n_dxdX(i, j, k) + 0.8333333333333333D-1 * (-0.1D1 * n_By(i, j, k) 
     #* (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) - 0.1D1
     # * n_Bz(i, j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c
     # * z(k))) + 0.8333333333333333D-1 * (d * exp(c * y(j)) - 0.1D1 * d
     # * exp(-0.1D1 * c * y(j))) * (n_At(i, j - 2, k) - 0.8D1 * n_At(i, 
     #j - 1, k) + 0.8D1 * n_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, k))
     # / hy * n_dydY(i, j, k) + 0.8333333333333333D-1 * (d * exp(c * z(k
     #)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (n_At(i, j, k - 2) - 0.
     #8D1 * n_At(i, j, k - 1) + 0.8D1 * n_At(i, j, k + 1) - 0.1D1 * n_At
     #(i, j, k + 2)) / hz * n_dzdZ(i, j, k)) * (n_Ay(i, j, k - 2) - 0.8D
     #1 * n_Ay(i, j, k - 1) + 0.8D1 * n_Ay(i, j, k + 1) - 0.1D1 * n_Ay(i
     #, j, k + 2)) / hz * n_dzdZ(i, j, k) + 0.8333333333333333D-1 * (n_B
     #y(i, j, k) * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j
     #))) + n_Bz(i, j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 
     #* c * z(k))) - 0.8333333333333333D-1 * (d * exp(c * y(j)) - 0.1D1 
     #* d * exp(-0.1D1 * c * y(j))) * (n_At(i, j - 2, k) - 0.8D1 * n_At(
     #i, j - 1, k) + 0.8D1 * n_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, 
     #k)) / hy * n_dydY(i, j, k) - 0.8333333333333333D-1 * (d * exp(c * 
     #z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (n_At(i, j, k - 2) -
     # 0.8D1 * n_At(i, j, k - 1) + 0.8D1 * n_At(i, j, k + 1) - 0.1D1 * n
     #_At(i, j, k + 2)) / hz * n_dzdZ(i, j, k)) * (n_Az(i, j - 2, k) - 0
     #.8D1 * n_Az(i, j - 1, k) + 0.8D1 * n_Az(i, j + 1, k) - 0.1D1 * n_A
     #z(i, j + 2, k)) / hy * n_dydY(i, j, k) + 0.1666666666666667D0 * (d
     # * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (n_phi2(i
     #, j, k) * e * n_At(i, j, k) + n_pi1(i, j, k)) * (n_phi1(i, j - 2, 
     #k) - 0.8D1 * n_phi1(i, j - 1, k) + 0.8D1 * n_phi1(i, j + 1, k) - 0
     #.1D1 * n_phi1(i, j + 2, k)) / hy * n_dydY(i, j, k) - 0.16666666666
     #66667D0 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j)))
     # * (n_phi2(i, j, k) * e * n_At(i, j, k) + n_pi1(i, j, k)) * (n_phi
     #1(i, j, k - 2) - 0.8D1 * n_phi1(i, j, k - 1) + 0.8D1 * n_phi1(i, j
     #, k + 1) - 0.1D1 * n_phi1(i, j, k + 2)) / hz * n_dzdZ(i, j, k) + 0
     #.1666666666666667D0 * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 
     #* c * z(k))) * (-0.1D1 * n_phi1(i, j, k) * e * n_At(i, j, k) + n_p
     #i2(i, j, k)) * (n_phi2(i, j - 2, k) - 0.8D1 * n_phi2(i, j - 1, k) 
     #+ 0.8D1 * n_phi2(i, j + 1, k) - 0.1D1 * n_phi2(i, j + 2, k)) / hy 
     #* n_dydY(i, j, k) - 0.1666666666666667D0 * (d * exp(c * y(j)) - 0.
     #1D1 * d * exp(-0.1D1 * c * y(j))) * (-0.1D1 * n_phi1(i, j, k) * e 
     #* n_At(i, j, k) + n_pi2(i, j, k)) * (n_phi2(i, j, k - 2) - 0.8D1 *
     # n_phi2(i, j, k - 1) + 0.8D1 * n_phi2(i, j, k + 1) - 0.1D1 * n_phi
     #2(i, j, k + 2)) / hz * n_dzdZ(i, j, k) + 0.8333333333333333D-1 * n
     #_Bx(i, j, k) * (n_Ax(i, j - 2, k) - 0.8D1 * n_Ax(i, j - 1, k) + 0.
     #8D1 * n_Ax(i, j + 1, k) - 0.1D1 * n_Ax(i, j + 2, k)) / hy * n_dydY
     #(i, j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k)
     #)) - 0.8333333333333333D-1 * n_Bx(i, j, k) * (n_Ax(i, j, k - 2) - 
     #0.8D1 * n_Ax(i, j, k - 1) + 0.8D1 * n_Ax(i, j, k + 1) - 0.1D1 * n_
     #Ax(i, j, k + 2)) / hz * n_dzdZ(i, j, k) * (d * exp(c * y(j)) - 0.1
     #D1 * d * exp(-0.1D1 * c * y(j))) - 0.8333333333333333D-1 * n_Bx(i,
     # j, k) * (n_Ay(i - 2, j, k) - 0.8D1 * n_Ay(i - 1, j, k) + 0.8D1 * 
     #n_Ay(i + 1, j, k) - 0.1D1 * n_Ay(i + 2, j, k)) / hx * n_dxdX(i, j,
     # k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) + 0
     #.8333333333333333D-1 * (n_Az(i - 2, j, k) - 0.8D1 * n_Az(i - 1, j,
     # k) + 0.8D1 * n_Az(i + 1, j, k) - 0.1D1 * n_Az(i + 2, j, k)) / hx 
     #* n_dxdX(i, j, k) * n_Bx(i, j, k) * (d * exp(c * y(j)) - 0.1D1 * d
     # * exp(-0.1D1 * c * y(j))) - 0.2D1 * ((d * exp(c * y(j)) - 0.1D1 *
     # d * exp(-0.1D1 * c * y(j))) * n_Az(i, j, k) - 0.1D1 * (d * exp(c 
     #* z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * n_Ay(i, j, k)) * e
     # * (e * (n_phi1(i, j, k) ** 2 + n_phi2(i, j, k) ** 2) * n_At(i, j,
     # k) + n_phi2(i, j, k) * n_pi1(i, j, k) - 0.1D1 * n_phi1(i, j, k) *
     # n_pi2(i, j, k))
      res(i,j,k)=qb
      end do
      end do
      end do
      END
