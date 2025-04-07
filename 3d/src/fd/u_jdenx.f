      subroutine u_jdenx(np1_At,np1_Ax,np1_Ay,np1_Az,np1_Bx,np1_By,np1_B
     &z,np1_dxdX,np1_dydY,np1_dzdZ,np1_phi1,np1_phi2,np1_pi1,np1_pi2,y,z
     &,Nx,Ny,Nz,c,d,e,hx,hy,hz,phys_bdy,res)
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
      real*8 np1_At(Nx,Ny,Nz)
      real*8 np1_Ax(Nx,Ny,Nz)
      real*8 np1_Ay(Nx,Ny,Nz)
      real*8 np1_Az(Nx,Ny,Nz)
      real*8 np1_Bx(Nx,Ny,Nz)
      real*8 np1_By(Nx,Ny,Nz)
      real*8 np1_Bz(Nx,Ny,Nz)
      real*8 np1_dxdX(Nx,Ny,Nz)
      real*8 np1_dydY(Nx,Ny,Nz)
      real*8 np1_dzdZ(Nx,Ny,Nz)
      real*8 np1_phi1(Nx,Ny,Nz)
      real*8 np1_phi2(Nx,Ny,Nz)
      real*8 np1_pi1(Nx,Ny,Nz)
      real*8 np1_pi2(Nx,Ny,Nz)
      real*8 y(Ny)
      real*8 z(Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (0.8333333333333333D-1 * (d * exp(c *
     # y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) * (np1_Ax(i, j, k - 2
     #) - 0.8D1 * np1_Ax(i, j, k - 1) + 0.8D1 * np1_Ax(i, j, k + 1) - 0.
     #1D1 * np1_Ax(i, j, k + 2)) / hz * np1_dzdZ(i, j, k) - 0.8333333333
     #333333D-1 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j)
     #)) * (np1_Az(i - 2, j, k) - 0.8D1 * np1_Az(i - 1, j, k) + 0.8D1 * 
     #np1_Az(i + 1, j, k) - 0.1D1 * np1_Az(i + 2, j, k)) / hx * np1_dxdX
     #(i, j, k) - 0.8333333333333333D-1 * (d * exp(c * z(k)) - 0.1D1 * d
     # * exp(-0.1D1 * c * z(k))) * (np1_Ax(i, j - 2, k) - 0.8D1 * np1_Ax
     #(i, j - 1, k) + 0.8D1 * np1_Ax(i, j + 1, k) - 0.1D1 * np1_Ax(i, j 
     #+ 2, k)) / hy * np1_dydY(i, j, k) + 0.8333333333333333D-1 * (d * e
     #xp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (np1_Ay(i - 2
     #, j, k) - 0.8D1 * np1_Ay(i - 1, j, k) + 0.8D1 * np1_Ay(i + 1, j, k
     #) - 0.1D1 * np1_Ay(i + 2, j, k)) / hx * np1_dxdX(i, j, k)) * (np1_
     #At(i - 2, j, k) - 0.8D1 * np1_At(i - 1, j, k) + 0.8D1 * np1_At(i +
     # 1, j, k) - 0.1D1 * np1_At(i + 2, j, k)) / hx * np1_dxdX(i, j, k) 
     #+ 0.8333333333333333D-1 * (-0.1D1 * np1_By(i, j, k) * (d * exp(c *
     # y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) - 0.1D1 * np1_Bz(i, j
     #, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) + 
     #0.8333333333333333D-1 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D
     #1 * c * y(j))) * (np1_At(i, j - 2, k) - 0.8D1 * np1_At(i, j - 1, k
     #) + 0.8D1 * np1_At(i, j + 1, k) - 0.1D1 * np1_At(i, j + 2, k)) / h
     #y * np1_dydY(i, j, k) + 0.8333333333333333D-1 * (d * exp(c * z(k))
     # - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (np1_At(i, j, k - 2) - 0.
     #8D1 * np1_At(i, j, k - 1) + 0.8D1 * np1_At(i, j, k + 1) - 0.1D1 * 
     #np1_At(i, j, k + 2)) / hz * np1_dzdZ(i, j, k)) * (np1_Ay(i, j, k -
     # 2) - 0.8D1 * np1_Ay(i, j, k - 1) + 0.8D1 * np1_Ay(i, j, k + 1) - 
     #0.1D1 * np1_Ay(i, j, k + 2)) / hz * np1_dzdZ(i, j, k) + 0.83333333
     #33333333D-1 * (np1_By(i, j, k) * (d * exp(c * y(j)) - 0.1D1 * d * 
     #exp(-0.1D1 * c * y(j))) + np1_Bz(i, j, k) * (d * exp(c * z(k)) - 0
     #.1D1 * d * exp(-0.1D1 * c * z(k))) - 0.8333333333333333D-1 * (d * 
     #exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) * (np1_At(i, j
     # - 2, k) - 0.8D1 * np1_At(i, j - 1, k) + 0.8D1 * np1_At(i, j + 1, 
     #k) - 0.1D1 * np1_At(i, j + 2, k)) / hy * np1_dydY(i, j, k) - 0.833
     #3333333333333D-1 * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c
     # * z(k))) * (np1_At(i, j, k - 2) - 0.8D1 * np1_At(i, j, k - 1) + 0
     #.8D1 * np1_At(i, j, k + 1) - 0.1D1 * np1_At(i, j, k + 2)) / hz * n
     #p1_dzdZ(i, j, k)) * (np1_Az(i, j - 2, k) - 0.8D1 * np1_Az(i, j - 1
     #, k) + 0.8D1 * np1_Az(i, j + 1, k) - 0.1D1 * np1_Az(i, j + 2, k)) 
     #/ hy * np1_dydY(i, j, k) + 0.1666666666666667D0 * (d * exp(c * z(k
     #)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) * (np1_phi2(i, j, k) * e 
     #* np1_At(i, j, k) + np1_pi1(i, j, k)) * (np1_phi1(i, j - 2, k) - 0
     #.8D1 * np1_phi1(i, j - 1, k) + 0.8D1 * np1_phi1(i, j + 1, k) - 0.1
     #D1 * np1_phi1(i, j + 2, k)) / hy * np1_dydY(i, j, k) - 0.166666666
     #6666667D0 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j)
     #)) * (np1_phi2(i, j, k) * e * np1_At(i, j, k) + np1_pi1(i, j, k)) 
     #* (np1_phi1(i, j, k - 2) - 0.8D1 * np1_phi1(i, j, k - 1) + 0.8D1 *
     # np1_phi1(i, j, k + 1) - 0.1D1 * np1_phi1(i, j, k + 2)) / hz * np1
     #_dzdZ(i, j, k) + 0.1666666666666667D0 * (d * exp(c * z(k)) - 0.1D1
     # * d * exp(-0.1D1 * c * z(k))) * (-0.1D1 * np1_phi1(i, j, k) * e *
     # np1_At(i, j, k) + np1_pi2(i, j, k)) * (np1_phi2(i, j - 2, k) - 0.
     #8D1 * np1_phi2(i, j - 1, k) + 0.8D1 * np1_phi2(i, j + 1, k) - 0.1D
     #1 * np1_phi2(i, j + 2, k)) / hy * np1_dydY(i, j, k) - 0.1666666666
     #666667D0 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))
     #) * (-0.1D1 * np1_phi1(i, j, k) * e * np1_At(i, j, k) + np1_pi2(i,
     # j, k)) * (np1_phi2(i, j, k - 2) - 0.8D1 * np1_phi2(i, j, k - 1) +
     # 0.8D1 * np1_phi2(i, j, k + 1) - 0.1D1 * np1_phi2(i, j, k + 2)) / 
     #hz * np1_dzdZ(i, j, k) + 0.8333333333333333D-1 * np1_Bx(i, j, k) *
     # (np1_Ax(i, j - 2, k) - 0.8D1 * np1_Ax(i, j - 1, k) + 0.8D1 * np1_
     #Ax(i, j + 1, k) - 0.1D1 * np1_Ax(i, j + 2, k)) / hy * np1_dydY(i, 
     #j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * c * z(k))) -
     # 0.8333333333333333D-1 * np1_Bx(i, j, k) * (np1_Ax(i, j, k - 2) - 
     #0.8D1 * np1_Ax(i, j, k - 1) + 0.8D1 * np1_Ax(i, j, k + 1) - 0.1D1 
     #* np1_Ax(i, j, k + 2)) / hz * np1_dzdZ(i, j, k) * (d * exp(c * y(j
     #)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) - 0.8333333333333333D-1 *
     # np1_Bx(i, j, k) * (np1_Ay(i - 2, j, k) - 0.8D1 * np1_Ay(i - 1, j,
     # k) + 0.8D1 * np1_Ay(i + 1, j, k) - 0.1D1 * np1_Ay(i + 2, j, k)) /
     # hx * np1_dxdX(i, j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.
     #1D1 * c * z(k))) + 0.8333333333333333D-1 * (np1_Az(i - 2, j, k) - 
     #0.8D1 * np1_Az(i - 1, j, k) + 0.8D1 * np1_Az(i + 1, j, k) - 0.1D1 
     #* np1_Az(i + 2, j, k)) / hx * np1_dxdX(i, j, k) * np1_Bx(i, j, k) 
     #* (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) - 0.2D1
     # * ((d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) * np1
     #_Az(i, j, k) - 0.1D1 * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1
     # * c * z(k))) * np1_Ay(i, j, k)) * e * (e * (np1_phi1(i, j, k) ** 
     #2 + np1_phi2(i, j, k) ** 2) * np1_At(i, j, k) + np1_phi2(i, j, k) 
     #* np1_pi1(i, j, k) - 0.1D1 * np1_phi1(i, j, k) * np1_pi2(i, j, k))
      res(i,j,k)=qb
      end do
      end do
      end do
      END
