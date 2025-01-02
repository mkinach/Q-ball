      subroutine u_eden_poly(np1_At,np1_Ax,np1_Ay,np1_Az,np1_Bx,np1_By,n
     &p1_Bz,np1_dxdX,np1_dydY,np1_dzdZ,np1_phi1,np1_phi2,np1_pi1,np1_pi2
     &,Nx,Ny,Nz,e,g,h,hx,hy,hz,m,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 e
      real*8 g
      real*8 h
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 m
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
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = e ** 2 * (np1_Az(i, j, k) ** 2 + np1_At(i, j, k) ** 2 + np1_A
     #x(i, j, k) ** 2 + np1_Ay(i, j, k) ** 2) * np1_phi2(i, j, k) ** 2 +
     # e ** 2 * (np1_Az(i, j, k) ** 2 + np1_At(i, j, k) ** 2 + np1_Ax(i,
     # j, k) ** 2 + np1_Ay(i, j, k) ** 2) * np1_phi1(i, j, k) ** 2 - 0.2
     #D1 * np1_phi1(i, j, k) * np1_pi2(i, j, k) * np1_At(i, j, k) * e + 
     #0.2D1 * np1_phi2(i, j, k) * np1_pi1(i, j, k) * np1_At(i, j, k) * e
     # - 0.8333333333333333D-1 * np1_Bz(i, j, k) * (np1_At(i, j, k - 2) 
     #- 0.8D1 * np1_At(i, j, k - 1) + 0.8D1 * np1_At(i, j, k + 1) - 0.1D
     #1 * np1_At(i, j, k + 2)) / hz * np1_dzdZ(i, j, k) - 0.833333333333
     #3333D-1 * np1_By(i, j, k) * (np1_At(i, j - 2, k) - 0.8D1 * np1_At(
     #i, j - 1, k) + 0.8D1 * np1_At(i, j + 1, k) - 0.1D1 * np1_At(i, j +
     # 2, k)) / hy * np1_dydY(i, j, k) - 0.8333333333333333D-1 * np1_Bx(
     #i, j, k) * (np1_At(i - 2, j, k) - 0.8D1 * np1_At(i - 1, j, k) + 0.
     #8D1 * np1_At(i + 1, j, k) - 0.1D1 * np1_At(i + 2, j, k)) / hx * np
     #1_dxdX(i, j, k) + 0.5000000000000000D0 * np1_By(i, j, k) ** 2 + 0.
     #5000000000000000D0 * np1_Bx(i, j, k) ** 2 + 0.5000000000000000D0 *
     # np1_Bz(i, j, k) ** 2 + np1_pi1(i, j, k) ** 2 + np1_pi2(i, j, k) *
     #* 2 - 0.5000000000000000D0 * g * (np1_phi1(i, j, k) ** 2 + np1_phi
     #2(i, j, k) ** 2) ** 2 + 0.3333333333333333D0 * h * (np1_phi1(i, j,
     # k) ** 2 + np1_phi2(i, j, k) ** 2) ** 3 + m ** 2 * (np1_phi1(i, j,
     # k) ** 2 + np1_phi2(i, j, k) ** 2) + 0.6944444444444444D-2 * (np1_
     #phi2(i - 2, j, k) - 0.8D1 * np1_phi2(i - 1, j, k) + 0.8D1 * np1_ph
     #i2(i + 1, j, k) - 0.1D1 * np1_phi2(i + 2, j, k)) ** 2 / hx ** 2 * 
     #np1_dxdX(i, j, k) ** 2 + 0.6944444444444444D-2 * (np1_phi2(i, j - 
     #2, k) - 0.8D1 * np1_phi2(i, j - 1, k) + 0.8D1 * np1_phi2(i, j + 1,
     # k) - 0.1D1 * np1_phi2(i, j + 2, k)) ** 2 / hy ** 2 * np1_dydY(i, 
     #j, k) ** 2 + 0.3472222222222222D-2 * (np1_Az(i - 2, j, k) - 0.8D1 
     #* np1_Az(i - 1, j, k) + 0.8D1 * np1_Az(i + 1, j, k) - 0.1D1 * np1_
     #Az(i + 2, j, k)) ** 2 / hx ** 2 * np1_dxdX(i, j, k) ** 2 + 0.34722
     #22222222222D-2 * (np1_Ay(i, j, k - 2) - 0.8D1 * np1_Ay(i, j, k - 1
     #) + 0.8D1 * np1_Ay(i, j, k + 1) - 0.1D1 * np1_Ay(i, j, k + 2)) ** 
     #2 / hz ** 2 * np1_dzdZ(i, j, k) ** 2 + 0.6944444444444444D-2 * (np
     #1_phi2(i, j, k - 2) - 0.8D1 * np1_phi2(i, j, k - 1) + 0.8D1 * np1_
     #phi2(i, j, k + 1) - 0.1D1 * np1_phi2(i, j, k + 2)) ** 2 / hz ** 2 
     #* np1_dzdZ(i, j, k) ** 2 + 0.3472222222222222D-2 * (np1_Az(i, j - 
     #2, k) - 0.8D1 * np1_Az(i, j - 1, k) + 0.8D1 * np1_Az(i, j + 1, k) 
     #- 0.1D1 * np1_Az(i, j + 2, k)) ** 2 / hy ** 2 * np1_dydY(i, j, k) 
     #** 2 + 0.3472222222222222D-2 * (np1_At(i - 2, j, k) - 0.8D1 * np1_
     #At(i - 1, j, k) + 0.8D1 * np1_At(i + 1, j, k) - 0.1D1 * np1_At(i +
     # 2, j, k)) ** 2 / hx ** 2 * np1_dxdX(i, j, k) ** 2 + 0.69444444444
     #44444D-2 * (np1_phi1(i, j, k - 2) - 0.8D1 * np1_phi1(i, j, k - 1) 
     #+ 0.8D1 * np1_phi1(i, j, k + 1) - 0.1D1 * np1_phi1(i, j, k + 2)) *
     #* 2 / hz ** 2 * np1_dzdZ(i, j, k) ** 2 + 0.6944444444444444D-2 * (
     #np1_phi1(i, j - 2, k) - 0.8D1 * np1_phi1(i, j - 1, k) + 0.8D1 * np
     #1_phi1(i, j + 1, k) - 0.1D1 * np1_phi1(i, j + 2, k)) ** 2 / hy ** 
     #2 * np1_dydY(i, j, k) ** 2 + 0.6944444444444444D-2 * (np1_phi1(i -
     # 2, j, k) - 0.8D1 * np1_phi1(i - 1, j, k) + 0.8D1 * np1_phi1(i + 1
     #, j, k) - 0.1D1 * np1_phi1(i + 2, j, k)) ** 2 / hx ** 2 * np1_dxdX
     #(i, j, k) ** 2 + 0.3472222222222222D-2 * (np1_Ay(i - 2, j, k) - 0.
     #8D1 * np1_Ay(i - 1, j, k) + 0.8D1 * np1_Ay(i + 1, j, k) - 0.1D1 * 
     #np1_Ay(i + 2, j, k)) ** 2 / hx ** 2 * np1_dxdX(i, j, k) ** 2 + 0.3
     #472222222222222D-2 * (np1_Ax(i, j - 2, k) - 0.8D1 * np1_Ax(i, j - 
     #1, k) + 0.8D1 * np1_Ax(i, j + 1, k) - 0.1D1 * np1_Ax(i, j + 2, k))
     # ** 2 / hy ** 2 * np1_dydY(i, j, k) ** 2 + 0.3472222222222222D-2 *
     # (np1_At(i, j - 2, k) - 0.8D1 * np1_At(i, j - 1, k) + 0.8D1 * np1_
     #At(i, j + 1, k) - 0.1D1 * np1_At(i, j + 2, k)) ** 2 / hy ** 2 * np
     #1_dydY(i, j, k) ** 2 + 0.3472222222222222D-2 * (np1_At(i, j, k - 2
     #) - 0.8D1 * np1_At(i, j, k - 1) + 0.8D1 * np1_At(i, j, k + 1) - 0.
     #1D1 * np1_At(i, j, k + 2)) ** 2 / hz ** 2 * np1_dzdZ(i, j, k) ** 2
     # + 0.3472222222222222D-2 * (np1_Ax(i, j, k - 2) - 0.8D1 * np1_Ax(i
     #, j, k - 1) + 0.8D1 * np1_Ax(i, j, k + 1) - 0.1D1 * np1_Ax(i, j, k
     # + 2)) ** 2 / hz ** 2 * np1_dzdZ(i, j, k) ** 2 + 0.166666666666666
     #7D0 * np1_phi2(i, j, k) * (np1_phi1(i, j - 2, k) - 0.8D1 * np1_phi
     #1(i, j - 1, k) + 0.8D1 * np1_phi1(i, j + 1, k) - 0.1D1 * np1_phi1(
     #i, j + 2, k)) / hy * np1_dydY(i, j, k) * np1_Ay(i, j, k) * e - 0.1
     #666666666666667D0 * np1_phi1(i, j, k) * (np1_phi2(i - 2, j, k) - 0
     #.8D1 * np1_phi2(i - 1, j, k) + 0.8D1 * np1_phi2(i + 1, j, k) - 0.1
     #D1 * np1_phi2(i + 2, j, k)) / hx * np1_dxdX(i, j, k) * np1_Ax(i, j
     #, k) * e + 0.1666666666666667D0 * np1_phi2(i, j, k) * (np1_phi1(i 
     #- 2, j, k) - 0.8D1 * np1_phi1(i - 1, j, k) + 0.8D1 * np1_phi1(i + 
     #1, j, k) - 0.1D1 * np1_phi1(i + 2, j, k)) / hx * np1_dxdX(i, j, k)
     # * np1_Ax(i, j, k) * e - 0.1666666666666667D0 * np1_phi1(i, j, k) 
     #* (np1_phi2(i, j, k - 2) - 0.8D1 * np1_phi2(i, j, k - 1) + 0.8D1 *
     # np1_phi2(i, j, k + 1) - 0.1D1 * np1_phi2(i, j, k + 2)) / hz * np1
     #_dzdZ(i, j, k) * np1_Az(i, j, k) * e + 0.1666666666666667D0 * np1_
     #phi2(i, j, k) * (np1_phi1(i, j, k - 2) - 0.8D1 * np1_phi1(i, j, k 
     #- 1) + 0.8D1 * np1_phi1(i, j, k + 1) - 0.1D1 * np1_phi1(i, j, k + 
     #2)) / hz * np1_dzdZ(i, j, k) * np1_Az(i, j, k) * e - 0.16666666666
     #66667D0 * np1_phi1(i, j, k) * (np1_phi2(i, j - 2, k) - 0.8D1 * np1
     #_phi2(i, j - 1, k) + 0.8D1 * np1_phi2(i, j + 1, k) - 0.1D1 * np1_p
     #hi2(i, j + 2, k)) / hy * np1_dydY(i, j, k) * np1_Ay(i, j, k) * e -
     # 0.6944444444444444D-2 * (np1_Ay(i - 2, j, k) - 0.8D1 * np1_Ay(i -
     # 1, j, k) + 0.8D1 * np1_Ay(i + 1, j, k) - 0.1D1 * np1_Ay(i + 2, j,
     # k)) / hx * np1_dxdX(i, j, k) * (np1_Ax(i, j - 2, k) - 0.8D1 * np1
     #_Ax(i, j - 1, k) + 0.8D1 * np1_Ax(i, j + 1, k) - 0.1D1 * np1_Ax(i,
     # j + 2, k)) / hy * np1_dydY(i, j, k) - 0.6944444444444444D-2 * (np
     #1_Az(i - 2, j, k) - 0.8D1 * np1_Az(i - 1, j, k) + 0.8D1 * np1_Az(i
     # + 1, j, k) - 0.1D1 * np1_Az(i + 2, j, k)) / hx * np1_dxdX(i, j, k
     #) * (np1_Ax(i, j, k - 2) - 0.8D1 * np1_Ax(i, j, k - 1) + 0.8D1 * n
     #p1_Ax(i, j, k + 1) - 0.1D1 * np1_Ax(i, j, k + 2)) / hz * np1_dzdZ(
     #i, j, k) - 0.6944444444444444D-2 * (np1_Az(i, j - 2, k) - 0.8D1 * 
     #np1_Az(i, j - 1, k) + 0.8D1 * np1_Az(i, j + 1, k) - 0.1D1 * np1_Az
     #(i, j + 2, k)) / hy * np1_dydY(i, j, k) * (np1_Ay(i, j, k - 2) - 0
     #.8D1 * np1_Ay(i, j, k - 1) + 0.8D1 * np1_Ay(i, j, k + 1) - 0.1D1 *
     # np1_Ay(i, j, k + 2)) / hz * np1_dzdZ(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
