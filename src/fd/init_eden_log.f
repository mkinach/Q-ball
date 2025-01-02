      subroutine init_eden_log(n_At,n_Ax,n_Ay,n_Az,n_Bx,n_By,n_Bz,n_dxdX
     &,n_dydY,n_dzdZ,n_phi1,n_phi2,n_pi1,n_pi2,Nx,Ny,Nz,beta,e,hx,hy,hz,
     &mu,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 beta
      real*8 e
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 mu
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
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = e ** 2 * (n_Az(i, j, k) ** 2 + n_At(i, j, k) ** 2 + n_Ax(i, j
     #, k) ** 2 + n_Ay(i, j, k) ** 2) * n_phi2(i, j, k) ** 2 + e ** 2 * 
     #(n_Az(i, j, k) ** 2 + n_At(i, j, k) ** 2 + n_Ax(i, j, k) ** 2 + n_
     #Ay(i, j, k) ** 2) * n_phi1(i, j, k) ** 2 + 0.1666666666666667D0 * 
     #n_phi2(i, j, k) * (n_phi1(i, j - 2, k) - 0.8D1 * n_phi1(i, j - 1, 
     #k) + 0.8D1 * n_phi1(i, j + 1, k) - 0.1D1 * n_phi1(i, j + 2, k)) / 
     #hy * n_dydY(i, j, k) * n_Ay(i, j, k) * e - 0.1666666666666667D0 * 
     #n_phi1(i, j, k) * (n_phi2(i - 2, j, k) - 0.8D1 * n_phi2(i - 1, j, 
     #k) + 0.8D1 * n_phi2(i + 1, j, k) - 0.1D1 * n_phi2(i + 2, j, k)) / 
     #hx * n_dxdX(i, j, k) * n_Ax(i, j, k) * e + 0.1666666666666667D0 * 
     #n_phi2(i, j, k) * (n_phi1(i - 2, j, k) - 0.8D1 * n_phi1(i - 1, j, 
     #k) + 0.8D1 * n_phi1(i + 1, j, k) - 0.1D1 * n_phi1(i + 2, j, k)) / 
     #hx * n_dxdX(i, j, k) * n_Ax(i, j, k) * e - 0.1666666666666667D0 * 
     #n_phi1(i, j, k) * (n_phi2(i, j, k - 2) - 0.8D1 * n_phi2(i, j, k - 
     #1) + 0.8D1 * n_phi2(i, j, k + 1) - 0.1D1 * n_phi2(i, j, k + 2)) / 
     #hz * n_dzdZ(i, j, k) * n_Az(i, j, k) * e + 0.1666666666666667D0 * 
     #n_phi2(i, j, k) * (n_phi1(i, j, k - 2) - 0.8D1 * n_phi1(i, j, k - 
     #1) + 0.8D1 * n_phi1(i, j, k + 1) - 0.1D1 * n_phi1(i, j, k + 2)) / 
     #hz * n_dzdZ(i, j, k) * n_Az(i, j, k) * e - 0.1666666666666667D0 * 
     #n_phi1(i, j, k) * (n_phi2(i, j - 2, k) - 0.8D1 * n_phi2(i, j - 1, 
     #k) + 0.8D1 * n_phi2(i, j + 1, k) - 0.1D1 * n_phi2(i, j + 2, k)) / 
     #hy * n_dydY(i, j, k) * n_Ay(i, j, k) * e - 0.6944444444444444D-2 *
     # (n_Ay(i - 2, j, k) - 0.8D1 * n_Ay(i - 1, j, k) + 0.8D1 * n_Ay(i +
     # 1, j, k) - 0.1D1 * n_Ay(i + 2, j, k)) / hx * n_dxdX(i, j, k) * (n
     #_Ax(i, j - 2, k) - 0.8D1 * n_Ax(i, j - 1, k) + 0.8D1 * n_Ax(i, j +
     # 1, k) - 0.1D1 * n_Ax(i, j + 2, k)) / hy * n_dydY(i, j, k) - 0.694
     #4444444444444D-2 * (n_Az(i - 2, j, k) - 0.8D1 * n_Az(i - 1, j, k) 
     #+ 0.8D1 * n_Az(i + 1, j, k) - 0.1D1 * n_Az(i + 2, j, k)) / hx * n_
     #dxdX(i, j, k) * (n_Ax(i, j, k - 2) - 0.8D1 * n_Ax(i, j, k - 1) + 0
     #.8D1 * n_Ax(i, j, k + 1) - 0.1D1 * n_Ax(i, j, k + 2)) / hz * n_dzd
     #Z(i, j, k) - 0.6944444444444444D-2 * (n_Az(i, j - 2, k) - 0.8D1 * 
     #n_Az(i, j - 1, k) + 0.8D1 * n_Az(i, j + 1, k) - 0.1D1 * n_Az(i, j 
     #+ 2, k)) / hy * n_dydY(i, j, k) * (n_Ay(i, j, k - 2) - 0.8D1 * n_A
     #y(i, j, k - 1) + 0.8D1 * n_Ay(i, j, k + 1) - 0.1D1 * n_Ay(i, j, k 
     #+ 2)) / hz * n_dzdZ(i, j, k) - 0.2D1 * n_phi1(i, j, k) * n_pi2(i, 
     #j, k) * n_At(i, j, k) * e + 0.2D1 * n_phi2(i, j, k) * n_pi1(i, j, 
     #k) * n_At(i, j, k) * e + 0.5000000000000000D0 * n_Bx(i, j, k) ** 2
     # + 0.5000000000000000D0 * n_Bz(i, j, k) ** 2 + 0.6944444444444444D
     #-2 * (n_phi2(i - 2, j, k) - 0.8D1 * n_phi2(i - 1, j, k) + 0.8D1 * 
     #n_phi2(i + 1, j, k) - 0.1D1 * n_phi2(i + 2, j, k)) ** 2 / hx ** 2 
     #* n_dxdX(i, j, k) ** 2 + 0.6944444444444444D-2 * (n_phi2(i, j - 2,
     # k) - 0.8D1 * n_phi2(i, j - 1, k) + 0.8D1 * n_phi2(i, j + 1, k) - 
     #0.1D1 * n_phi2(i, j + 2, k)) ** 2 / hy ** 2 * n_dydY(i, j, k) ** 2
     # + 0.3472222222222222D-2 * (n_Az(i - 2, j, k) - 0.8D1 * n_Az(i - 1
     #, j, k) + 0.8D1 * n_Az(i + 1, j, k) - 0.1D1 * n_Az(i + 2, j, k)) *
     #* 2 / hx ** 2 * n_dxdX(i, j, k) ** 2 + 0.3472222222222222D-2 * (n_
     #Ay(i, j, k - 2) - 0.8D1 * n_Ay(i, j, k - 1) + 0.8D1 * n_Ay(i, j, k
     # + 1) - 0.1D1 * n_Ay(i, j, k + 2)) ** 2 / hz ** 2 * n_dzdZ(i, j, k
     #) ** 2 + 0.6944444444444444D-2 * (n_phi2(i, j, k - 2) - 0.8D1 * n_
     #phi2(i, j, k - 1) + 0.8D1 * n_phi2(i, j, k + 1) - 0.1D1 * n_phi2(i
     #, j, k + 2)) ** 2 / hz ** 2 * n_dzdZ(i, j, k) ** 2 + 0.34722222222
     #22222D-2 * (n_Az(i, j - 2, k) - 0.8D1 * n_Az(i, j - 1, k) + 0.8D1 
     #* n_Az(i, j + 1, k) - 0.1D1 * n_Az(i, j + 2, k)) ** 2 / hy ** 2 * 
     #n_dydY(i, j, k) ** 2 + 0.3472222222222222D-2 * (n_At(i - 2, j, k) 
     #- 0.8D1 * n_At(i - 1, j, k) + 0.8D1 * n_At(i + 1, j, k) - 0.1D1 * 
     #n_At(i + 2, j, k)) ** 2 / hx ** 2 * n_dxdX(i, j, k) ** 2 + 0.69444
     #44444444444D-2 * (n_phi1(i, j, k - 2) - 0.8D1 * n_phi1(i, j, k - 1
     #) + 0.8D1 * n_phi1(i, j, k + 1) - 0.1D1 * n_phi1(i, j, k + 2)) ** 
     #2 / hz ** 2 * n_dzdZ(i, j, k) ** 2 + 0.6944444444444444D-2 * (n_ph
     #i1(i, j - 2, k) - 0.8D1 * n_phi1(i, j - 1, k) + 0.8D1 * n_phi1(i, 
     #j + 1, k) - 0.1D1 * n_phi1(i, j + 2, k)) ** 2 / hy ** 2 * n_dydY(i
     #, j, k) ** 2 + 0.6944444444444444D-2 * (n_phi1(i - 2, j, k) - 0.8D
     #1 * n_phi1(i - 1, j, k) + 0.8D1 * n_phi1(i + 1, j, k) - 0.1D1 * n_
     #phi1(i + 2, j, k)) ** 2 / hx ** 2 * n_dxdX(i, j, k) ** 2 + 0.34722
     #22222222222D-2 * (n_Ay(i - 2, j, k) - 0.8D1 * n_Ay(i - 1, j, k) + 
     #0.8D1 * n_Ay(i + 1, j, k) - 0.1D1 * n_Ay(i + 2, j, k)) ** 2 / hx *
     #* 2 * n_dxdX(i, j, k) ** 2 + 0.3472222222222222D-2 * (n_Ax(i, j - 
     #2, k) - 0.8D1 * n_Ax(i, j - 1, k) + 0.8D1 * n_Ax(i, j + 1, k) - 0.
     #1D1 * n_Ax(i, j + 2, k)) ** 2 / hy ** 2 * n_dydY(i, j, k) ** 2 + 0
     #.3472222222222222D-2 * (n_At(i, j - 2, k) - 0.8D1 * n_At(i, j - 1,
     # k) + 0.8D1 * n_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, k)) ** 2 
     #/ hy ** 2 * n_dydY(i, j, k) ** 2 + 0.3472222222222222D-2 * (n_At(i
     #, j, k - 2) - 0.8D1 * n_At(i, j, k - 1) + 0.8D1 * n_At(i, j, k + 1
     #) - 0.1D1 * n_At(i, j, k + 2)) ** 2 / hz ** 2 * n_dzdZ(i, j, k) **
     # 2 + 0.3472222222222222D-2 * (n_Ax(i, j, k - 2) - 0.8D1 * n_Ax(i, 
     #j, k - 1) + 0.8D1 * n_Ax(i, j, k + 1) - 0.1D1 * n_Ax(i, j, k + 2))
     # ** 2 / hz ** 2 * n_dzdZ(i, j, k) ** 2 - 0.1D1 * mu ** 2 * (n_phi1
     #(i, j, k) ** 2 + n_phi2(i, j, k) ** 2) * log(beta ** 2 * (n_phi1(i
     #, j, k) ** 2 + n_phi2(i, j, k) ** 2) + 0.1000000000000000D-99) + n
     #_pi2(i, j, k) ** 2 + n_pi1(i, j, k) ** 2 + 0.5000000000000000D0 * 
     #n_By(i, j, k) ** 2 - 0.8333333333333333D-1 * n_Bz(i, j, k) * (n_At
     #(i, j, k - 2) - 0.8D1 * n_At(i, j, k - 1) + 0.8D1 * n_At(i, j, k +
     # 1) - 0.1D1 * n_At(i, j, k + 2)) / hz * n_dzdZ(i, j, k) - 0.833333
     #3333333333D-1 * n_By(i, j, k) * (n_At(i, j - 2, k) - 0.8D1 * n_At(
     #i, j - 1, k) + 0.8D1 * n_At(i, j + 1, k) - 0.1D1 * n_At(i, j + 2, 
     #k)) / hy * n_dydY(i, j, k) - 0.8333333333333333D-1 * n_Bx(i, j, k)
     # * (n_At(i - 2, j, k) - 0.8D1 * n_At(i - 1, j, k) + 0.8D1 * n_At(i
     # + 1, j, k) - 0.1D1 * n_At(i + 2, j, k)) / hx * n_dxdX(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
