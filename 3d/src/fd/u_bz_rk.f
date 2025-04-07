      subroutine u_bz_rk(n_Az,n_dxdX,n_dxdX2,n_dydY,n_dydY2,n_dzdZ,n_dzd
     &Z2,n_phi1,n_phi2,np1_Bz,x,y,z,Nx,Ny,Nz,e,hx,hy,hz,myzero,phys_bdy,
     &res)
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
      real*8 myzero
      real*8 n_Az(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dxdX2(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dydY2(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 n_dzdZ2(Nx,Ny,Nz)
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 np1_Bz(Nx,Ny,Nz)
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = -0.1D1 * (0.2D1 * n_Az(i, j, k) * n_phi1(i, j, k) ** 2 + 0.2D
     #1 * n_Az(i, j, k) * n_phi2(i, j, k) ** 2) * e ** 2 - 0.1D1 * (-0.1
     #666666666666667D0 * n_phi1(i, j, k) * (n_phi2(i, j, k - 2) - 0.8D1
     # * n_phi2(i, j, k - 1) + 0.8D1 * n_phi2(i, j, k + 1) - 0.1D1 * n_p
     #hi2(i, j, k + 2)) / hz * n_dzdZ(i, j, k) + 0.1666666666666667D0 * 
     #n_phi2(i, j, k) * (n_phi1(i, j, k - 2) - 0.8D1 * n_phi1(i, j, k - 
     #1) + 0.8D1 * n_phi1(i, j, k + 1) - 0.1D1 * n_phi1(i, j, k + 2)) / 
     #hz * n_dzdZ(i, j, k)) * e + 0.8333333333333333D-1 * (-0.1D1 * n_Az
     #(i, j, k - 2) + 0.16D2 * n_Az(i, j, k - 1) - 0.30D2 * n_Az(i, j, k
     #) + 0.16D2 * n_Az(i, j, k + 1) - 0.1D1 * n_Az(i, j, k + 2)) / hz *
     #* 2 * n_dzdZ(i, j, k) ** 2 + 0.8333333333333333D-1 * (n_Az(i, j, k
     # - 2) - 0.8D1 * n_Az(i, j, k - 1) + 0.8D1 * n_Az(i, j, k + 1) - 0.
     #1D1 * n_Az(i, j, k + 2)) / hz * n_dzdZ2(i, j, k) + 0.8333333333333
     #333D-1 * (-0.1D1 * n_Az(i - 2, j, k) + 0.16D2 * n_Az(i - 1, j, k) 
     #- 0.30D2 * n_Az(i, j, k) + 0.16D2 * n_Az(i + 1, j, k) - 0.1D1 * n_
     #Az(i + 2, j, k)) / hx ** 2 * n_dxdX(i, j, k) ** 2 + 0.833333333333
     #3333D-1 * (n_Az(i - 2, j, k) - 0.8D1 * n_Az(i - 1, j, k) + 0.8D1 *
     # n_Az(i + 1, j, k) - 0.1D1 * n_Az(i + 2, j, k)) / hx * n_dxdX2(i, 
     #j, k) + 0.8333333333333333D-1 * (-0.1D1 * n_Az(i, j - 2, k) + 0.16
     #D2 * n_Az(i, j - 1, k) - 0.30D2 * n_Az(i, j, k) + 0.16D2 * n_Az(i,
     # j + 1, k) - 0.1D1 * n_Az(i, j + 2, k)) / hy ** 2 * n_dydY(i, j, k
     #) ** 2 + 0.8333333333333333D-1 * (n_Az(i, j - 2, k) - 0.8D1 * n_Az
     #(i, j - 1, k) + 0.8D1 * n_Az(i, j + 1, k) - 0.1D1 * n_Az(i, j + 2,
     # k)) / hy * n_dydY2(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=1, Nx, 1
      do j=1, 1, 1
      do k=1, Nz, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=1, Nx, 1
      do j=Ny, Ny, 1
      do k=1, Nz, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(5) .eq. 1) then
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, 1, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(6) .eq. 1) then
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=Nz, Nz, 1
      qb = np1_Bz(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      END
