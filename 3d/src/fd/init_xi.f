      subroutine init_xi(n_chi,n_dxdX,n_dydY,n_dzdZ,x,y,z,Nx,Ny,Nz,c,d,h
     &x,hy,hz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 c
      real*8 d
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 n_chi(Nx,Ny,Nz)
      real*8 n_dxdX(Nx,Ny,Nz)
      real*8 n_dydY(Nx,Ny,Nz)
      real*8 n_dzdZ(Nx,Ny,Nz)
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.10000000000D11 * (0.8333333333333333D-1 * (n_chi(i - 2, j, 
     #k) - 0.8D1 * n_chi(i - 1, j, k) + 0.8D1 * n_chi(i + 1, j, k) - 0.1
     #D1 * n_chi(i + 2, j, k)) / hx * n_dxdX(i, j, k) * (d * exp(c * x(i
     #)) - 0.1D1 * d * exp(-0.1D1 * c * x(i))) + 0.8333333333333333D-1 *
     # (n_chi(i, j - 2, k) - 0.8D1 * n_chi(i, j - 1, k) + 0.8D1 * n_chi(
     #i, j + 1, k) - 0.1D1 * n_chi(i, j + 2, k)) / hy * n_dydY(i, j, k) 
     #* (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j))) + 0.833
     #3333333333333D-1 * (n_chi(i, j, k - 2) - 0.8D1 * n_chi(i, j, k - 1
     #) + 0.8D1 * n_chi(i, j, k + 1) - 0.1D1 * n_chi(i, j, k + 2)) / hz 
     #* n_dzdZ(i, j, k) * (d * exp(c * z(k)) - 0.1D1 * d * exp(-0.1D1 * 
     #c * z(k))) + n_chi(i, j, k)) * (0.1000000000000000D21 * (d * exp(c
     # * x(i)) - 0.1D1 * d * exp(-0.1D1 * c * x(i))) ** 2 + 0.1000000000
     #000000D21 * (d * exp(c * y(j)) - 0.1D1 * d * exp(-0.1D1 * c * y(j)
     #)) ** 2 + 0.1000000000000000D21 * (d * exp(c * z(k)) - 0.1D1 * d *
     # exp(-0.1D1 * c * z(k))) ** 2 + 0.1D1) ** (-0.1D1 / 0.2D1)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
