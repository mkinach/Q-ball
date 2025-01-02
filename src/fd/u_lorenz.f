      subroutine u_lorenz(np1_Ax,np1_Ay,np1_Az,np1_Bt,np1_dxdX,np1_dydY,
     &np1_dzdZ,Nx,Ny,Nz,hx,hy,hz,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hx
      real*8 hy
      real*8 hz
      real*8 np1_Ax(Nx,Ny,Nz)
      real*8 np1_Ay(Nx,Ny,Nz)
      real*8 np1_Az(Nx,Ny,Nz)
      real*8 np1_Bt(Nx,Ny,Nz)
      real*8 np1_dxdX(Nx,Ny,Nz)
      real*8 np1_dydY(Nx,Ny,Nz)
      real*8 np1_dzdZ(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = -0.1D1 * np1_Bt(i, j, k) + 0.8333333333333333D-1 * (np1_Ax(i 
     #- 2, j, k) - 0.8D1 * np1_Ax(i - 1, j, k) + 0.8D1 * np1_Ax(i + 1, j
     #, k) - 0.1D1 * np1_Ax(i + 2, j, k)) / hx * np1_dxdX(i, j, k) + 0.8
     #333333333333333D-1 * (np1_Ay(i, j - 2, k) - 0.8D1 * np1_Ay(i, j - 
     #1, k) + 0.8D1 * np1_Ay(i, j + 1, k) - 0.1D1 * np1_Ay(i, j + 2, k))
     # / hy * np1_dydY(i, j, k) + 0.8333333333333333D-1 * (np1_Az(i, j, 
     #k - 2) - 0.8D1 * np1_Az(i, j, k - 1) + 0.8D1 * np1_Az(i, j, k + 1)
     # - 0.1D1 * np1_Az(i, j, k + 2)) / hz * np1_dzdZ(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
