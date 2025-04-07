      subroutine u_elecx(np1_At,np1_Bx,np1_dxdX,Nx,Ny,Nz,hx,phys_bdy,res
     &)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hx
      real*8 np1_At(Nx,Ny,Nz)
      real*8 np1_Bx(Nx,Ny,Nz)
      real*8 np1_dxdX(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (np1_At(i - 2, j, k) - 0.8D1 * np1_At
     #(i - 1, j, k) + 0.8D1 * np1_At(i + 1, j, k) - 0.1D1 * np1_At(i + 2
     #, j, k)) / hx * np1_dxdX(i, j, k) - 0.1D1 * np1_Bx(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
