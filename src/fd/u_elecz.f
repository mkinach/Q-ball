      subroutine u_elecz(np1_At,np1_Bz,np1_dzdZ,Nx,Ny,Nz,hz,phys_bdy,res
     &)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 hz
      real*8 np1_At(Nx,Ny,Nz)
      real*8 np1_Bz(Nx,Ny,Nz)
      real*8 np1_dzdZ(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.8333333333333333D-1 * (np1_At(i, j, k - 2) - 0.8D1 * np1_At
     #(i, j, k - 1) + 0.8D1 * np1_At(i, j, k + 1) - 0.1D1 * np1_At(i, j,
     # k + 2)) / hz * np1_dzdZ(i, j, k) - 0.1D1 * np1_Bz(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
