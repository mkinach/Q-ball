      subroutine u_modphi(np1_phi1,np1_phi2,Nx,Ny,Nz,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 np1_phi1(Nx,Ny,Nz)
      real*8 np1_phi2(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = sqrt(np1_phi1(i, j, k) ** 2 + np1_phi2(i, j, k) ** 2)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
