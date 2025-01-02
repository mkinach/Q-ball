      subroutine u_qden(np1_At,np1_modphi,np1_phi1,np1_phi2,np1_pi1,np1_
     &pi2,Nx,Ny,Nz,e,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 e
      real*8 np1_At(Nx,Ny,Nz)
      real*8 np1_modphi(Nx,Ny,Nz)
      real*8 np1_phi1(Nx,Ny,Nz)
      real*8 np1_phi2(Nx,Ny,Nz)
      real*8 np1_pi1(Nx,Ny,Nz)
      real*8 np1_pi2(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = -0.2D1 * np1_phi1(i, j, k) * np1_pi2(i, j, k) + 0.2D1 * np1_p
     #hi2(i, j, k) * np1_pi1(i, j, k) + 0.2D1 * e * np1_At(i, j, k) * np
     #1_modphi(i, j, k) ** 2
      res(i,j,k)=qb
      end do
      end do
      end do
      END
