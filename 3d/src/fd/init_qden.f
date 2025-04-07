      subroutine init_qden(n_At,n_modphi,n_phi1,n_phi2,n_pi1,n_pi2,Nx,Ny
     &,Nz,e,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 e
      real*8 n_At(Nx,Ny,Nz)
      real*8 n_modphi(Nx,Ny,Nz)
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 n_pi1(Nx,Ny,Nz)
      real*8 n_pi2(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = -0.2D1 * n_phi1(i, j, k) * n_pi2(i, j, k) + 0.2D1 * n_phi2(i,
     # j, k) * n_pi1(i, j, k) + 0.2D1 * e * n_At(i, j, k) * n_modphi(i, 
     #j, k) ** 2
      res(i,j,k)=qb
      end do
      end do
      end do
      END
