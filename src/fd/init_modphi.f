      subroutine init_modphi(n_phi1,n_phi2,Nx,Ny,Nz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 n_phi1(Nx,Ny,Nz)
      real*8 n_phi2(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = sqrt(n_phi1(i, j, k) ** 2 + n_phi2(i, j, k) ** 2)
      res(i,j,k)=qb
      end do
      end do
      end do
      END
