      subroutine init_eem(n_elecx,n_elecy,n_elecz,n_magx,n_magy,n_magz,N
     &x,Ny,Nz,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 n_elecx(Nx,Ny,Nz)
      real*8 n_elecy(Nx,Ny,Nz)
      real*8 n_elecz(Nx,Ny,Nz)
      real*8 n_magx(Nx,Ny,Nz)
      real*8 n_magy(Nx,Ny,Nz)
      real*8 n_magz(Nx,Ny,Nz)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = 0.5D0 * n_elecx(i, j, k) ** 2 + 0.5D0 * n_elecy(i, j, k) ** 2
     # + 0.5D0 * n_elecz(i, j, k) ** 2 + 0.5D0 * n_magx(i, j, k) ** 2 + 
     #0.5D0 * n_magy(i, j, k) ** 2 + 0.5D0 * n_magz(i, j, k) ** 2
      res(i,j,k)=qb
      end do
      end do
      end do
      END
