      subroutine u_eem(np1_elecx,np1_elecy,np1_elecz,np1_magx,np1_magy,n
     &p1_magz,Nx,Ny,Nz,phys_bdy,res)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 np1_elecx(Nx,Ny,Nz)
      real*8 np1_elecy(Nx,Ny,Nz)
      real*8 np1_elecz(Nx,Ny,Nz)
      real*8 np1_magx(Nx,Ny,Nz)
      real*8 np1_magy(Nx,Ny,Nz)
      real*8 np1_magz(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = 0.5D0 * np1_elecx(i, j, k) ** 2 + 0.5D0 * np1_elecy(i, j, k) 
     #** 2 + 0.5D0 * np1_elecz(i, j, k) ** 2 + 0.5D0 * np1_magx(i, j, k)
     # ** 2 + 0.5D0 * np1_magy(i, j, k) ** 2 + 0.5D0 * np1_magz(i, j, k)
     # ** 2
      res(i,j,k)=qb
      end do
      end do
      end do
      END
