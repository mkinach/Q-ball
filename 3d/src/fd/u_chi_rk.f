      subroutine u_chi_rk(n_xi,np1_chi,x,y,z,Nx,Ny,Nz,myzero,phys_bdy,re
     &s)
      implicit none
      integer i
      integer j
      integer k
      integer Nx
      integer Ny
      integer Nz
      real*8 myzero
      real*8 n_xi(Nx,Ny,Nz)
      real*8 np1_chi(Nx,Ny,Nz)
      real*8 x(Nx)
      real*8 y(Ny)
      real*8 z(Nz)
      integer phys_bdy(6)
      real*8 res(Nx,Ny,Nz)
      real*8 qb
      do i=3, Nx-2, 1
      do j=3, Ny-2, 1
      do k=3, Nz-2, 1
      qb = n_xi(i, j, k)
      res(i,j,k)=qb
      end do
      end do
      end do
      if (phys_bdy(1) .eq. 1) then
      do i=1, 1, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(2) .eq. 1) then
      do i=Nx, Nx, 1
      do j=1, Ny, 1
      do k=1, Nz, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(3) .eq. 1) then
      do i=1, Nx, 1
      do j=1, 1, 1
      do k=1, Nz, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(4) .eq. 1) then
      do i=1, Nx, 1
      do j=Ny, Ny, 1
      do k=1, Nz, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(5) .eq. 1) then
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=1, 1, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      if (phys_bdy(6) .eq. 1) then
      do i=1, Nx, 1
      do j=1, Ny, 1
      do k=Nz, Nz, 1
      qb = np1_chi(i, j, k) - 0.1D1 * myzero * x(i) * y(j) * z(k)
      res(i,j,k)=qb
      end do
      end do
      end do
      endif
      END
