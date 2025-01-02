cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     w2f: write multigrid data to binary file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine w2f(V,Nx,Ny,Nz)

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz)

      integer i, j, k
      integer unum

      unum=20

      ! open the file
      open(unit=unum, file='V.bin', status='replace',
     &     form='unformatted', access='stream')
      
      ! write the array
      do k = 1, Nz
      do j = 1, Ny
      do i = 1, Nx
        write(unum) V(i, j, k)
      end do
      end do
      end do

      ! close the file
      close(unum)

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     rff: read multigrid data from binary file
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rff(V,Nx,Ny,Nz)

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz)

      integer i, j, k
      integer unum

      unum=20

      ! open the file
      open(unit=unum, file='V.bin', status='old',
     &     form='unformatted', access='stream')
       
      ! read the array
      do k = 1, Nz
      do j = 1, Ny
      do i = 1, Nx
        read(unum) V(i, j, k)
      end do
      end do
      end do

      ! close the file
      close(unum)

      return
      end subroutine
