cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     interp: interpolates a value at a point via the Neville algorithm.
c             Note that the interpolation is hard-coded to be
c             fourth-order. Input arguments are a list of 5 axis values,
c             a list of 5 known function values at those points, and a
c             point at which the interpolation is desired
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function neville(x, f, p)

      implicit none

      real*8, intent(in), dimension(0:4) :: x
      real*8, intent(in), dimension(0:4) :: f
      real*8, intent(in) :: p

      real*8, dimension(0:4,0:4) :: Q = 0d0
      integer i, j

      do i=0,4
        Q(i,0) = f(i)
      end do

      ! Neville algorithm to fourth-order
      do i=1,4
      do j=1,i
        Q(i,j)=(p-x(i-j))*Q(i,j-1)-(p-x(i))*Q(i-1,j-1)
        Q(i,j)=Q(i,j)/(x(i)-x(i-j))

        ! prevent underflow
        if ( Q(i,j) .lt. 1d-50 ) then
          neville = 0d0
        end if
      end do
      end do

c     write(*,*) x(0), Q(0,0)
c     write(*,*) x(1), Q(1,0), Q(1,1)
c     write(*,*) x(2), Q(2,0), Q(2,1), Q(2,2)
c     write(*,*) x(3), Q(3,0), Q(3,1), Q(3,2), Q(3,3)
c     write(*,*) x(4), Q(4,0), Q(4,1), Q(4,2), Q(4,3), Q(4,4)

      neville = Q(4,4)

      return
      end function neville
