cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     etotcalc: manually integrates the total energy over the domain
c               using an O(h^2) trapezoidal rule
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function etotcalc(cmask,Eden_np1,Eden_n,g1_NR,g1_NX,R,X,
     &dR,dX,gw,base,c,d)

      implicit none

      include 'globals.inc'
      include 'cmask.inc'

      real, parameter :: pi = 3.14159265359d0

      integer i, j
      integer g1_NR,g1_NX
      real*8  cmask(g1_NR,g1_NX)
      real*8  Eden_np1(g1_NR,g1_NX)
      real*8  Eden_n(g1_NR,g1_NX)
      real*8  R(*), X(*)
      real*8  dR, dX
      integer gw(4)
      integer base
      real*8  Etotsum

      real*8  c, d
      real*8, allocatable :: jac(:,:)

      allocate(jac(g1_NR,g1_NX))

      Etotsum = 0.0d0

      do i=1, g1_NR
      do j=1, g1_NX

        ! linear coordinates
c        jac(i,j) = R(i)

        ! compactified coordinates
        jac(i,j) =
     #  c ** 2 * d ** 3 * (exp(c * X(j)) + exp(-c * X(j))) * (exp(c 
     #* R(i)) + exp(-c * R(i))) * (exp(c * R(i)) - exp(-c * R(i)))

      end do
      end do

      do i=1+gw(1), g1_NR-1-gw(2)
      do j=1+gw(3), g1_NX-1-gw(4)

        ! when integrating over base level only, ignore cmask
        if ( base .eq. 1 ) then

          Etotsum = Etotsum + dR*dX/4.0d0*
     &  (  jac(i,j)*Eden_np1(i,j)     + jac(i+1,j)*Eden_np1(i+1,j)
     &   + jac(i,j+1)*Eden_np1(i,j+1) + jac(i+1,j+1)*Eden_np1(i+1,j+1) )

        else

          if ( ( cmask(i,j)     .eq. CMASK_off ) .and.
     &         ( cmask(i,j+1)   .ne. CMASK_on )  .and.
     &         ( cmask(i+1,j)   .ne. CMASK_on )  .and.
     &         ( cmask(i+1,j+1) .ne. CMASK_on ) ) then

            Etotsum = Etotsum + dR*dX/4.0d0*
     &  (  jac(i,j)*Eden_np1(i,j)     + jac(i+1,j)*Eden_np1(i+1,j)
     &   + jac(i,j+1)*Eden_np1(i,j+1) + jac(i+1,j+1)*Eden_np1(i+1,j+1) )

          end if

        end if

      end do
      end do

      etotcalc = 2d0*pi*Etotsum

      return
      end function etotcalc
