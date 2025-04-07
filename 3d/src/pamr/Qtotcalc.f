cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     qtotcalc: manually integrates the total charge over the domain
c               using an O(h^2) trapezoidal rule
c
c               Note that this code is distinct from Etotcalc.f due to
c               added functionality for integration in the half-volume
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function qtotcalc(cmask,Qden_np1,Qden_n,jac_n,Nx,Ny,Nz,x,y,
     &z,dx,dy,dz,gw,base,half,c,d)

      implicit none

      integer CMASK_ON,CMASK_OFF
      parameter (CMASK_ON=1,CMASK_OFF=0)

      integer i, j, k
      integer Nx, Ny, Nz
      real*8  cmask(Nx,Ny,Nz)
      real*8  Qden_np1(Nx,Ny,Nz)
      real*8  Qden_n(Nx,Ny,Nz)
      real*8  jac_n(Nx,Ny,Nz)
      real*8  x(*), y(*), z(*)
      real*8  dx, dy, dz
      integer gw(6)
      integer base
      integer half
      real*8  c
      real*8  d
      real*8  Qtotsum

      Qtotsum = 0.0d0

      do i=1+gw(1), Nx-1-gw(2)
      do j=1+gw(3), Ny-1-gw(4)
      do k=1+gw(5), Nz-1-gw(6)

        ! for half-volume integration
        if ( ( half .eq.  1 ) .and. ( z(k) .le. 0d0 ) ) cycle
        if ( ( half .eq. -1 ) .and. ( z(k) .gt. 0d0 ) ) cycle

        ! when integrating over base level only, ignore cmask
        if ( base .eq. 1 ) then

          Qtotsum = Qtotsum + dx*dy*dz/8.0d0*
     &( jac_n(i,j,k)*abs(Qden_np1(i,j,k))
     &+ jac_n(i+1,j,k)*abs(Qden_np1(i+1,j,k))
     &+ jac_n(i,j+1,k)*abs(Qden_np1(i,j+1,k))
     &+ jac_n(i+1,j+1,k)*abs(Qden_np1(i+1,j+1,k))
     &+ jac_n(i,j,k+1)*abs(Qden_np1(i,j,k+1))
     &+ jac_n(i+1,j,k+1)*abs(Qden_np1(i+1,j,k+1))
     &+ jac_n(i,j+1,k+1)*abs(Qden_np1(i,j+1,k+1))
     &+ jac_n(i+1,j+1,k+1)*abs(Qden_np1(i+1,j+1,k+1)) )

c          Qtotsum = Qtotsum + dx*dy*dz/8.0d0*
c     &( jac_n(i,j,k)*Qden_np1(i,j,k)
c     &+ jac_n(i+1,j,k)*Qden_np1(i+1,j,k)
c     &+ jac_n(i,j+1,k)*Qden_np1(i,j+1,k)
c     &+ jac_n(i+1,j+1,k)*Qden_np1(i+1,j+1,k)
c     &+ jac_n(i,j,k+1)*Qden_np1(i,j,k+1)
c     &+ jac_n(i+1,j,k+1)*Qden_np1(i+1,j,k+1)
c     &+ jac_n(i,j+1,k+1)*Qden_np1(i,j+1,k+1)
c     &+ jac_n(i+1,j+1,k+1)*Qden_np1(i+1,j+1,k+1) )

        else

          if ( ( cmask(i,j,k)       .eq. CMASK_off )  .and.
     &         ( cmask(i,j+1,k)     .ne. CMASK_on )   .and.
     &         ( cmask(i+1,j,k)     .ne. CMASK_on )   .and.
     &         ( cmask(i+1,j+1,k)   .ne. CMASK_on )   .and.
     &         ( cmask(i,j,k+1)     .ne. CMASK_on )   .and.
     &         ( cmask(i+1,j,k+1)   .ne. CMASK_on )   .and.
     &         ( cmask(i,j+1,k+1)   .ne. CMASK_on )   .and.
     &         ( cmask(i+1,j+1,k+1) .ne. CMASK_on ) )  then

            Qtotsum = Qtotsum + dx*dy*dz/8.0d0*
     &( jac_n(i,j,k)*abs(Qden_np1(i,j,k))
     &+ jac_n(i+1,j,k)*abs(Qden_np1(i+1,j,k))
     &+ jac_n(i,j+1,k)*abs(Qden_np1(i,j+1,k))
     &+ jac_n(i+1,j+1,k)*abs(Qden_np1(i+1,j+1,k))
     &+ jac_n(i,j,k+1)*abs(Qden_np1(i,j,k+1))
     &+ jac_n(i+1,j,k+1)*abs(Qden_np1(i+1,j,k+1))
     &+ jac_n(i,j+1,k+1)*abs(Qden_np1(i,j+1,k+1))
     &+ jac_n(i+1,j+1,k+1)*abs(Qden_np1(i+1,j+1,k+1)) )

c            Qtotsum = Qtotsum + dx*dy*dz/8.0d0*
c     &( jac_n(i,j,k)*Qden_np1(i,j,k)
c     &+ jac_n(i+1,j,k)*Qden_np1(i+1,j,k)
c     &+ jac_n(i,j+1,k)*Qden_np1(i,j+1,k)
c     &+ jac_n(i+1,j+1,k)*Qden_np1(i+1,j+1,k)
c     &+ jac_n(i,j,k+1)*Qden_np1(i,j,k+1)
c     &+ jac_n(i+1,j,k+1)*Qden_np1(i+1,j,k+1)
c     &+ jac_n(i,j+1,k+1)*Qden_np1(i,j+1,k+1)
c     &+ jac_n(i+1,j+1,k+1)*Qden_np1(i+1,j+1,k+1) )

          end if

        end if

      end do
      end do
      end do

      qtotcalc = Qtotsum

      return
      end function qtotcalc
