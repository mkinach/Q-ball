c=======================================================================
c     Code to manually initialize some grid functions. Note that RNPL
c     supplies the routine header, declaration of some routine
c     arguments, and RETURN and END statements.
c=======================================================================

      integer     Nx, Ny      ! a more convenient spacing variable
      integer     i, j
      integer     ini_type
      integer     whichqball  ! for telling apart multiple Q-balls

c-----------------------------------------------------------------------
c     Some declarations for LSODA
c-----------------------------------------------------------------------
      
      integer        iargc
      real*8         r8arg
      real*8         rmax,       tol
      external       fcn,        jac 
      integer        neq 
      parameter    ( neq = 2 )         ! the equations are in fcn.f
      real*8         ysol(neq)
      real*8         r,          rout  ! axis variables
      integer        itol
      real*8         rtol,       atol
      integer        itask,      istate,     iopt
      integer        lrw 
      parameter    ( lrw = 22 + neq * 16 )
      real*8         rwork(lrw)
      integer        liw 
      parameter    ( liw = 20 + neq )
      integer        iwork(liw)
      integer        jt  
      include       'fcn.inc'

c-----------------------------------------------------------------------
c     Some declarations for the interpolation & boost
c-----------------------------------------------------------------------

      real*8, allocatable :: phi_shape(:)      ! stores interpolated phi 
      real*8, allocatable :: phi_shape_log(:)  ! stores linearized log of phi for exp tail fit 
      real*8, allocatable :: radial_axis(:)    ! stores radial axis values (0, 0.02, ...)
      real*8, allocatable :: x_p(:)            ! stores Lorentz contracted x-axis values
      real*8, allocatable :: phi_cntrct(:,:)   ! stores Lorentz contracted phi_shape
      real*8, allocatable :: d_phi_cntrct(:,:) ! stores derivative of phi_cntrct

      ! for binary boosted qball case; stores qball 1 or 2 data
      real*8, allocatable :: phi1_n1(:,:)
      real*8, allocatable :: phi2_n1(:,:)
      real*8, allocatable :: pi1_n1(:,:)
      real*8, allocatable :: pi2_n1(:,:)
      real*8, allocatable :: phi1_n2(:,:)
      real*8, allocatable :: phi2_n2(:,:)
      real*8, allocatable :: pi1_n2(:,:)
      real*8, allocatable :: pi2_n2(:,:)

      real*8         drout          ! LSODA integration step size
      real*8         radius(1)      ! radial distance to qball center
      real*8         phi_interp(1)  ! stores interpolated phi value at r=radius(1)
      integer        gridsize       ! Nx*Ny
      integer        length         ! precomputed length of arrays along radial axis
      real*8         dx_p           ! Lorentz contracted dx for boosted qball
      
      real*8         sigma0         ! initial qball height
      real*8         v              ! velocity of the Lorentz boost
      real*8         gamma          ! Lorentz factor calculated from v

      real*8         slope              ! slope for linear log approximation of phi
      real*8         exppt1, exppt2     ! y-points for slope
      real*8         exppt1_r, exppt2_r ! x-points for slope
      integer        expfit_flag        ! flag to collect x-points and y-points
       
c-----------------------------------------------------------------------
c     Model parameters
c-----------------------------------------------------------------------

      sigma0 = 1.375
      w      = 0.7739293

      B      = 0.5d0  ! defined for use in fcn.f via common block
      A      = 1d0      
     
      v      = 0.9    ! velocity of boosted Q-ball (c=1)
      gamma  = 1d0/sqrt(1d0-v**2) 

c-----------------------------------------------------------------------
c     LSODA parameters
c-----------------------------------------------------------------------

      tol    = 1.0d-10
      itol   = 1 
      rtol   = tol 
      atol   = tol 
      itask  = 1 
      iopt   = 0 
      jt     = 2
      
      ysol(1) = sigma0  ! initial condition for fcn.f
      ysol(2) = 0

      ! set max integration range. This is made large enough so
      ! the Q-ball is always interpolated over the entire domain
      if ( (xmax-xmin) .ge. (ymax-ymin) ) then
        rmax = sqrt(2d0)*(xmax-xmin)*gamma
      else
        rmax = sqrt(2d0)*(ymax-ymin)*gamma
      end if 
      
      r      = 0.0d0
      drout  = 0.02     ! integration step size   

c-----------------------------------------------------------------------
c     Some new variables for convenience
c-----------------------------------------------------------------------
 
      Nx       = g1_Nx  ! rename for convenience
      Ny       = g1_Ny  ! rename for convenience
      gridsize = Nx*Ny
      length   = int(rmax/drout)+2  ! length of arrays along radial axis

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Select initial data type:
c
c        1  Single static Q-ball
c        2  Single boosted Q-ball
c        3  Two identical boosted Q-balls
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ini_type = 2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
cccccccccccccc          SINGLE STATIC Q-BALL           ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      if ( ini_type .eq. 1 ) then

      allocate(phi_shape(length))
      allocate(phi_shape_log(length))
      allocate(radial_axis(length))

      !-------------------------
      !  Numerically integrate 
      !-------------------------

      istate  = 1
      i       = 1
      radial_axis(1) = r
      phi_shape(1)   = sigma0
      expfit_flag    = 0

      do while( r .lt. rmax )
         i      = i + 1 
         rout   = r + drout
         call lsoda(fcn,neq,ysol,r,rout,
     &              itol,rtol,atol,itask,
     &              istate,iopt,rwork,lrw,iwork,liw,jac,jt)

         if( istate .lt. 0 ) then
            write(0,*) 'shape:  Error return ', istate,
     &                 ' from LSODA '
            write(0,*) 'shape: Current interval ',
     &                 r, r + drout
            stop
         end if

         radial_axis(i) = r       ! axis values 0, 0.02, 0.04, ...
         phi_shape(i)   = ysol(1) ! phi function value at each r

         if ( r.gt.16 .and. expfit_flag.eq.0) then
            exppt1 = log(ysol(1))
            exppt1_r = r
            expfit_flag = 1
         endif
         if ( r.gt.17 .and. expfit_flag.eq.1) then
            exppt2 = log(ysol(1))
            exppt2_r = r
            expfit_flag = -1
         endif
        
      end do

      !----------------------------
      !  Prepare exponential tail
      !----------------------------

      ! calculate slope
      slope = (exppt1-exppt2)/(exppt1_r-exppt2_r) 

      j = 1 
      do while( j .lt. i ) ! i is length of non-zero elements of radial_axis
         phi_shape_log(j) = slope*(radial_axis(j)-exppt1_r)+exppt1
         j = j + 1
      enddo

      !-------------------------
      !  Interpolate
      !-------------------------

      do i=1,Nx   ! interpolate the static solution over the entire 2D plane
      do j=1,Ny

         radius(1) = sqrt((x(i)-cx)**2+(y(j)-cy)**2)
        
         if ( radius(1) .lt. 18 ) then
            call dvinqn(phi_shape,radial_axis,phi_interp(1),radius(1),
     &               length,1,phi_shape(1),phi_shape(length),4)
            phi1_n(i,j) = phi_interp(1)
         else
            call dvinqn(phi_shape_log,radial_axis,phi_interp(1),
     &        radius(1),length,1,phi_shape_log(1),
     &        phi_shape_log(length),4)
            phi1_n(i,j) = exp(phi_interp(1)) 
         endif
         
         pi2_n(i,j) = w*phi1_n(i,j);

      end do
      end do

      phi2_n(:,:) = 0d0
      pi1_n(:,:)  = 0d0


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccc          SINGLE BOOSTED Q-BALL          ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      elseif ( ini_type .eq. 2 ) then
     
      allocate(phi_shape(length))
      allocate(phi_shape_log(length))
      allocate(x_p(length))
      allocate(radial_axis(length))
      allocate(phi_cntrct(Nx,Ny))
      allocate(d_phi_cntrct(Nx,Ny))

      !-------------------------
      !  Numerically integrate 
      !-------------------------

      istate  = 1
      i       = 1
      radial_axis(1) = r
      phi_shape(1)   = sigma0
      expfit_flag    = 0

      do while( r .lt. rmax )
         i      = i + 1 
         rout   = r + drout
         call lsoda(fcn,neq,ysol,r,rout,
     &              itol,rtol,atol,itask,
     &              istate,iopt,rwork,lrw,iwork,liw,jac,jt)

         if( istate .lt. 0 ) then
            write(0,*) 'shape:  Error return ', istate,
     &                 ' from LSODA '
            write(0,*) 'shape: Current interval ',
     &                 r, r + drout
            stop
         end if

         radial_axis(i) = r      
         phi_shape(i)   = ysol(1)

         if ( r.gt.16 .and. expfit_flag.eq.0) then 
            exppt1 = log(ysol(1))
            exppt1_r = r
            expfit_flag = 1
         endif
         if ( r.gt.17 .and. expfit_flag.eq.1) then
            exppt2 = log(ysol(1))
            exppt2_r = r
            expfit_flag = -1
         endif

      end do

      !----------------------------
      !  Prepare exponential tail
      !----------------------------

      ! calculate slope
      slope = (exppt1-exppt2)/(exppt1_r-exppt2_r) 

      j = 1 
      do while( j .lt. i ) ! i is length of non-zero elements of radial_axis
         phi_shape_log(j) = slope*(radial_axis(j)-exppt1_r)+exppt1
         j = j + 1
      enddo

      !-------------------------
      !  Interpolate
      !-------------------------

      dx_p = gamma*dx   ! rescale dx according to Lorentz contraction

      do i=1,Nx   ! rescale the x-axis according to Lorentz contraction
         x_p(i) = gamma*x(i)
      enddo

      do i=1,Nx   ! interpolate the contracted solution over the entire 2D plane
      do j=1,Ny

         radius(1) = sqrt((x_p(i)-cx*gamma)**2+(y(j)-cy)**2)

         if ( radius(1) .lt. 18 ) then
            call dvinqn(phi_shape,radial_axis,phi_interp(1),radius(1),
     &                  length,1,phi_shape(1),phi_shape(length),4)
            phi_cntrct(i,j) = phi_interp(1)
         else
            call dvinqn(phi_shape_log,radial_axis,phi_interp(1),
     &        radius(1),length,1,phi_shape_log(1),
     &        phi_shape_log(length),4)
            phi_cntrct(i,j) = exp(phi_interp(1)) 
         endif

      end do
      end do

      !------------------------
      !  Calculate derivative
      !------------------------

      do j=1,Ny   ! uses Crank-Nicolson 

         i=1
         d_phi_cntrct(i,j) = (-3*phi_cntrct(i,j) + 4*phi_cntrct(i+1,j)
     &                        - phi_cntrct(i+2,j)) / (2 * dx_p); 
        
         do i=2,Nx-1  
         d_phi_cntrct(i,j) = (phi_cntrct(i+1,j) - phi_cntrct(i-1,j))
     &                        / (2 * dx_p);
         enddo

         i=Nx
         d_phi_cntrct(i,j) = (3*phi_cntrct(i,j) - 4*phi_cntrct(i-1,j)
     &                        + phi_cntrct(i-2,j)) / (2 * dx_p); 

      enddo

      !---------------------------
      !  Output boosted solution
      !---------------------------

      do i=1,Nx 
      do j=1,Ny

         phi1_n(i,j) = phi_cntrct(i,j)*cos(w*gamma*v*(x(i)-cx))
         phi2_n(i,j) = phi_cntrct(i,j)*sin(w*gamma*v*(x(i)-cx))
         pi1_n(i,j)  = d_phi_cntrct(i,j)*cos(w*gamma*v*
     &                 (x(i)-cx))*gamma*v - w*gamma*phi2_n(i,j)
         pi2_n(i,j)  = d_phi_cntrct(i,j)*sin(w*gamma*v*
     &                 (x(i)-cx))*gamma*v + w*gamma*phi1_n(i,j)

      enddo
      enddo


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccccccccc      TWO IDENTICAL BOOSTED Q-BALLS      ccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      

      elseif ( ini_type .eq. 3 ) then
     
      allocate(phi_shape(length))
      allocate(phi_shape_log(length))
      allocate(x_p(length))
      allocate(radial_axis(length))
      allocate(phi_cntrct(Nx,Ny))
      allocate(d_phi_cntrct(Nx,Ny))

      allocate(phi1_n1(Nx,Ny))
      allocate(phi2_n1(Nx,Ny))
      allocate(pi1_n1(Nx,Ny))
      allocate(pi2_n1(Nx,Ny))
      allocate(phi1_n2(Nx,Ny))
      allocate(phi2_n2(Nx,Ny))
      allocate(pi1_n2(Nx,Ny))
      allocate(pi2_n2(Nx,Ny))

      !-------------------------
      !  Numerically integrate 
      !-------------------------

      istate  = 1
      i       = 1
      radial_axis(1) = r
      phi_shape(1)   = sigma0
      expfit_flag    = 0

      do while( r .lt. rmax )
         i      = i + 1 
         rout   = r + drout
         call lsoda(fcn,neq,ysol,r,rout,
     &              itol,rtol,atol,itask,
     &              istate,iopt,rwork,lrw,iwork,liw,jac,jt)

         if( istate .lt. 0 ) then
            write(0,*) 'shape:  Error return ', istate,
     &                 ' from LSODA '
            write(0,*) 'shape: Current interval ',
     &                 r, r + drout
            stop
         end if

         radial_axis(i) = r       ! these 2 lines keep track of
         phi_shape(i)   = ysol(1) ! axis values and phi function values

         if ( r.gt.16 .and. expfit_flag.eq.0) then 
            exppt1 = log(ysol(1))
            exppt1_r = r
            expfit_flag = 1
         endif
         if ( r.gt.17 .and. expfit_flag .eq. 1) then
            exppt2 = log(ysol(1))
            exppt2_r = r
            expfit_flag = -1
         endif

      end do

      !----------------------------
      !  Prepare exponential tail
      !----------------------------

      ! calculate slope
      slope = (exppt1-exppt2)/(exppt1_r-exppt2_r) 

      j = 1 
      do while( j .lt. i ) ! i is length of non-zero elements of radial_axis
         phi_shape_log(j) = slope*(radial_axis(j)-exppt1_r)+exppt1
         j = j + 1
      enddo

      !-------------------------
      !  Interpolate
      !-------------------------

      dx_p = gamma*dx   ! rescale dx according to Lorentz contraction

      do i=1,Nx   ! rescale the x-axis according to Lorentz contraction
         x_p(i) = gamma*x(i)
      enddo

      do whichqball=1,2  ! superpose Q-balls at different positions
         if ( whichqball .eq. 1 ) then
            cx = cx1
            cy = cy1
         else
            cx = cx2
            cy = cy2
            v  = -v  ! sends the second Q-ball towards the first
         endif
      
      do i=1,Nx ! interpolate the contracted solution over the entire 2D plane
      do j=1,Ny

         radius(1) = sqrt((x_p(i)-cx*gamma)**2+(y(j)-cy)**2)

         if ( radius(1) .lt. 18 ) then
            call dvinqn(phi_shape,radial_axis,phi_interp(1),radius(1),
     &                  length,1,phi_shape(1),phi_shape(length),4)
            phi_cntrct(i,j) = phi_interp(1)
         else
            call dvinqn(phi_shape_log,radial_axis,phi_interp(1),
     &        radius(1),length,1,phi_shape_log(1),
     &        phi_shape_log(length),4)
            phi_cntrct(i,j) = exp(phi_interp(1)) 
         endif
         
      end do
      end do

      !------------------------
      !  Calculate derivative
      !------------------------

      do j=1,Ny  ! use Crank-Nicolson 

         i=1
         d_phi_cntrct(i,j) = (-3*phi_cntrct(i,j) + 4*phi_cntrct(i+1,j)
     &                        - phi_cntrct(i+2,j)) / (2 * dx_p); 
        
         do i=2,Nx-1  
         d_phi_cntrct(i,j) = (phi_cntrct(i+1,j) - phi_cntrct(i-1,j))
     &                        / (2 * dx_p);
         enddo

         i=Nx
         d_phi_cntrct(i,j) = (3*phi_cntrct(i,j) - 4*phi_cntrct(i-1,j)
     &                        + phi_cntrct(i-2,j)) / (2 * dx_p); 

      end do

      !---------------------------
      !  Output boosted solution
      !---------------------------

      if ( whichqball .eq. 1 ) then      
         do i=1,Nx 
         do j=1,Ny
   
            phi1_n1(i,j) = phi_cntrct(i,j)*
     &                     cos(w*gamma*v*(x(i)-cx))
            phi2_n1(i,j) = phi_cntrct(i,j)*
     &                     sin(w*gamma*v*(x(i)-cx))
            pi1_n1(i,j)  = d_phi_cntrct(i,j)*cos(w*gamma*v*
     &                     (x(i)-cx))*gamma*v - w*gamma*phi2_n1(i,j)
            pi2_n1(i,j)  = d_phi_cntrct(i,j)*sin(w*gamma*v*
     &                     (x(i)-cx))*gamma*v + w*gamma*phi1_n1(i,j)
         end do
         end do
      else
         do i=1,Nx 
         do j=1,Ny
   
            phi1_n2(i,j) = phi_cntrct(i,j)*
     &                     cos(w*gamma*v*(x(i)-cx))
            phi2_n2(i,j) = phi_cntrct(i,j)*
     &                     sin(w*gamma*v*(x(i)-cx))
            pi1_n2(i,j)  = d_phi_cntrct(i,j)*cos(w*gamma*v*
     &                     (x(i)-cx))*gamma*v - w*gamma*phi2_n2(i,j)
            pi2_n2(i,j)  = d_phi_cntrct(i,j)*sin(w*gamma*v*
     &                     (x(i)-cx))*gamma*v + w*gamma*phi1_n2(i,j)
         end do
         end do
 
      end if

      end do

      do i=1,Nx   ! combine the two Q-balls calculated above
      do j=1,Ny
         phi1_n(i,j) = phi1_n1(i,j) + phi1_n2(i,j)
         phi2_n(i,j) = phi2_n1(i,j) + phi2_n2(i,j)
         pi1_n(i,j)  = pi1_n1(i,j)  + pi1_n2(i,j)
         pi2_n(i,j)  = pi2_n1(i,j)  + pi2_n2(i,j)
      end do
      end do

      end if
