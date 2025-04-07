cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init_qball: manually initialize some grid functions for gauged
c                 Q-ball evolution. Note that RNPL supplies the routine
c                 header, declaration of some routine arguments, and
c                 RETURN/END statements
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer     i, j, k
      integer     NR, NX
      integer     ini_type
      integer     unum, stat

      real*8        pi
      parameter   ( pi = 3.14159265358979d0 )

      character(len=100) :: datafile = "../../run/input/initdata"

      real*8, allocatable :: radial_axis(:)
      real*8, allocatable :: Q_shape(:)     ! stores precomputed shape from file
      real*8, allocatable :: A_shape(:)
      real*8, allocatable :: phic(:,:)      ! stores Lorentz-contracted shape (for boosted case)
      real*8, allocatable :: Ac(:,:)
      real*8, allocatable :: dphic(:,:)     ! stores derivative of contracted shape (for boosted case)
      real*8, allocatable :: dAc(:,:)

      real*8, allocatable :: phi1_n1(:,:)
      real*8, allocatable :: phi2_n1(:,:)
      real*8, allocatable :: pi1_n1(:,:)
      real*8, allocatable :: pi2_n1(:,:)
      real*8, allocatable :: A_t_n1(:,:)
      real*8, allocatable :: A_z_n1(:,:)
      real*8, allocatable :: B_t_n1(:,:)
      real*8, allocatable :: B_z_n1(:,:)

      real*8, allocatable :: phi1_n2(:,:)
      real*8, allocatable :: phi2_n2(:,:)
      real*8, allocatable :: pi1_n2(:,:)
      real*8, allocatable :: pi2_n2(:,:)
      real*8, allocatable :: A_t_n2(:,:)
      real*8, allocatable :: A_z_n2(:,:)
      real*8, allocatable :: B_t_n2(:,:)
      real*8, allocatable :: B_z_n2(:,:)

      real*8 :: neville                     ! interpolation function in interp.f
      real*8, dimension(0:4) :: xx=0, yy=0  ! storage for passing to neville()
      real*8        rr, drr, Q, A
      integer       length       ! precomputed length of arrays along radial axis
      integer       rrindex      ! array index closest to a point along radial axis
      integer       qb           ! label for individual Q-balls
      real*8        gamma, v     ! Lorentz factor and velocity
      real*8        phase        ! phase factor used in case of two boosted Q-balls
      logical       anti         ! toggle for Q-ball/anti-Q-ball collision
      real*8        w2           ! storage for anti-Q-ball frequency
      real*8, allocatable :: Ac2(:,:)       ! storage for anti-Q-ball Ac data

      ! parameters for Gaussian initial data
      real*8, allocatable, dimension(:,:) :: u, uold, f
      real*8        cw           ! overrelaxation parameter
      integer       citer        ! current iteration
      integer       cmaxiter     ! max # of iterations
      real*8        cerr         ! l-infinity norm
      real*8        ctol         ! error tol (stopping criteria)
      real*8        res          ! residual

      NR = g1_NR
      NX = g1_NX

      ! user-defined parameters
      v     = 0.6d0
      gamma = 1d0/sqrt(1d0-v**2)
      phase = 0d0
      anti  = .FALSE.

      ! user-defined initial data type
      !  1 - single static Q-ball
      !  2 - single boosted Q-ball
      !  3 - two boosted Q-balls
      !  4 - Gaussian initial data
      ini_type = 1

      ! retrieve number of points in the data file
      length = 0
      unum = 95
      open(unit=unum, file=datafile, status='old',
     &     action='read')
      do
        read(unum,*,iostat=stat)
        if (stat .ne. 0) exit
        length = length + 1
      end do
      close(unit=unum, status='keep', iostat=stat)

      ! read initial data from file
      allocate(radial_axis(length))
      allocate(Q_shape(length))
      allocate(A_shape(length))
      stat = 0
      open(unit=unum, file=datafile, status='old',
     &     action='read', iostat=stat)
      do i=1,length
        read (unum,*) rr, Q, A
        radial_axis(i)=rr
        Q_shape(i)=Q
        A_shape(i)=A
        if (stat .ne. 0) exit
      end do
      close(unit=unum, status='keep', iostat=stat)

      ! assumes initial data is output on uniformly-spaced intervals
      drr = radial_axis(2)-radial_axis(1)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          SINGLE STATIC Q-BALL           cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (ini_type .eq. 1) then

      ! interpolate
      do i=1,NR
      do j=1,NX

        rr = sqrt((pR_n(i,j)-cR)**2d0+(zX_n(i,j)-cX)**2d0)

        ! retrieve array index closest to rr
        rrindex=idnint(rr/drr)+1
        if ( rrindex+2 .ge. length ) then
          STOP 'Interpolation error (length exceeded)'
        end if
        if ( rrindex .le. 2 ) then  ! fix for interpolation near 0
          rrindex = 3
        end if

        xx(0)=radial_axis(rrindex-2)
        xx(1)=radial_axis(rrindex-1)
        xx(2)=radial_axis(rrindex)
        xx(3)=radial_axis(rrindex+1)
        xx(4)=radial_axis(rrindex+2)

        yy(0)=Q_shape(rrindex-2)
        yy(1)=Q_shape(rrindex-1)
        yy(2)=Q_shape(rrindex)
        yy(3)=Q_shape(rrindex+1)
        yy(4)=Q_shape(rrindex+2)
        phi1_n(i,j)=neville(xx,yy,rr)

        yy(0)=A_shape(rrindex-2)
        yy(1)=A_shape(rrindex-1)
        yy(2)=A_shape(rrindex)
        yy(3)=A_shape(rrindex+1)
        yy(4)=A_shape(rrindex+2)
        A_t_n(i,j)=neville(xx,yy,rr)

      end do
      end do

      ! specify initial data
      phi2_n(:,:) = 0d0
      pi1_n(:,:)  = 0d0
      do i=1,NR
      do j=1,NX
        pi2_n(i,j) = w*phi1_n(i,j)
      end do
      end do
      Atilde_p_n(:,:)   = 0d0
      Atilde_phi_n(:,:) = 0d0
      A_z_n(:,:)        = 0d0
      B_t_n(:,:)        = 0d0
      Btilde_p_n(:,:)   = 0d0
      Btilde_phi_n(:,:) = 0d0
      B_z_n(:,:)        = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          SINGLE BOOSTED Q-BALL          cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 2 ) then

      allocate(phic(NR,NX))
      allocate(dphic(NR,NX))
      allocate(Ac(NR,NX))
      allocate(dAc(NR,NX))

      ! interpolate
      do i=1,NR
      do j=1,NX

        rr = sqrt((pR_n(i,j)-cR)**2d0+(gamma*(zX_n(i,j)-cX))**2d0)

        ! retrieve array index closest to rr
        rrindex=idnint(rr/drr)+1
        if ( rrindex+2 .ge. length ) then
          STOP 'Interpolation error (length exceeded)'
        end if
        if ( rrindex .le. 2 ) then  ! fix for interpolation near 0
          rrindex = 3
        end if

        xx(0)=radial_axis(rrindex-2)
        xx(1)=radial_axis(rrindex-1)
        xx(2)=radial_axis(rrindex)
        xx(3)=radial_axis(rrindex+1)
        xx(4)=radial_axis(rrindex+2)

        yy(0)=Q_shape(rrindex-2)
        yy(1)=Q_shape(rrindex-1)
        yy(2)=Q_shape(rrindex)
        yy(3)=Q_shape(rrindex+1)
        yy(4)=Q_shape(rrindex+2)
        phic(i,j)=neville(xx,yy,rr)

        yy(0)=A_shape(rrindex-2)
        yy(1)=A_shape(rrindex-1)
        yy(2)=A_shape(rrindex)
        yy(3)=A_shape(rrindex+1)
        yy(4)=A_shape(rrindex+2)
        Ac(i,j)=neville(xx,yy,rr)

      end do
      end do

      ! compute derivatives using O(h^2) differencing
      do i=1,NR

        j=1
        dphic(i,j) = dXdz_n(i,j)*(-3*phic(i,j) + 4*phic(i,j+1)
     &               - phic(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dphic(i,j) = dXdz_n(i,j)*(phic(i,j+1) - phic(i,j-1))
     &               / (2 * gamma * dX)
        end do

        j=NX
        dphic(i,j) = dXdz_n(i,j)*(3*phic(i,j) - 4*phic(i,j-1)
     &               + phic(i,j-2)) / (2 * gamma * dX)

      end do

      do i=1,NR

        j=1
        dAc(i,j) = dXdz_n(i,j)*(-3*Ac(i,j) + 4*Ac(i,j+1)
     &             - Ac(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dAc(i,j) = dXdz_n(i,j)*(Ac(i,j+1) - Ac(i,j-1))
     &             / (2 * gamma * dX)
        end do

        j=NX
        dAc(i,j) = dXdz_n(i,j)*(3*Ac(i,j) - 4*Ac(i,j-1)
     &             + Ac(i,j-2)) / (2 * gamma * dX)

      end do

      ! specify initial data
      do i=1,NR
      do j=1,NX
        phi1_n(i,j) = phic(i,j)*cos(w*gamma*v*(zX_n(i,j)-cX))
        phi2_n(i,j) = phic(i,j)*sin(w*gamma*v*(zX_n(i,j)-cX))
        pi1_n(i,j)  = dphic(i,j)*cos(w*gamma*v*
     &                (zX_n(i,j)-cX))*gamma*v - w*gamma*phi2_n(i,j)
        pi2_n(i,j)  = dphic(i,j)*sin(w*gamma*v*
     &                (zX_n(i,j)-cX))*gamma*v + w*gamma*phi1_n(i,j)
      end do
      end do

      do i=1,NR
      do j=1,NX
        A_t_n(i,j) =   gamma*Ac(i,j)
        A_z_n(i,j) = v*gamma*Ac(i,j)
        B_t_n(i,j) = v*(gamma**2.0d0)*dAc(i,j)
        B_z_n(i,j) = (v**2.0d0)*(gamma**2.0d0)*dAc(i,j)
      end do
      end do

      Atilde_p_n(:,:)   = 0d0
      Atilde_phi_n(:,:) = 0d0
      Btilde_p_n(:,:)   = 0d0
      Btilde_phi_n(:,:) = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc           TWO BOOSTED Q-BALLS           cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 3 ) then

      allocate(phic(NR,NX))
      allocate(dphic(NR,NX))
      allocate(Ac(NR,NX))
      allocate(Ac2(NR,NX))
      allocate(dAc(NR,NX))

      ! Q-ball #1 setup
      allocate(phi1_n1(NR,NX))
      allocate(phi2_n1(NR,NX))
      allocate(pi1_n1(NR,NX))
      allocate(pi2_n1(NR,NX))
      allocate(A_t_n1(NR,NX))
      allocate(A_z_n1(NR,NX))
      allocate(B_t_n1(NR,NX))
      allocate(B_z_n1(NR,NX))

      cR = cR1
      cX = cX1

      ! interpolate
      do i=1,NR
      do j=1,NX

        rr = sqrt((pR_n(i,j)-cR)**2d0+(gamma*(zX_n(i,j)-cX))**2d0)

        ! retrieve array index closest to rr
        rrindex=idnint(rr/drr)+1
        if ( rrindex+2 .ge. length ) then
          STOP 'Interpolation error (length exceeded)'
        end if
        if ( rrindex .le. 2 ) then  ! fix for interpolation near 0
          rrindex = 3
        end if

        xx(0)=radial_axis(rrindex-2)
        xx(1)=radial_axis(rrindex-1)
        xx(2)=radial_axis(rrindex)
        xx(3)=radial_axis(rrindex+1)
        xx(4)=radial_axis(rrindex+2)

        yy(0)=Q_shape(rrindex-2)
        yy(1)=Q_shape(rrindex-1)
        yy(2)=Q_shape(rrindex)
        yy(3)=Q_shape(rrindex+1)
        yy(4)=Q_shape(rrindex+2)
        phic(i,j)=neville(xx,yy,rr)

        yy(0)=A_shape(rrindex-2)
        yy(1)=A_shape(rrindex-1)
        yy(2)=A_shape(rrindex)
        yy(3)=A_shape(rrindex+1)
        yy(4)=A_shape(rrindex+2)
        Ac(i,j)=neville(xx,yy,rr)

      end do
      end do

      ! compute derivatives using O(h^2) differencing
      do i=1,NR

        j=1
        dphic(i,j) = dXdz_n(i,j)*(-3*phic(i,j) + 4*phic(i,j+1)
     &               - phic(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dphic(i,j) = dXdz_n(i,j)*(phic(i,j+1) - phic(i,j-1))
     &               / (2 * gamma * dX)
        end do

        j=NX
        dphic(i,j) = dXdz_n(i,j)*(3*phic(i,j) - 4*phic(i,j-1)
     &               + phic(i,j-2)) / (2 * gamma * dX)

      end do

      do i=1,NR

        j=1
        dAc(i,j) = dXdz_n(i,j)*(-3*Ac(i,j) + 4*Ac(i,j+1)
     &             - Ac(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dAc(i,j) = dXdz_n(i,j)*(Ac(i,j+1) - Ac(i,j-1))
     &             / (2 * gamma * dX)
        end do

        j=NX
        dAc(i,j) = dXdz_n(i,j)*(3*Ac(i,j) - 4*Ac(i,j-1)
     &             + Ac(i,j-2)) / (2 * gamma * dX)

      end do

      ! specify initial data
      do i=1,NR
      do j=1,NX
        phi1_n1(i,j) = phic(i,j)*cos(w*gamma*v*(zX_n(i,j)-cX))
        phi2_n1(i,j) = phic(i,j)*sin(w*gamma*v*(zX_n(i,j)-cX))
        pi1_n1(i,j)  = dphic(i,j)*cos(w*gamma*v*
     &                 (zX_n(i,j)-cX))*gamma*v - w*gamma*phi2_n1(i,j)
        pi2_n1(i,j)  = dphic(i,j)*sin(w*gamma*v*
     &                 (zX_n(i,j)-cX))*gamma*v + w*gamma*phi1_n1(i,j)
      end do
      end do

      do i=1,NR
      do j=1,NX
        A_t_n1(i,j) =   gamma*Ac(i,j)
        A_z_n1(i,j) = v*gamma*Ac(i,j)
        B_t_n1(i,j) = v*(gamma**2.0d0)*dAc(i,j)
        B_z_n1(i,j) = (v**2.0d0)*(gamma**2.0d0)*dAc(i,j)
      end do
      end do

      ! Q-ball #2 setup
      allocate(phi1_n2(NR,NX))
      allocate(phi2_n2(NR,NX))
      allocate(pi1_n2(NR,NX))
      allocate(pi2_n2(NR,NX))
      allocate(A_t_n2(NR,NX))
      allocate(A_z_n2(NR,NX))
      allocate(B_t_n2(NR,NX))
      allocate(B_z_n2(NR,NX))

      cR = cR2
      cX = cX2

      v  = -v

      ! interpolate
      do i=1,NR
      do j=1,NX

        rr = sqrt((pR_n(i,j)-cR)**2d0+(gamma*(zX_n(i,j)-cX))**2d0)

        ! retrieve array index closest to rr
        rrindex=idnint(rr/drr)+1
        if ( rrindex+2 .ge. length ) then
          STOP 'Interpolation error (length exceeded)'
        end if
        if ( rrindex .le. 2 ) then  ! fix for interpolation near 0
          rrindex = 3
        end if

        xx(0)=radial_axis(rrindex-2)
        xx(1)=radial_axis(rrindex-1)
        xx(2)=radial_axis(rrindex)
        xx(3)=radial_axis(rrindex+1)
        xx(4)=radial_axis(rrindex+2)

        yy(0)=Q_shape(rrindex-2)
        yy(1)=Q_shape(rrindex-1)
        yy(2)=Q_shape(rrindex)
        yy(3)=Q_shape(rrindex+1)
        yy(4)=Q_shape(rrindex+2)
        phic(i,j)=neville(xx,yy,rr)

        yy(0)=A_shape(rrindex-2)
        yy(1)=A_shape(rrindex-1)
        yy(2)=A_shape(rrindex)
        yy(3)=A_shape(rrindex+1)
        yy(4)=A_shape(rrindex+2)
        Ac(i,j)=neville(xx,yy,rr)

      end do
      end do

      ! construct anti-Q-ball
      if (anti .eqv. .TRUE.) then
        w2 = -1d0*w
        do i=1,NR
        do j=1,NX
          Ac2(i,j)=-1d0*Ac(i,j)
        end do
        end do
      else
        w2 = w
        do i=1,NR
        do j=1,NX
          Ac2(i,j)=Ac(i,j)
        end do
        end do
      end if

      ! compute derivatives using O(h^2) differencing
      do i=1,NR

        j=1
        dphic(i,j) = dXdz_n(i,j)*(-3*phic(i,j) + 4*phic(i,j+1)
     &               - phic(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dphic(i,j) = dXdz_n(i,j)*(phic(i,j+1) - phic(i,j-1))
     &               / (2 * gamma * dX)
        end do

        j=NX
        dphic(i,j) = dXdz_n(i,j)*(3*phic(i,j) - 4*phic(i,j-1)
     &               + phic(i,j-2)) / (2 * gamma * dX)

      end do

      do i=1,NR

        j=1
        dAc(i,j) = dXdz_n(i,j)*(-3*Ac2(i,j) + 4*Ac2(i,j+1)
     &             - Ac2(i,j+2)) / (2 * gamma * dX)

        do j=2,NX-1
        dAc(i,j) = dXdz_n(i,j)*(Ac2(i,j+1) - Ac2(i,j-1))
     &             / (2 * gamma * dX)
        end do

        j=NX
        dAc(i,j) = dXdz_n(i,j)*(3*Ac2(i,j) - 4*Ac2(i,j-1)
     &             + Ac2(i,j-2)) / (2 * gamma * dX)

      end do

      ! specify initial data
      do i=1,NR
      do j=1,NX
        phi1_n2(i,j) = phic(i,j)*cos(w2*gamma*v*(zX_n(i,j)-cX)+phase)
        phi2_n2(i,j) = phic(i,j)*sin(w2*gamma*v*(zX_n(i,j)-cX)+phase)
        pi1_n2(i,j)  = dphic(i,j)*cos(w2*gamma*v*(zX_n(i,j)-cX)+phase)
     &                                *gamma*v - w2*gamma*phi2_n2(i,j)
        pi2_n2(i,j)  = dphic(i,j)*sin(w2*gamma*v*(zX_n(i,j)-cX)+phase)
     &                                *gamma*v + w2*gamma*phi1_n2(i,j)
      end do
      end do

      do i=1,NR
      do j=1,NX
         A_t_n2(i,j) =   gamma*Ac2(i,j)
         A_z_n2(i,j) = v*gamma*Ac2(i,j)
         B_t_n2(i,j) = v*(gamma**2.0d0)*dAc(i,j)
         B_z_n2(i,j) = (v**2.0d0)*(gamma**2.0d0)*dAc(i,j)
      end do
      end do

      ! combine Q-ball #1 and Q-ball #2 data
      do i=1,NR
      do j=1,NX
        phi1_n(i,j) = phi1_n1(i,j) + phi1_n2(i,j)
        phi2_n(i,j) = phi2_n1(i,j) + phi2_n2(i,j)
        pi1_n(i,j)  = pi1_n1(i,j)  + pi1_n2(i,j)
        pi2_n(i,j)  = pi2_n1(i,j)  + pi2_n2(i,j)
        A_t_n(i,j)  = A_t_n1(i,j)  + A_t_n2(i,j)
        A_z_n(i,j)  = A_z_n1(i,j)  + A_z_n2(i,j)
        B_t_n(i,j)  = B_t_n1(i,j)  + B_t_n2(i,j)
        B_z_n(i,j)  = B_z_n1(i,j)  + B_z_n2(i,j)
      end do
      end do

      Atilde_p_n(:,:)   = 0d0
      Atilde_phi_n(:,:) = 0d0
      Btilde_p_n(:,:)   = 0d0
      Btilde_phi_n(:,:) = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          GAUSSIAN INITIAL DATA          cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 4 ) then

      allocate( u(NR,NX) )
      allocate( uold(NR,NX) )
      allocate( f(NR,NX) )

      ! relaxation parameters
      cmaxiter=5000000
      cw=0.275d0
c      ctol=1.0d-10
c      ctol=1.0d-11
      ctol=5.0d-13

      write(*,*)
      write(*,101) cmaxiter
      write(*,102) ctol
      write(*,103) cw
      write(*,*)
 101  format( 'Maximum iterations:        ', I13 )
 102  format( 'Specified tolerance:       ', E15.8 )
 103  format( 'Overrelaxation parameter:  ', E15.8 )

      ! set arbitrary field values
      do i=1,NR
      do j=1,NX

        ! off-axis
        pi1_n(i,j)      =  0d0
        pi2_n(i,j)      =  0d0
        Atilde_p_n(i,j) =  0.3d0*exp(-1.0d0*( ( pR_n(i,j)-13d0)**2.0d0
     &                     + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (6.0d0))
        A_z_n(i,j)      =  0.1d0*exp(-1.0d0*( ( pR_n(i,j)-13d0)**2.0d0
     &                     + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (6.5d0))

c        ! on-axis
c        pi1_n(i,j)      =  0.0d0
c        pi2_n(i,j)      =  0.0d0
c        Atilde_p_n(i,j) =  0.3d0*exp(-1.0d0*( ( pR_n(i,j)-0d0)**2.0d0
c     &                     + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (12.0d0))
c        A_z_n(i,j)      =  0.1d0*exp(-1.0d0*( ( pR_n(i,j)-0d0)**2.0d0
c     &                     + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (19.0d0))

        Atilde_phi_n(i,j) =  0d0
        Btilde_phi_n(i,j) =  0d0

      end do
      end do

      ! set up the initial guess and source term
      do i=1,NR
      do j=1,NX

        ! off-axis
        phi1_n(i,j) = 0.25d0*exp(-1.0d0*( ( pR_n(i,j)-12.5d0)**2.0d0
     &                + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (5.5d0))
        phi2_n(i,j) = 0.31d0*exp(-1.0d0*( ( pR_n(i,j)-13d0)**2.0d0
     &                + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (6.0d0))
        A_t_n(i,j)  = 0.4d0*exp(-1.0d0*( ( pR_n(i,j)-12.5d0)**2.0d0
     &                + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (5.0d0))

c        ! on-axis
c        phi1_n(i,j) = 0.25d0*exp(-1.0d0*( ( pR_n(i,j)-0d0)**2.0d0
c     &                + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (4.5d0))
c        phi2_n(i,j) = 0.31d0*exp(-1.0d0*( ( pR_n(i,j)-0d0)**2.0d0
c     &                + ( zR_n(i,j)-1.0d0)**2.0d0 ) / (5.5d0))
c        A_t_n(i,j)  = 0.4d0*exp(-1.0d0*( ( pR_n(i,j)-0d0)**2.0d0
c     &                + ( zX_n(i,j)-1.0d0)**2.0d0 ) / (5.0d0))

        f(i,j) = 2.0d0*e*
     &  ( phi2_n(i,j)*pi1_n(i,j) - phi1_n(i,j)*pi2_n(i,j) )
     &  + 2.0d0*e**2.0d0*A_t_n(i,j)*
     &  ( phi1_n(i,j)**2.0d0 + phi2_n(i,j)**2.0d0 )

        ! initial guess
        u(i,j)=f(i,j)
        uold(i,j)=f(i,j)

      end do
      end do

      ! solve the elliptic constraint using SOR
      do citer=1,cmaxiter
      uold=u

      ! boundary condition at R=0 (i.e. i=1)
      i=1
      do j=2,NX-1

      ! quadratic fit
      u(i,j) = -u(i + 2, j) / 0.3D1 + 0.4D1 / 0.3D1 * u(i + 1, j)
      end do

      ! boundary condition at R=Rmax (i.e. i=NR)
      i=NR
      do j=1,NX

      u(i,j)=-dRdp_n(i, j)* (pR_n(i, j) ** 2 * u(i - 2, j) - 4 * pR_n(i,
     # j) ** 2 * u(i - 1, j) + zX_n(i, j) ** 2 * u(i - 2, j) - 4 * zX_n(
     #i, j) ** 2 * u(i - 1, j)) / (3 * pR_n(i, j) ** 2 * dRdp_n(i, j) +
     #3 * zX_n(i, j) ** 2 * dRdp_n(i, j) + 2 * pR_n(i, j) * dR)

      end do

      ! interior points
      do i=2,NR-1
      do j=2,NX-1

      u(i,j)=u(i, j) + cw * (-((dRdp2_n(i, j) * dR - 2 * dRdp_n(i, j) **
     #2) * pR_n(i, j) + dR * dRdp_n(i, j)) * dX ** 2 * u(i - 1, j) + ((d
     #Rdp2_n(i, j) * dR + 2 * dRdp_n(i, j) ** 2) * pR_n(i, j) + dR * dRd
     #p_n(i, j)) * dX ** 2 * u(i + 1, j) - (dR ** 2 * (dXdz2_n(i, j) * d
     #X - 2 * dXdz_n(i, j) ** 2) * u(i, j - 1) - dR ** 2 * (dXdz2_n(i, j
     #) * dX + 2 * dXdz_n(i, j) ** 2) * u(i, j + 1) + 4 * dX ** 2 * u(i,
     # j) * dRdp_n(i, j) ** 2 + 2 * dR ** 2 * (dX ** 2 * f(i, j) + 2 * u
     #(i, j) * dXdz_n(i, j) ** 2)) * pR_n(i, j)) / pR_n(i, j) / (dR ** 2
     # * dXdz_n(i, j) ** 2 + dX ** 2 * dRdp_n(i, j) ** 2)

      end do
      end do

      ! boundary condition at X=Xmin (i.e. j=1)
      j=1
      do i=1,NR

      u(i,j)=-dXdz_n(i,j) * (u(i, j + 2) * pR_n(i, j) ** 2 + u(i, j + 2)
     # * zX_n(i, j) ** 2 - 4 * u(i, j + 1) * pR_n(i, j) ** 2 - 4 * u(i, 
     #j + 1) * zX_n(i, j) ** 2) / (3 * dXdz_n(i, j) * pR_n(i, j) ** 2 + 
     #3 * dXdz_n(i, j) * zX_n(i, j) ** 2 - 2 * zX_n(i, j) * dX)

      end do

      ! boundary condition at X=Xmax (i.e. j=NX)
      j=NX
      do i=1,NR

      u(i,j)=-dXdz_n(i,j)* (u(i, j - 2) * pR_n(i, j) ** 2 + u(i, j - 2)
     # * zX_n(i, j) ** 2 - 4 * u(i, j - 1) * pR_n(i, j) ** 2 - 4 * u(i, 
     #j - 1) * zX_n(i, j) ** 2) / (3 * dXdz_n(i, j) * pR_n(i, j) ** 2 + 
     #3 * dXdz_n(i, j) * zX_n(i, j) ** 2 + 2 * zX_n(i, j) * dX)

      end do

      ! stop when L2-norm of error is less than tolerance
      cerr = 0.0d0
      do i=1,NR
      do j=1,NX
        cerr = cerr + (uold(i,j)-u(i,j))**2.0d0
      end do
      end do
      if ( sqrt(cerr)/(NR*NX) < ctol ) exit
      if ( mod(citer,1000) .EQ. 0 ) write(*,*) citer, sqrt(cerr)/(NR*NX)

      end do  ! end main iterative do loop

      ! compute relevant quantities using O(h^2) differencing
      do i=2,NR-1
      do j=2,NX-1

        B_t_n(i,j) =
     # dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (Atilde_p_n(i + 1,
     # j) - Atilde_p_n(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 + (A_z_n(i
     #, j + 1) - A_z_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1

        Btilde_p_n(i,j) =
     # ((A_t_n(i + 1, j) - A_t_n(i - 1, j)) / dR * dRdp_n(i, j) / 0
     #.2D1 - (u(i + 1, j) - u(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1) / 
     #pR_n(i, j)

        B_z_n(i,j) =
     # (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (u(i, j + 1) - u(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1

      end do
      end do

      ! [1:1][2:NX-1]
      i=1
      do j=2,NX-1

      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (-0.3D1 / 0.2D1 *
     #dble(Atilde_p_n(i, j)) + dble(2 * Atilde_p_n(i + 1, j)) - Atilde_p
     #_n(i + 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (A_z_n(i, j + 1) - A_z
     #_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
      B_t_n(i,j)=res

      res = dble((2 * A_t_n(i, j) - 5 * A_t_n(i + 1, j) + 4 * A_t_n(i + 
     #2, j) - A_t_n(i + 3, j)) / dR ** 2 * dRdp_n(i, j) ** 2) + (-0.3D1 
     #/ 0.2D1 * dble(A_t_n(i, j)) + dble(2 * A_t_n(i + 1, j)) - dble(A_t
     #_n(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2 * u(i,
     # j) - 5 * u(i + 1, j) + 4 * u(i + 2, j) - u(i + 3, j)) / dR ** 2 *
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(u(i, j)) + dble(2 * u
     #(i + 1, j)) - dble(u(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j
     #)
      Btilde_p_n(i,j)=res

      res = (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (u(i, j + 1) - u(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
      B_z_n(i,j)=res

      end do

      ! [2:NR-1][NX:NX]
      j=NX
      do i=2,NR-1

      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (Atilde_p_n(i + 1,
     # j) - Atilde_p_n(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 + (0.3D1 /
     # 0.2D1 * A_z_n(i, j) - dble(2 * A_z_n(i, j - 1)) + A_z_n(i, j - 2)
     # / 0.2D1) / dX * dXdz_n(i, j)
      B_t_n(i,j)=res

      res = ((A_t_n(i + 1, j) - A_t_n(i - 1, j)) / dR * dRdp_n(i, j) / 0
     #.2D1 - (u(i + 1, j) - u(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1) / 
     #pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1)) + A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1 / 0.2D1 * u(i
     #, j) - dble(2 * u(i, j - 1)) + u(i, j - 2) / 0.2D1) / dX * dXdz_n(
     #i, j)
      B_z_n(i,j)=res

      end do

      ! [NR:NR][2:NX-1]
      i=NR
      do j=2,NX-1

      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (0.3D1 / 0.2D1 * d
     #ble(Atilde_p_n(i, j)) - dble(2 * Atilde_p_n(i - 1, j)) + Atilde_p_
     #n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (A_z_n(i, j + 1) - A_z_
     #n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
      B_t_n(i,j)=res

      res = ((0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i - 1, j)) + 
     #A_t_n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) - (0.3D1 / 0.2D1 * u(
     #i, j) - dble(2 * u(i - 1, j)) + u(i - 2, j) / 0.2D1) / dR * dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (u(i, j + 1) - u(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
      B_z_n(i,j)=res

      end do

      ! [2:NR-1][1:1]
      j=1
      do i=2,NR-1

      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (Atilde_p_n(i + 1,
     # j) - Atilde_p_n(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 + (-0.3D1 
     #/ 0.2D1 * A_z_n(i, j) + dble(2 * A_z_n(i, j + 1)) - A_z_n(i, j + 2
     #) / 0.2D1) / dX * dXdz_n(i, j)
      B_t_n(i,j)=res

      res = ((A_t_n(i + 1, j) - A_t_n(i - 1, j)) / dR * dRdp_n(i, j) / 0
     #.2D1 - (u(i + 1, j) - u(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1) / 
     #pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1)) - 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1 * u
     #(i, j) + dble(2 * u(i, j + 1)) - u(i, j + 2) / 0.2D1) / dX * dXdz_
     #n(i, j)
      B_z_n(i,j)=res

      end do

      ! [1:1][1:1]
      i=1
      j=1
      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (-0.3D1 / 0.2D1 *
     #dble(Atilde_p_n(i, j)) + dble(2 * Atilde_p_n(i + 1, j)) - Atilde_p
     #_n(i + 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (-0.3D1 / 0.2D1 * A_z_
     #n(i, j) + dble(2 * A_z_n(i, j + 1)) - A_z_n(i, j + 2) / 0.2D1) / d
     #X * dXdz_n(i, j)
      B_t_n(i,j)=res

      res = dble((2 * A_t_n(i, j) - 5 * A_t_n(i + 1, j) + 4 * A_t_n(i + 
     #2, j) - A_t_n(i + 3, j)) / dR ** 2 * dRdp_n(i, j) ** 2) + (-0.3D1 
     #/ 0.2D1 * dble(A_t_n(i, j)) + dble(2 * A_t_n(i + 1, j)) - dble(A_t
     #_n(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2 * u(i,
     # j) - 5 * u(i + 1, j) + 4 * u(i + 2, j) - u(i + 3, j)) / dR ** 2 *
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(u(i, j)) + dble(2 * u
     #(i + 1, j)) - dble(u(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j
     #)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1)) - 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1 * u
     #(i, j) + dble(2 * u(i, j + 1)) - u(i, j + 2) / 0.2D1) / dX * dXdz_
     #n(i, j)
      B_z_n(i,j)=res

      ! [1:1][NX:NX]
      i=1
      j=NX
      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (-0.3D1 / 0.2D1 *
     #dble(Atilde_p_n(i, j)) + dble(2 * Atilde_p_n(i + 1, j)) - Atilde_p
     #_n(i + 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (0.3D1 / 0.2D1 * A_z_n
     #(i, j) - dble(2 * A_z_n(i, j - 1)) + A_z_n(i, j - 2) / 0.2D1) / dX
     # * dXdz_n(i, j)
      B_t_n(i,j)=res

      res = dble((2 * A_t_n(i, j) - 5 * A_t_n(i + 1, j) + 4 * A_t_n(i +
     #2, j) - A_t_n(i + 3, j)) / dR ** 2 * dRdp_n(i, j) ** 2) + (-0.3D1
     #/ 0.2D1 * dble(A_t_n(i, j)) + dble(2 * A_t_n(i + 1, j)) - dble(A_t
     #_n(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2 * u(i,
     # j) - 5 * u(i + 1, j) + 4 * u(i + 2, j) - u(i + 3, j)) / dR ** 2 *
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(u(i, j)) + dble(2 * u
     #(i + 1, j)) - dble(u(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j
     #)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1)) + A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1 / 0.2D1 * u(i
     #, j) - dble(2 * u(i, j - 1)) + u(i, j - 2) / 0.2D1) / dX * dXdz_n(
     #i, j)
      B_z_n(i,j)=res

      ! [NR:NR][NX:NX]
      i=NR
      j=NX
      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (0.3D1 / 0.2D1 * d
     #ble(Atilde_p_n(i, j)) - dble(2 * Atilde_p_n(i - 1, j)) + Atilde_p_
     #n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (0.3D1 / 0.2D1 * A_z_n(
     #i, j) - dble(2 * A_z_n(i, j - 1)) + A_z_n(i, j - 2) / 0.2D1) / dX 
     #* dXdz_n(i, j)
      B_t_n(i,j)=res

      res = ((0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i - 1, j)) + 
     #A_t_n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) - (0.3D1 / 0.2D1 * u(
     #i, j) - dble(2 * u(i - 1, j)) + u(i - 2, j) / 0.2D1) / dR * dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1)) + A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1 / 0.2D1 * u(i
     #, j) - dble(2 * u(i, j - 1)) + u(i, j - 2) / 0.2D1) / dX * dXdz_n(
     #i, j)
      B_z_n(i,j)=res

      ! [NR:NR][1:1]
      i=NR
      j=1
      res = dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (0.3D1 / 0.2D1 * d
     #ble(Atilde_p_n(i, j)) - dble(2 * Atilde_p_n(i - 1, j)) + Atilde_p_
     #n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) + (-0.3D1 / 0.2D1 * A_z_n
     #(i, j) + dble(2 * A_z_n(i, j + 1)) - A_z_n(i, j + 2) / 0.2D1) / dX
     # * dXdz_n(i, j)
      B_t_n(i,j)=res

      res = ((0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i - 1, j)) + 
     #A_t_n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) - (0.3D1 / 0.2D1 * u(
     #i, j) - dble(2 * u(i - 1, j)) + u(i - 2, j) / 0.2D1) / dR * dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1)) - 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1 * u
     #(i, j) + dble(2 * u(i, j + 1)) - u(i, j + 2) / 0.2D1) / dX * dXdz_
     #n(i, j)
      B_z_n(i,j)=res

      end if
