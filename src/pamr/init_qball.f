cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init_qball: manually initialize some grid functions for Q-ball
c                 evolution
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_qball(Xx_n,dxdX_n,dxdX2_n,Yy_n,dydY_n,dydY2_n,
     &Zz_n,dzdZ_n,dzdZ2_n,phi1_np1,phi1_n,phi2_np1,phi2_n,pi1_np1,pi1_n,
     &pi2_np1,pi2_n,at_np1,at_n,ax_np1,ax_n,ay_np1,ay_n,az_np1,az_n,
     &bt_np1,bt_n,bx_np1,bx_n,by_np1,by_n,bz_np1,bz_n,Nx,Ny,Nz,x,y,z,dx,
     &dy,dz,xmax,xmin,ymax,ymin,zmax,zmin,cx,cy,cz,cx1,cx2,cy1,cy2,cz1,
     &cz2,w,v,phase,anti,ini_type,e)
      use interp_mod
      
      implicit none

      integer i, j, k
      integer Nx
      integer Ny
      integer Nz
      real*8  Xx_n(Nx,Ny,Nz)
      real*8  dxdX_n(Nx,Ny,Nz)
      real*8  dxdX2_n(Nx,Ny,Nz)
      real*8  Yy_n(Nx,Ny,Nz)
      real*8  dydY_n(Nx,Ny,Nz)
      real*8  dydY2_n(Nx,Ny,Nz)
      real*8  Zz_n(Nx,Ny,Nz)
      real*8  dzdZ_n(Nx,Ny,Nz)
      real*8  dzdZ2_n(Nx,Ny,Nz)
      real*8  phi1_np1(Nx,Ny,Nz)
      real*8  phi1_n(Nx,Ny,Nz)
      real*8  phi2_np1(Nx,Ny,Nz)
      real*8  phi2_n(Nx,Ny,Nz)
      real*8  pi1_np1(Nx,Ny,Nz)
      real*8  pi1_n(Nx,Ny,Nz)
      real*8  pi2_np1(Nx,Ny,Nz)
      real*8  pi2_n(Nx,Ny,Nz)
      real*8  at_np1(Nx,Ny,Nz)
      real*8  at_n(Nx,Ny,Nz)
      real*8  ax_np1(Nx,Ny,Nz)
      real*8  ax_n(Nx,Ny,Nz)
      real*8  ay_np1(Nx,Ny,Nz)
      real*8  ay_n(Nx,Ny,Nz)
      real*8  az_np1(Nx,Ny,Nz)
      real*8  az_n(Nx,Ny,Nz)
      real*8  bt_np1(Nx,Ny,Nz)
      real*8  bt_n(Nx,Ny,Nz)
      real*8  bx_np1(Nx,Ny,Nz)
      real*8  bx_n(Nx,Ny,Nz)
      real*8  by_np1(Nx,Ny,Nz)
      real*8  by_n(Nx,Ny,Nz)
      real*8  bz_np1(Nx,Ny,Nz)
      real*8  bz_n(Nx,Ny,Nz)
      real*8  x(*)
      real*8  y(*)
      real*8  z(*)
      real*8  dx
      real*8  dy
      real*8  dz
      real*8  xmax
      real*8  xmin
      real*8  ymin
      real*8  ymax
      real*8  zmin
      real*8  zmax
      real*8  cx
      real*8  cy
      real*8  cz
      real*8  cx1
      real*8  cx2
      real*8  cy1
      real*8  cy2
      real*8  cz1
      real*8  cz2
      real*8  w, w2
      real*8  v, v2
      real*8  phase
      integer anti
      integer ini_type
      real*8  e

      real*8      pi
      parameter ( pi = 3.14159265358979d0 )

      integer unum, stat  ! used to read/open/close a file
      character(len=100) :: datafile = "../../run/input/initdata"

      real*8, allocatable :: radial_axis(:)
      real*8, allocatable :: Q_shape(:)     ! stores pre-computed shape from file
      real*8, allocatable :: A_shape(:)
      real*8, allocatable :: phic(:,:,:)    ! stores Lorentz-contracted shape (for boosted case)
      real*8, allocatable :: Ac(:,:,:)
      real*8, allocatable :: dphic(:,:,:)   ! stores derivative of contracted shape (for boosted case)
      real*8, allocatable :: dAc(:,:,:)

      real*8, allocatable :: phi1_n1(:,:,:)
      real*8, allocatable :: phi2_n1(:,:,:)
      real*8, allocatable :: pi1_n1(:,:,:)
      real*8, allocatable :: pi2_n1(:,:,:)
      real*8, allocatable :: at_n1(:,:,:)
      real*8, allocatable :: az_n1(:,:,:)
      real*8, allocatable :: bt_n1(:,:,:)
      real*8, allocatable :: bz_n1(:,:,:)

      real*8, allocatable :: phi1_n2(:,:,:)
      real*8, allocatable :: phi2_n2(:,:,:)
      real*8, allocatable :: pi1_n2(:,:,:)
      real*8, allocatable :: pi2_n2(:,:,:)
      real*8, allocatable :: at_n2(:,:,:)
      real*8, allocatable :: az_n2(:,:,:)
      real*8, allocatable :: bt_n2(:,:,:)
      real*8, allocatable :: bz_n2(:,:,:)

      real*8, dimension(0:4) :: xx=0, yy=0  ! storage for passing to neville()
      real*8  rr, drr, Q, A
      integer length       ! pre-computed length of radial axis
      integer rrindex      ! array index closest to a point along radial axis
      real*8  gamma        ! Lorentz factor
      real*8  phasef

      gamma  = 1d0/sqrt(1d0-v**2)
      phasef = phase*pi

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
     &     action='read',iostat=stat)
      do i=1,length
       read (unum,*) rr, Q, A
       radial_axis(i)=rr
       Q_shape(i)=Q
       A_shape(i)=A
       if (stat .ne. 0) exit
      end do
      close(unit=unum, status='keep', iostat=stat)

      ! assumes initial data is output on uniformly-spaced intervals
      drr=radial_axis(2)-radial_axis(1)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          SINGLE STATIC Q-BALL           cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (ini_type .eq. 1) then

      ! interpolate
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz

         rr = sqrt((Xx_n(i,j,k)-cx)**2d0+(Yy_n(i,j,k)-cy)**2d0+
     &             (Zz_n(i,j,k)-cz)**2d0)

         ! retrieve array index closest to rr
         rrindex=idnint(rr/drr)+1
         if ( rrindex+2 .ge. length ) then
           !STOP 'Interpolation error (length exceeded)'
           phi1_n(i,j,k) = 0d0
           at_n(i,j,k) = 0d0
           cycle
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
         phi1_n(i,j,k)=neville(xx,yy,rr)

         yy(0)=A_shape(rrindex-2)
         yy(1)=A_shape(rrindex-1)
         yy(2)=A_shape(rrindex)
         yy(3)=A_shape(rrindex+1)
         yy(4)=A_shape(rrindex+2)
         at_n(i,j,k)=neville(xx,yy,rr)

      end do
      end do
      end do

      ! specify initial data
      phi2_n(:,:,:) = 0d0
      pi1_n(:,:,:)  = 0d0
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
        pi2_n(i,j,k) = w*phi1_n(i,j,k)
      end do
      end do
      end do
      ax_n(:,:,:) = 0d0
      ay_n(:,:,:) = 0d0
      az_n(:,:,:) = 0d0
      bx_n(:,:,:) = 0d0
      by_n(:,:,:) = 0d0
      bz_n(:,:,:) = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          SINGLE BOOSTED Q-BALL          cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 2 ) then

      allocate(phic(Nx,Ny,Nz))
      allocate(dphic(Nx,Ny,Nz))
      allocate(Ac(Nx,Ny,Nz))
      allocate(dAc(Nx,Ny,Nz))

      ! interpolate
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz

         rr = sqrt((Xx_n(i,j,k)-cx)**2d0+(Yy_n(i,j,k)-cy)**2d0+
     &             (gamma*(Zz_n(i,j,k)-cz))**2d0)

         ! retrieve array index closest to rr
         rrindex=idnint(rr/drr)+1
         if ( rrindex+2 .ge. length ) then
           !STOP 'Interpolation error (length exceeded)'
           phic(i,j,k) = 0d0
           cycle
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

         phic(i,j,k)=neville(xx,yy,rr)

         yy(0)=A_shape(rrindex-2)
         yy(1)=A_shape(rrindex-1)
         yy(2)=A_shape(rrindex)
         yy(3)=A_shape(rrindex+1)
         yy(4)=A_shape(rrindex+2)

         Ac(i,j,k)=neville(xx,yy,rr)

      end do
      end do
      end do

      ! compute derivatives using O(h^4) differencing
      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(-(25*phic(i,j,k))/12 +
     &4*phic(i,j,k+1) - 3*phic(i,j,k+2) + (4*phic(i,j,k+3))/3 -
     &phic(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(phic(i,j,k-2)/12 - (2*phic(i,j,
     &k-1))/3 + (2*phic(i,j,k+1))/3 - phic(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dphic(i,j,k) = dzdZ_n(i,j,k)*((25*phic(i,j,k))/12 -
     &4*phic(i,j,k-1) + 3*phic(i,j,k-2) - (4*phic(i,j,k-3))/3 +
     &phic(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(-(25*Ac(i,j,k))/12 +
     &4*Ac(i,j,k+1) - 3*Ac(i,j,k+2) + (4*Ac(i,j,k+3))/3 -
     &Ac(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(Ac(i,j,k-2)/12 - (2*Ac(i,j,
     &k-1))/3 + (2*Ac(i,j,k+1))/3 - Ac(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dAc(i,j,k) = dzdZ_n(i,j,k)*((25*Ac(i,j,k))/12 -
     &4*Ac(i,j,k-1) + 3*Ac(i,j,k-2) - (4*Ac(i,j,k-3))/3 +
     &Ac(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      ! specify initial data
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         phi1_n(i,j,k) = phic(i,j,k)*cos(w*gamma*v*(Zz_n(i,j,k)-cz))
         phi2_n(i,j,k) = phic(i,j,k)*sin(w*gamma*v*(Zz_n(i,j,k)-cz))
         pi1_n(i,j,k)  = dphic(i,j,k)*cos(w*gamma*v*
     &                 (Zz_n(i,j,k)-cz))*gamma*v - w*gamma*phi2_n(i,j,k)
         pi2_n(i,j,k)  = dphic(i,j,k)*sin(w*gamma*v*
     &                 (Zz_n(i,j,k)-cz))*gamma*v + w*gamma*phi1_n(i,j,k)
      end do
      end do
      end do

      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         at_n(i,j,k) =   gamma*Ac(i,j,k)
         az_n(i,j,k) = v*gamma*Ac(i,j,k)
         bt_n(i,j,k) = v*(gamma**2.0d0)*dAc(i,j,k)
         bz_n(i,j,k) = (v**2.0d0)*(gamma**2.0d0)*dAc(i,j,k)
      enddo
      enddo
      enddo

      ax_n(:,:,:) = 0d0
      ay_n(:,:,:) = 0d0
      bx_n(:,:,:) = 0d0
      by_n(:,:,:) = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc           TWO BOOSTED Q-BALLS           cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 3 ) then

      allocate(phic(Nx,Ny,Nz))
      allocate(dphic(Nx,Ny,Nz))
      allocate(Ac(Nx,Ny,Nz))
      allocate(dAc(Nx,Ny,Nz))

      ! Q-ball #1 setup
      allocate(phi1_n1(Nx,Ny,Nz))
      allocate(phi2_n1(Nx,Ny,Nz))
      allocate(pi1_n1(Nx,Ny,Nz))
      allocate(pi2_n1(Nx,Ny,Nz))
      allocate(at_n1(Nx,Ny,Nz))
      allocate(az_n1(Nx,Ny,Nz))
      allocate(bt_n1(Nx,Ny,Nz))
      allocate(bz_n1(Nx,Ny,Nz))

      cx = cx1
      cy = cy1
      cz = cz1

      ! interpolate
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz

         rr = sqrt((Xx_n(i,j,k)-cx)**2d0+(Yy_n(i,j,k)-cy)**2d0+
     &             (gamma*(Zz_n(i,j,k)-cz))**2d0)

         ! retrieve array index closest to rr
         rrindex=idnint(rr/drr)+1
         if ( rrindex+2 .ge. length ) then
           !STOP 'Interpolation error (length exceeded)'
           phic(i,j,k) = 0d0
           cycle
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

         phic(i,j,k)=neville(xx,yy,rr)

         yy(0)=A_shape(rrindex-2)
         yy(1)=A_shape(rrindex-1)
         yy(2)=A_shape(rrindex)
         yy(3)=A_shape(rrindex+1)
         yy(4)=A_shape(rrindex+2)

         Ac(i,j,k)=neville(xx,yy,rr)

      end do
      end do
      end do

      ! compute derivatives using O(h^4) differencing
      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(-(25*phic(i,j,k))/12 +
     &4*phic(i,j,k+1) - 3*phic(i,j,k+2) + (4*phic(i,j,k+3))/3 -
     &phic(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(phic(i,j,k-2)/12 - (2*phic(i,j,
     &k-1))/3 + (2*phic(i,j,k+1))/3 - phic(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dphic(i,j,k) = dzdZ_n(i,j,k)*((25*phic(i,j,k))/12 -
     &4*phic(i,j,k-1) + 3*phic(i,j,k-2) - (4*phic(i,j,k-3))/3 +
     &phic(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(-(25*Ac(i,j,k))/12 +
     &4*Ac(i,j,k+1) - 3*Ac(i,j,k+2) + (4*Ac(i,j,k+3))/3 -
     &Ac(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(Ac(i,j,k-2)/12 - (2*Ac(i,j,
     &k-1))/3 + (2*Ac(i,j,k+1))/3 - Ac(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dAc(i,j,k) = dzdZ_n(i,j,k)*((25*Ac(i,j,k))/12 -
     &4*Ac(i,j,k-1) + 3*Ac(i,j,k-2) - (4*Ac(i,j,k-3))/3 +
     &Ac(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      ! specify initial data
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         phi1_n1(i,j,k) = phic(i,j,k)*cos(w*gamma*v*(Zz_n(i,j,k)-cz))
         phi2_n1(i,j,k) = phic(i,j,k)*sin(w*gamma*v*(Zz_n(i,j,k)-cz))
         pi1_n1(i,j,k)  = dphic(i,j,k)*cos(w*gamma*v*
     &                (Zz_n(i,j,k)-cz))*gamma*v - w*gamma*phi2_n1(i,j,k)
         pi2_n1(i,j,k)  = dphic(i,j,k)*sin(w*gamma*v*
     &                (Zz_n(i,j,k)-cz))*gamma*v + w*gamma*phi1_n1(i,j,k)
      end do
      end do
      end do

      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         at_n1(i,j,k) =   gamma*Ac(i,j,k)
         az_n1(i,j,k) = v*gamma*Ac(i,j,k)
         bt_n1(i,j,k) = v*(gamma**2.0d0)*dAc(i,j,k)
         bz_n1(i,j,k) = (v**2.0d0)*(gamma**2.0d0)*dAc(i,j,k)
      enddo
      enddo
      enddo

      ! Q-ball #2 setup
      allocate(phi1_n2(Nx,Ny,Nz))
      allocate(phi2_n2(Nx,Ny,Nz))
      allocate(pi1_n2(Nx,Ny,Nz))
      allocate(pi2_n2(Nx,Ny,Nz))
      allocate(at_n2(Nx,Ny,Nz))
      allocate(az_n2(Nx,Ny,Nz))
      allocate(bt_n2(Nx,Ny,Nz))
      allocate(bz_n2(Nx,Ny,Nz))

      cx = cx2
      cy = cy2
      cz = cz2

      v2 = -v

      ! interpolate
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz

         rr = sqrt((Xx_n(i,j,k)-cx)**2d0+(Yy_n(i,j,k)-cy)**2d0+
     &             (gamma*(Zz_n(i,j,k)-cz))**2d0)

         ! retrieve array index closest to rr
         rrindex=idnint(rr/drr)+1
         if ( rrindex+2 .ge. length ) then
           !STOP 'Interpolation error (length exceeded)'
           phic(i,j,k) = 0d0
           cycle
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

         phic(i,j,k)=neville(xx,yy,rr)

         yy(0)=A_shape(rrindex-2)
         yy(1)=A_shape(rrindex-1)
         yy(2)=A_shape(rrindex)
         yy(3)=A_shape(rrindex+1)
         yy(4)=A_shape(rrindex+2)

         Ac(i,j,k)=neville(xx,yy,rr)

      end do
      end do
      end do

      ! Construct anti-Q-ball
      if ( anti .eq. 1 ) then
        w2 = -1d0*w
        do i=1,Nx
        do j=1,Ny
        do k=1,Nz
         Ac(i,j,k)=-1d0*Ac(i,j,k)
        end do
        end do
        end do
      else
        w2 = w
      end if
      
      ! compute derivatives using O(h^4) differencing
      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(-(25*phic(i,j,k))/12 +
     &4*phic(i,j,k+1) - 3*phic(i,j,k+2) + (4*phic(i,j,k+3))/3 -
     &phic(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dphic(i,j,k) = dzdZ_n(i,j,k)*(phic(i,j,k-2)/12 - (2*phic(i,j,
     &k-1))/3 + (2*phic(i,j,k+1))/3 - phic(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dphic(i,j,k) = dzdZ_n(i,j,k)*((25*phic(i,j,k))/12 -
     &4*phic(i,j,k-1) + 3*phic(i,j,k-2) - (4*phic(i,j,k-3))/3 +
     &phic(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      do i=1,Nx
      do j=1,Ny
        do k=1,2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(-(25*Ac(i,j,k))/12 +
     &4*Ac(i,j,k+1) - 3*Ac(i,j,k+2) + (4*Ac(i,j,k+3))/3 -
     &Ac(i,j,k+4)/4)/(dz*gamma)
        end do

        do k=3,Nz-2
          dAc(i,j,k) = dzdZ_n(i,j,k)*(Ac(i,j,k-2)/12 - (2*Ac(i,j,
     &k-1))/3 + (2*Ac(i,j,k+1))/3 - Ac(i,j,k+2)/12)/(dz*gamma)
        end do

        do k=Nz-1,Nz
          dAc(i,j,k) = dzdZ_n(i,j,k)*((25*Ac(i,j,k))/12 -
     &4*Ac(i,j,k-1) + 3*Ac(i,j,k-2) - (4*Ac(i,j,k-3))/3 +
     &Ac(i,j,k-4)/4)/(dz*gamma)
        end do
      end do
      end do

      ! specify initial data
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         phi1_n2(i,j,k) = phic(i,j,k)*cos(w2*gamma*v2*(Zz_n(i,j,k)-cz)
     &                    +phasef)
         phi2_n2(i,j,k) = phic(i,j,k)*sin(w2*gamma*v2*(Zz_n(i,j,k)-cz)
     &                    +phasef)
         pi1_n2(i,j,k) = dphic(i,j,k)*cos(w2*gamma*v2*(Zz_n(i,j,k)-cz)
     &                   +phasef)*gamma*v2 - w2*gamma*phi2_n2(i,j,k)
         pi2_n2(i,j,k) = dphic(i,j,k)*sin(w2*gamma*v2*(Zz_n(i,j,k)-cz)
     &                   +phasef)*gamma*v2 + w2*gamma*phi1_n2(i,j,k)
      end do
      end do
      end do
        
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         at_n2(i,j,k) =   gamma*Ac(i,j,k)
         az_n2(i,j,k) = v2*gamma*Ac(i,j,k)
         bt_n2(i,j,k) = v2*(gamma**2.0d0)*dAc(i,j,k)
         bz_n2(i,j,k) = (v2**2.0d0)*(gamma**2.0d0)*dAc(i,j,k)
      enddo
      enddo
      enddo

      ! combine Q-ball #1 and Q-ball #2 data
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
        phi1_n(i,j,k) = phi1_n1(i,j,k) + phi1_n2(i,j,k)
        phi2_n(i,j,k) = phi2_n1(i,j,k) + phi2_n2(i,j,k)
        pi1_n(i,j,k)  = pi1_n1(i,j,k)  + pi1_n2(i,j,k)
        pi2_n(i,j,k)  = pi2_n1(i,j,k)  + pi2_n2(i,j,k)
        at_n(i,j,k)   = at_n1(i,j,k)  + at_n2(i,j,k)
        az_n(i,j,k)   = az_n1(i,j,k)  + az_n2(i,j,k)
        bt_n(i,j,k)   = bt_n1(i,j,k)  + bt_n2(i,j,k)
        bz_n(i,j,k)   = bz_n1(i,j,k)  + bz_n2(i,j,k)
      end do
      end do
      end do

      ax_n(:,:,:) = 0d0
      ay_n(:,:,:) = 0d0
      bx_n(:,:,:) = 0d0
      by_n(:,:,:) = 0d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccc          GAUSSIAN INITIAL DATA          cccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      elseif ( ini_type .eq. 4 ) then

      ! set arbitrary field values
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz

      at_n(i,j,k) =
     #0.3D0 * exp(-0.4000000000D0 * Xx_n(i, j, k) ** 2 - 0.40000000
     #00D0 * Yy_n(i, j, k) ** 2 - 0.4000000000D0 * Zz_n(i, j, k) ** 2)

      ax_n(i,j,k) =
     #0.25D0 * exp(-0.4347826087D0 * Xx_n(i, j, k) ** 2 - 0.434782
     #6087D0 * Yy_n(i, j, k) ** 2 - 0.4347826087D0 * Zz_n(i, j, k) ** 2)

      ay_n(i,j,k) =
     #0.35D0 * exp(-0.3846153846D0 * Xx_n(i, j, k) ** 2 - 0.384615
     #3846D0 * Yy_n(i, j, k) ** 2 - 0.3846153846D0 * Zz_n(i, j, k) ** 2)

      az_n(i,j,k) =
     #0.3D0 * exp(-0.4166666667D0 * Xx_n(i, j, k) ** 2 - 0.4166666
     #667D0 * Yy_n(i, j, k) ** 2 - 0.4166666667D0 * Zz_n(i, j, k) ** 2)

      bt_n(i,j,k) =
     #-0.2173913044D0 * Xx_n(i, j, k) * exp(-0.4347826087D0 * Xx_n
     #(i, j, k) ** 2 - 0.4347826087D0 * Yy_n(i, j, k) ** 2 - 0.434782608
     #7D0 * Zz_n(i, j, k) ** 2) - 0.2692307692D0 * Yy_n(i, j, k) * exp(-
     #0.3846153846D0 * Xx_n(i, j, k) ** 2 - 0.3846153846D0 * Yy_n(i, j, 
     #k) ** 2 - 0.3846153846D0 * Zz_n(i, j, k) ** 2) - 0.2500000000D0 * 
     #Zz_n(i, j, k) * exp(-0.4166666667D0 * Xx_n(i, j, k) ** 2 - 0.41666
     #66667D0 * Yy_n(i, j, k) ** 2 - 0.4166666667D0 * Zz_n(i, j, k) ** 2
     #)

      bx_n(i,j,k) =
     #-0.2400000000D0 * Xx_n(i, j, k) * exp(-0.4000000000D0 * Xx_n
     #(i, j, k) ** 2 - 0.4000000000D0 * Yy_n(i, j, k) ** 2 - 0.400000000
     #0D0 * Zz_n(i, j, k) ** 2)

      by_n(i,j,k) =
     #-0.2400000000D0 * Yy_n(i, j, k) * exp(-0.4000000000D0 * Xx_n
     #(i, j, k) ** 2 - 0.4000000000D0 * Yy_n(i, j, k) ** 2 - 0.400000000
     #0D0 * Zz_n(i, j, k) ** 2)

      bz_n(i,j,k) =
     #-0.2400000000D0 * Zz_n(i, j, k) * exp(-0.4000000000D0 * Xx_n
     #(i, j, k) ** 2 - 0.4000000000D0 * Yy_n(i, j, k) ** 2 - 0.400000000
     #0D0 * Zz_n(i, j, k) ** 2)

      phi1_n(i,j,k) = 0d0

      phi2_n(i,j,k) = 0d0

      pi1_n(i,j,k) =
     #0.325D0 * exp(-0.4000000000D0 * (Xx_n(i, j, k) - 0.20D1) **
     #2 - 0.4000000000D0 * Yy_n(i, j, k) ** 2 - 0.4000000000D0 * Zz_n(i,
     # j, k) ** 2)

      pi2_n(i,j,k) =
     #0.275D0 * exp(-0.4000000000D0 * (Xx_n(i, j, k) - 0.20D1) **
     # 2 - 0.4000000000D0 * Yy_n(i, j, k) ** 2 - 0.4000000000D0 * Zz_n(i
     #, j, k) ** 2)
 
      end do
      end do
      end do
 
      end if

      end
