      module globals
      implicit none

      integer counter
      integer N1, N2, N3
      real*8 Q, sigma, pi
      parameter ( Q = 1.0d0,   sigma = 1.0d0 )
      parameter ( pi = 3.14159265358979d0 )
      parameter ( N1 = 257, N2 = 257, N3 = 257 )
c      parameter ( N1 = 129, N2 = 129, N3 = 129 )
c      parameter ( N1 = 65, N2 = 65, N3 = 65 )

      end module

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     initguess: initializes a guess for the multigrid solve
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initguess(V,x,y,z,Nx,Ny,Nz,Xx,Yy,Zz,At)
      use globals

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz)
      real*8 X(Nx), Y(Ny), Z(Nz)
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 At(Nx,Ny,Nz)
      real*8 rand

      integer i, j, k

      counter = 0

      ! Initial guess
      do i=2,Nx-1
      do j=2,Ny-1
      do k=2,Nz-1
c        call RANDOM_NUMBER(rand)
c        V(i,j,k)=rand/4d0
c        V(i,j,k)=0.25d0
        V(i,j,k)=At(i,j,k)
      end do
      end do
      end do
      
      ! note that we set boundary conditions at two outermost points to
      ! facilitate fourth-order stencil used for defect correction

c     boundary condition @ x=xmin
      do j=1,Ny
      do k=1,Nz
        V(1:2,j,k) = At(1:2,j,k)
      end do
      end do

c     boundary condition @ x=xmax
      do j=1,Ny
      do k=1,Nz
        V(Nx-1:Nx,j,k) = At(Nx-1:Nx,j,k) 
      end do
      end do

c     boundary condition @ y=ymin
      do i=1,Nx
      do k=1,Nz
        V(i,1:2,k) = At(i,1:2,k)
      end do
      end do

c     boundary condition @ y=ymax
      do i=1,Nx
      do k=1,Nz
        V(i,Ny-1:Ny,k) = At(i,Ny-1:Ny,k)
      end do
      end do

c     boundary condition @ z=zmin
      do i=1,Nx
      do j=1,Ny
        V(i,j,1:2) = At(i,j,1:2)
      end do
      end do

c     boundary condition @ z=zmax
      do i=1,Nx
      do j=1,Ny
        V(i,j,Nz-1:Nz) = At(i,j,Nz-1:Nz)
      end do
      end do

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init_qball2: initializes gauged Q-ball fields using data from the
c                  multigrid solve; reads in data using rff() function
c                  and interpolates using interpd() function
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_qball2(V,x,y,z,Nx,Ny,Nz,dx,dy,dz,Xx,Yy,Zz,
     &dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,
     &At,Ax,Ay,Az,Bt,Bx,By,Bz,
     &xmin,xmax,ymin,ymax,zmin,zmax)
      use globals
      use interp_mod

      implicit none

      integer Nx, Ny, Nz
      real*8 dx, dy, dz
      real*8 V(Nx,Ny,Nz)
      real*8 x(Nx), y(Ny), z(Nz)
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 dxdX(Nx,Ny,Nz), dxdX2(Nx,Ny,Nz)
      real*8 dydY(Nx,Ny,Nz), dydY2(Nx,Ny,Nz)
      real*8 dzdZ(Nx,Ny,Nz), dzdZ2(Nx,Ny,Nz)
      real*8 At(Nx,Ny,Nz), Bt(Nx,Ny,Nz)
      real*8 Ax(Nx,Ny,Nz), Bx(Nx,Ny,Nz)
      real*8 Ay(Nx,Ny,Nz), By(Nx,Ny,Nz)
      real*8 Az(Nx,Ny,Nz), Bz(Nx,Ny,Nz)
      real*8 xmin,xmax,ymin,ymax,zmin,zmax

      integer i, j, k
      real*8 u(N1,N2,N3)
      real*8 d1, d2, d3

      call rff(u,N1,N2,N3)  ! N1, N2, N3 needs to match actual size used in V.bin

      d1 = (xmax-xmin)/real(N1-1)
      d2 = (ymax-ymin)/real(N2-1)
      d3 = (zmax-zmin)/real(N3-1)

      ! interpolate multigrid solution
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
        V(i,j,k) = interpd(N1, N2, N3, d1, d2, d3, xmin, xmax,
     &                     ymin, ymax, zmin, zmax, x(i), y(j), z(k),
     &                     u)
      end do
      end do
      end do

      ! initialize fields to fourth-order using interpolated values
      do i=3,Nx-2 
      do j=3,Ny-2
      do k=3,Nz-2

      Bt(i, j, k) = 
     #dble((dzdZ(i, j, k) * dx * dy * Az(i, j, k - 2) - dzdZ(i, j,
     # k) * dx * dy * Az(i, j, k + 2) + 8 * dzdZ(i, j, k) * dx * dy * Az
     #(i, j, k + 1) - 8 * dzdZ(i, j, k) * dx * dy * Az(i, j, k - 1) + dy
     #dY(i, j, k) * dx * dz * Ay(i, j - 2, k) - dydY(i, j, k) * dx * dz 
     #* Ay(i, j + 2, k) + 8 * dydY(i, j, k) * dx * dz * Ay(i, j + 1, k) 
     #- 8 * dydY(i, j, k) * dx * dz * Ay(i, j - 1, k) + dxdX(i, j, k) * 
     #dy * dz * Ax(i - 2, j, k) - dxdX(i, j, k) * dy * dz * Ax(i + 2, j,
     # k) + 8 * dxdX(i, j, k) * dy * dz * Ax(i + 1, j, k) - 8 * dxdX(i, 
     #j, k) * dy * dz * Ax(i - 1, j, k)) / dx / dy / dz) / 0.12D2

      Bx(i, j, k) =
     #dble(dxdX(i, j, k) * (At(i - 2, j, k) - 8 * At(i - 1, j, k)
     # + 8 * At(i + 1, j, k) - At(i + 2, j, k) - V(i - 2, j, k) + 8 * V(
     #i - 1, j, k) - 8 * V(i + 1, j, k) + V(i + 2, j, k)) / dx) / 0.12D2

      By(i, j, k) = 
     #-dble(dydY(i, j, k) * (-At(i, j - 2, k) + 8 * At(i, j - 1, 
     #k) - 8 * At(i, j + 1, k) + At(i, j + 2, k) + V(i, j - 2, k) - 8 * 
     #V(i, j - 1, k) + 8 * V(i, j + 1, k) - V(i, j + 2, k)) / dy) / 0.12
     #D2

      Bz(i, j, k) = 
     #dble(dzdZ(i, j, k) * (At(i, j, k - 2) - 8 * At(i, j, k - 1)
     # + 8 * At(i, j, k + 1) - At(i, j, k + 2) - V(i, j, k - 2) + 8 * V(
     #i, j, k - 1) - 8 * V(i, j, k + 1) + V(i, j, k + 2)) / dz) / 0.12D2

      end do
      end do
      end do

      ! hacky fix for AMR boundaries
      do i=2,1,-1
      do j=2,Ny-1
      do k=2,Nz-1
        Bt(i,j,k) =  Bt(i+1,j,k)
        Bx(i,j,k) =  Bx(i+1,j,k)
        By(i,j,k) =  By(i+1,j,k)
        Bz(i,j,k) =  Bz(i+1,j,k)
      end do
      end do
      end do
      
      do i=Nx-1,Nx
      do j=2,Ny-1
      do k=2,Nz-1
        Bt(i,j,k) =  Bt(i-1,j,k)
        Bx(i,j,k) =  Bx(i-1,j,k)
        By(i,j,k) =  By(i-1,j,k)
        Bz(i,j,k) =  Bz(i-1,j,k)
      end do
      end do
      end do
      
      do j=2,1,-1
      do i=2,Nx-1
      do k=2,Nz-1
        Bt(i,j,k) =  Bt(i,j+1,k)
        Bx(i,j,k) =  Bx(i,j+1,k)
        By(i,j,k) =  By(i,j+1,k)
        Bz(i,j,k) =  Bz(i,j+1,k)
      end do
      end do
      end do
      
      do j=Ny-1,Ny
      do i=2,Nx-1
      do k=2,Nz-1
        Bt(i,j,k) =  Bt(i,j-1,k)
        Bx(i,j,k) =  Bx(i,j-1,k)
        By(i,j,k) =  By(i,j-1,k)
        Bz(i,j,k) =  Bz(i,j-1,k)
      end do
      end do
      end do
      
      do k=2,1,-1
      do i=2,Nx-1
      do j=2,Ny-1
        Bt(i,j,k) =  Bt(i,j,k+1)
        Bx(i,j,k) =  Bx(i,j,k+1)
        By(i,j,k) =  By(i,j,k+1)
        Bz(i,j,k) =  Bz(i,j,k+1)
      end do
      end do
      end do
      
      do k=Nz-1,Nz
      do i=2,Nx-1
      do j=2,Ny-1
        Bt(i,j,k) =  Bt(i,j,k-1)
        Bx(i,j,k) =  Bx(i,j,k-1)
        By(i,j,k) =  By(i,j,k-1)
        Bz(i,j,k) =  Bz(i,j,k-1)
      end do
      end do
      end do

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     lop: differential operator L(V)=0 computed where cmask=CMASK_ON
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lop(LV,V,defect,cmask,x,y,z,Nx,Ny,Nz,dx,dy,dz,Xx,Yy,Zz,
     &dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,e,phi1,phi2,pi1,pi2,modphi,At)
      use globals

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz)
      real*8 defect(Nx,Ny,Nz)
      real*8 cmask(Nx,Ny,Nz), LV(Nx,Ny,Nz)
      real*8 x(Nx), y(Ny), z(Nz)
      real*8 dx, dy, dz
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 dXdx(Nx,Ny,Nz), dXdx2(Nx,Ny,Nz)
      real*8 dYdy(Nx,Ny,Nz), dYdy2(Nx,Ny,Nz)
      real*8 dZdz(Nx,Ny,Nz), dZdz2(Nx,Ny,Nz)

      integer i,j,k
      real*8 f(Nx,Ny,Nz)

      real*8 e
      real*8 phi1(Nx,Ny,Nz), phi2(Nx,Ny,Nz)
      real*8 pi1(Nx,Ny,Nz), pi2(Nx,Ny,Nz)
      real*8 modphi(Nx,Ny,Nz), At(Nx,Ny,Nz)

      include 'cmask.inc'

      ! source term
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         f(i,j,k)=2d0*e*(phi2(i,j,k)*pi1(i,j,k)-phi1(i,j,k)*pi2(i,j,k))
     &          +2d0*e**2d0*At(i,j,k)*modphi(i,j,k)**2d0
      end do
      end do
      end do

c     interior points
      do i=2,Nx-1
      do j=2,Ny-1
      do k=2,Nz-1
        if (cmask(i,j,k) .eq. CMASK_ON) then
          LV(i,j,k)=
     #dble(1 / dx ** 2 * (V(i - 1, j, k) - 2 * V(i, j, k) + V(i + 1
     #, j, k)) * dxdX(i, j, k) ** 2) + dble((-V(i - 1, j, k) + V(i + 1, 
     #j, k)) / dx * dxdX2(i, j, k)) / 0.2D1 + dble(1 / dy ** 2 * (V(i, j
     # - 1, k) - 2 * V(i, j, k) + V(i, j + 1, k)) * dydY(i, j, k) ** 2) 
     #+ dble((-V(i, j - 1, k) + V(i, j + 1, k)) / dy * dydY2(i, j, k)) /
     # 0.2D1 + dble(1 / dz ** 2 * (V(i, j, k - 1) - 2 * V(i, j, k) + V(i
     #, j, k + 1)) * dzdZ(i, j, k) ** 2) + dble((-V(i, j, k - 1) + V(i, 
     #j, k + 1)) / dz * dzdZ2(i, j, k)) / 0.2D1 - defect(i, j, k)
        end if
      end do
      end do
      end do

c      write(*,*) counter, "lop:     ", Nx, Ny, Nz

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     residual: residual L[V]-rhs computed where cmask=CMASK_ON 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine residual(res,rhs,V,defect,cmask,x,y,z,norm,Nx,Ny,Nz,
     &dx,dy,dz,Xx,Yy,Zz,dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,e,
     &phi1,phi2,pi1,pi2,modphi,At)
      use globals

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz)
      real*8 defect(Nx,Ny,Nz)
      real*8 cmask(Nx,Ny,Nz), res(Nx,Ny,Nz), rhs(Nx,Ny,Nz)
      real*8 x(Nx), y(Ny), z(Nz), norm
      real*8 dx, dy, dz
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 dxdX(Nx,Ny,Nz), dxdX2(Nx,Ny,Nz)
      real*8 dydY(Nx,Ny,Nz), dydY2(Nx,Ny,Nz)
      real*8 dzdZ(Nx,Ny,Nz), dzdZ2(Nx,Ny,Nz)
      integer i, j, k, sum

      real*8 e
      real*8 phi1(Nx,Ny,Nz), phi2(Nx,Ny,Nz)
      real*8 pi1(Nx,Ny,Nz), pi2(Nx,Ny,Nz)
      real*8 modphi(Nx,Ny,Nz), At(Nx,Ny,Nz)

      include 'cmask.inc'

      call lop(res,V,defect,cmask,x,y,z,Nx,Ny,Nz,dx,dy,dz,Xx,Yy,Zz,
     &dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,e,phi1,phi2,pi1,pi2,modphi,At)

      norm = 0d0
      sum  = 0

      do i=2,Nx-1
      do j=2,Ny-1
      do k=2,Nz-1
        if (cmask(i,j,k) .eq. CMASK_ON) then
          res(i,j,k) = res(i,j,k)-rhs(i,j,k)
          norm       = norm+res(i,j,k)**2
          sum        = sum+1
        end if
      end do
      end do
      end do

      norm = sqrt(norm/sum)

c      write(*,*) counter, "residual:", Nx, Ny, Nz, norm

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     relax: applies alternating-direction zebra line relaxation using
c            LAPACK 'dgtsv' solver
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine relax(V,V_rhs,defect,cmask,phys_bdy,norm,x,y,z,
     &Nx,Ny,Nz,dx,dy,dz,Xx,Yy,Zz,dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,
     &e,phi1,phi2,pi1,pi2,modphi,At)
      use globals

      implicit none

      integer Nx, Ny, Nz
      real*8 V(Nx,Ny,Nz), V_rhs(Nx,Ny,Nz), cmask(Nx,Ny,Nz)
      real*8 defect(Nx,Ny,Nz)
      integer phys_bdy(6)
      real*8 norm
      real*8 dx, dy, dz
      real*8 x(Nx), y(Ny), z(Nz)
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 dxdX(Nx,Ny,Nz), dxdX2(Nx,Ny,Nz)
      real*8 dydY(Nx,Ny,Nz), dydY2(Nx,Ny,Nz)
      real*8 dzdZ(Nx,Ny,Nz), dzdZ2(Nx,Ny,Nz)

      integer i, j, k, l
      real*8, dimension (Nx) :: d, du, dl, rhs
      integer nrhs, info
      real*8, dimension (Nx, Ny, Nz) :: f
      real*8 res
      integer sum

      real*8 e
      real*8 phi1(Nx,Ny,Nz), phi2(Nx,Ny,Nz)
      real*8 pi1(Nx,Ny,Nz), pi2(Nx,Ny,Nz)
      real*8 modphi(Nx,Ny,Nz), At(Nx,Ny,Nz)

      include 'cmask.inc'

      norm = 0d0
      sum  = 0

      ! source term
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
        f(i,j,k)=2d0*e*(phi2(i,j,k)*pi1(i,j,k)-phi1(i,j,k)*pi2(i,j,k))
     &          +2d0*e**2d0*At(i,j,k)*modphi(i,j,k)**2d0
      end do
      end do
      end do

c     Z-DIRECTION LINE RELAXATION
c     set up tridiagonal system; note that indexing on 
c     lower diagonal is always (j-1) when implementing the 
c     jth equation
      do i=3,Nx-2
      do j=3,Ny-2
        d(1)   = 1.0d0
        du(1)  = 0.0d0
        rhs(1) = V(i,j,1)
        dl(1)  = 0.0d0
        d(2)   = 1.0d0
        du(2)  = 0.0d0
        rhs(2) = V(i,j,2)
        do k=3,Nz-2  ! using 7-point stencil
          dl(k-1) =
     #0.1D1 / dz ** 2 * dzdZ(i, j, k) ** 2 - 0.1D1 / dz * dzdZ2(i,
     # j, k) / 0.2D1
          d(k)    =
     #-2 / dx ** 2 * dxdX(i, j, k) ** 2 - 2 / dy ** 2 * dydY(i, j
     #, k) ** 2 - 2 / dz ** 2 * dzdZ(i, j, k) ** 2
          du(k)   =
     #0.1D1 / dz ** 2 * dzdZ(i, j, k) ** 2 + 0.1D1 / dz * dzdZ2(i
     #, j, k) / 0.2D1
          rhs(k)  = 
     #-(V(i - 1, j, k) + V(i + 1, j, k)) * dxdX(i, j, k) ** 2 / dx
     # ** 2 - (-V(i - 1, j, k) + V(i + 1, j, k)) / dx * dxdX2(i, j, k) /
     # 0.2D1 - (V(i, j - 1, k) + V(i, j + 1, k)) * dydY(i, j, k) ** 2 / 
     #dy ** 2 - (-V(i, j - 1, k) + V(i, j + 1, k)) / dy * dydY2(i, j, k)
     # / 0.2D1 + defect(i, j, k)
     # + V_rhs(i,j,k)
        end do
        dl(Nz-2) = 0.0d0
        d(Nz-1)  = 1.0d0
        du(Nz-1) = 0.0d0
        rhs(Nz-1)  = V(i,j,Nz-1)
        dl(Nz-1) = 0.0d0
        d(Nz)    = 1.0d0
        rhs(Nz)  = V(i,j,Nz)
  
        nrhs = 1
        call dgtsv( Nz, nrhs, dl, d, du, rhs, Nz, info )
  
        do k=1,Nz
          V(i,j,k)=rhs(k)
        end do

      end do
      end do

c     Y-DIRECTION LINE RELAXATION
c     set up tridiagonal system; note that indexing on 
c     lower diagonal is always (j-1) when implementing the 
c     jth equation
      do i=3,Nx-2
      do k=3,Ny-2
        d(1)   = 1.0d0
        du(1)  = 0.0d0
        rhs(1) = V(i,1,k)
        dl(1)  = 0.0d0
        d(2)   = 1.0d0
        du(2)  = 0.0d0
        rhs(2) = V(i,2,k)
        do j=3,Ny-2  ! using 7-point stencil
          dl(j-1) =
     #0.1D1 / dy ** 2 * dydY(i, j, k) ** 2 - 0.1D1 / dy * dydY2(i, 
     #j, k) / 0.2D1
          d(j)    =
     #-2 / dx ** 2 * dxdX(i, j, k) ** 2 - 2 / dy ** 2 * dydY(i, j,
     # k) ** 2 - 2 / dz ** 2 * dzdZ(i, j, k) ** 2
          du(j)   =
     #0.1D1 / dy ** 2 * dydY(i, j, k) ** 2 + 0.1D1 / dy * dydY2(i,
     # j, k) / 0.2D1
          rhs(j)  = 
     #-0.1D1 / dx ** 2 * (V(i - 1, j, k) + V(i + 1, j, k)) * dxdX(
     #i, j, k) ** 2 - (-V(i - 1, j, k) + V(i + 1, j, k)) / dx * dxdX2(i,
     # j, k) / 0.2D1 - 0.1D1 / dz ** 2 * (V(i, j, k - 1) + V(i, j, k + 1
     #)) * dzdZ(i, j, k) ** 2 - (-V(i, j, k - 1) + V(i, j, k + 1)) / dz 
     #* dzdZ2(i, j, k) / 0.2D1 + defect(i, j, k)
     # + V_rhs(i,j,k)
        end do
        dl(Ny-2) = 0.0d0
        d(Ny-1)  = 1.0d0
        du(Ny-1) = 0.0d0
        rhs(Ny-1)  = V(i,Ny-1,k)
        dl(Ny-1) = 0.0d0
        d(Ny)    = 1.0d0
        rhs(Ny)  = V(i,Ny,k)
  
        nrhs = 1
        call dgtsv( Ny, nrhs, dl, d, du, rhs, Ny, info )
  
        do j=1,Ny
          V(i,j,k)=rhs(j)
        end do

      end do
      end do

c     X-DIRECTION LINE RELAXATION
c     set up tridiagonal system; note that indexing on 
c     lower diagonal is always (j-1) when implementing the 
c     jth equation
      do j=3,Ny-2
      do k=3,Nz-2
        d(1)   = 1.0d0
        du(1)  = 0.0d0
        rhs(1) = V(1,j,k)
        dl(1)  = 0.0d0
        d(2)   = 1.0d0
        du(2)  = 0.0d0
        rhs(2) = V(2,j,k)
        do i=3,Nx-2  ! using 7-point stencil
          dl(i-1) =
     #0.1D1 / dx ** 2 * dxdX(i, j, k) ** 2 - 0.1D1 / dx * dxdX2(i,
     #j, k) / 0.2D1
          d(i)    =
     #-2 / dx ** 2 * dxdX(i, j, k) ** 2 - 2 / dy ** 2 * dydY(i, j,
     # k) ** 2 - 2 / dz ** 2 * dzdZ(i, j, k) ** 2
          du(i)   =
     #0.1D1 / dx ** 2 * dxdX(i, j, k) ** 2 + 0.1D1 / dx * dxdX2(i,
     # j, k) / 0.2D1
          rhs(i)  = 
     #-0.1D1 / dy ** 2 * (V(i, j - 1, k) + V(i, j + 1, k)) * dydY(
     #i, j, k) ** 2 - (-V(i, j - 1, k) + V(i, j + 1, k)) / dy * dydY2(i,
     # j, k) / 0.2D1 - 0.1D1 / dz ** 2 * (V(i, j, k - 1) + V(i, j, k + 1
     #)) * dzdZ(i, j, k) ** 2 - (-V(i, j, k - 1) + V(i, j, k + 1)) / dz 
     #* dzdZ2(i, j, k) / 0.2D1 + defect(i, j, k)
     # + V_rhs(i,j,k)
        end do
        dl(Nx-2)  = 0.0d0
        d(Nx-1)   = 1.0d0
        du(Nx-1)  = 0.0d0
        rhs(Nx-1) = V(Nx-1,j,k)
        dl(Nx-1)  = 0.0d0
        d(Nx)     = 1.0d0
        rhs(Nx)   = V(Nx,j,k)
  
        nrhs = 1
        call dgtsv( Nx, nrhs, dl, d, du, rhs, Nx, info )
  
        do i=1,Nx
          V(i,j,k)=rhs(i)
        end do

      end do
      end do

c     compute residual
      do i=2,Nx-1
      do j=2,Ny-1
      do k=2,Nz-1
        if (cmask(i,j,k) .eq. CMASK_ON) then
          res =
     # dble(1 / dx ** 2 * (V(i - 1, j, k) - 2 * V(i, j, k) + V(i +
     # 1, j, k)) * dxdX(i, j, k) ** 2) + dble((-V(i - 1, j, k) + V(i + 1
     #, j, k)) / dx * dxdX2(i, j, k)) / 0.2D1 + dble(1 / dy ** 2 * (V(i,
     # j - 1, k) - 2 * V(i, j, k) + V(i, j + 1, k)) * dydY(i, j, k) ** 2
     #) + dble((-V(i, j - 1, k) + V(i, j + 1, k)) / dy * dydY2(i, j, k))
     # / 0.2D1 + dble(1 / dz ** 2 * (V(i, j, k - 1) - 2 * V(i, j, k) + V
     #(i, j, k + 1)) * dzdZ(i, j, k) ** 2) + dble((-V(i, j, k - 1) + V(i
     #, j, k + 1)) / dz * dzdZ2(i, j, k)) / 0.2D1 - defect(i, j, k)
     #-V_rhs(i,j,k)
          norm = norm + res**2
          sum  = sum + 1
        end if
      end do
      end do
      end do

      norm = sqrt(norm/sum)

      counter = counter + 1

c      write(*,*) counter, "relax:   ", Nx, Ny, Nz, norm

      return
      end subroutine

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     defect: compute the defect correction (see Sec. 5.4.1 of
c             Trottenberg et al., "Multigrid", Elsevier, 2000)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dc(V,defect,x,y,z,dx,dy,dz,Nx,Ny,Nz,Xx,Yy,Zz,
     &dxdX,dxdX2,dydY,dydY2,dzdZ,dzdZ2,e,phi1,phi2,pi1,pi2,modphi,At)
      use globals

      implicit none

      integer Nx, Ny, Nz
      real*8 dx, dy, dz
      real*8 V(Nx,Ny,Nz), defect(Nx,Ny,Nz)
      real*8 f(Nx,Ny,Nz)
      real*8 x(Nx), y(Ny), z(Nz)
      real*8 Xx(Nx,Ny,Nz), Yy(Nx,Ny,Nz), Zz(Nx,Ny,Nz)
      real*8 dxdX(Nx,Ny,Nz), dxdX2(Nx,Ny,Nz)
      real*8 dydY(Nx,Ny,Nz), dydY2(Nx,Ny,Nz)
      real*8 dzdZ(Nx,Ny,Nz), dzdZ2(Nx,Ny,Nz)

      real*8 e
      real*8 phi1(Nx,Ny,Nz), phi2(Nx,Ny,Nz)
      real*8 pi1(Nx,Ny,Nz), pi2(Nx,Ny,Nz)
      real*8 modphi(Nx,Ny,Nz), At(Nx,Ny,Nz)

      real*8 d4(Nx,Ny,Nz)
      integer i, j, k, l

      ! source term
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
         f(i,j,k)=2d0*e*(phi2(i,j,k)*pi1(i,j,k)-phi1(i,j,k)*pi2(i,j,k))
     &          +2d0*e**2d0*At(i,j,k)*modphi(i,j,k)**2d0
      end do
      end do
      end do

      ! construct defect term using fourth-order differences
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
       d4(i,j,k) = f(i,j,k)
      end do
      end do
      end do

      do i=3,Nx-2
      do j=3,Ny-2
      do k=3,Nz-2
       d4(i,j,k) =
     #f(i, j, k) - dble((-V(i - 2, j, k) + 16 * V(i - 1, j, k) - 30
     # * V(i, j, k) + 16 * V(i + 1, j, k) - V(i + 2, j, k)) / dx ** 2 * 
     #dxdX(i, j, k) ** 2) / 0.12D2 - dble((V(i - 2, j, k) - 8 * V(i - 1,
     # j, k) + 8 * V(i + 1, j, k) - V(i + 2, j, k)) / dx * dxdX2(i, j, k
     #)) / 0.12D2 - dble((-V(i, j - 2, k) + 16 * V(i, j - 1, k) - 30 * V
     #(i, j, k) + 16 * V(i, j + 1, k) - V(i, j + 2, k)) / dy ** 2 * dydY
     #(i, j, k) ** 2) / 0.12D2 - dble((V(i, j - 2, k) - 8 * V(i, j - 1, 
     #k) + 8 * V(i, j + 1, k) - V(i, j + 2, k)) / dy * dydY2(i, j, k)) /
     # 0.12D2 - dble((-V(i, j, k - 2) + 16 * V(i, j, k - 1) - 30 * V(i, 
     #j, k) + 16 * V(i, j, k + 1) - V(i, j, k + 2)) / dz ** 2 * dzdZ(i, 
     #j, k) ** 2) / 0.12D2 - dble((V(i, j, k - 2) - 8 * V(i, j, k - 1) +
     # 8 * V(i, j, k + 1) - V(i, j, k + 2)) / dz * dzdZ2(i, j, k)) / 0.1
     #2D2
      end do
      end do
      end do

      ! construct the new RHS
      do i=1,Nx
      do j=1,Ny
      do k=1,Nz
       defect(i,j,k) = f(i,j,k)
      end do
      end do
      end do

      do i=2,Nx-1
      do j=2,Ny-1
      do k=2,Nz-1
       defect(i,j,k) = d4(i,j,k) +
     #dble(1 / dx ** 2 * (V(i - 1, j, k) - 2 * V(i, j, k) + V(i + 1
     #, j, k)) * dxdX(i, j, k) ** 2) + dble((-V(i - 1, j, k) + V(i + 1,
     #j, k)) / dx * dxdX2(i, j, k)) / 0.2D1 + dble(1 / dy ** 2 * (V(i, j
     # - 1, k) - 2 * V(i, j, k) + V(i, j + 1, k)) * dydY(i, j, k) ** 2)
     #+ dble((-V(i, j - 1, k) + V(i, j + 1, k)) / dy * dydY2(i, j, k)) /
     # 0.2D1 + dble(1 / dz ** 2 * (V(i, j, k - 1) - 2 * V(i, j, k) + V(i
     #, j, k + 1)) * dzdZ(i, j, k) ** 2) + dble((-V(i, j, k - 1) + V(i,
     #j, k + 1)) / dz * dzdZ2(i, j, k)) / 0.2D1
      end do
      end do
      end do

      return
      end subroutine
