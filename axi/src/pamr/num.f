      module globals
      implicit none

      integer counter
      integer N1, N2
      real*8 Q, sigma, pi
      parameter ( Q = 1.0d0,   sigma = 1.0d0 )
      parameter ( pi = 3.14159265358979d0 )
c      parameter ( N1 = 1025, N2 = 2049 )
      parameter ( N1 = 2049, N2 = 4097 )

      end module globals

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     initguess: initializes a guess for the multigrid solve
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine initguess(V,R,X,NR,NX,pR,zX,A_t)
      use globals

      implicit none

      integer NR,NX
      real*8 V(NR,NX), R(NR), X(NX), pR(NR,NX), zX(NR,NX)
      real*8 A_t(NR,NX)
      real*8 rand

      integer i, j

      counter = 0

      ! Initial guess
      do i=1,NR
      do j=1,NX
c        call RANDOM_NUMBER(rand)
c        V(i,j)=rand/4d0
c        V(i,j)=0.25d0
        V(i,j)=A_t(i,j)
      end do
      end do

      return
      end subroutine initguess

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     init_qball2: initializes gauged Q-ball fields using data from the
c                  multigrid solve; reads in data using rff() function
c                  and interpolates using bicubic() function
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_qball2(V_n,R,X,NR,NX,pR_n,zX_n,dR,dX,
     &dRdp_n,dRdp2_n,dXdz_n,dXdz2_n,A_t_n,Atilde_p_n,A_z_n,
     &B_t_n,Btilde_p_n,B_z_n,x1min,x1max,x2min,x2max)
      use globals

      implicit none

      integer NR,NX
      real*8 dR,dX
      real*8 V_n(NR,NX), R(NR), X(NX), pR_n(NR,NX), zX_n(NR,NX)
      real*8 dRdp_n(NR,NX),dRdp2_n(NR,NX),dXdz_n(NR,NX),dXdz2_n(NR,NX)
      real*8 A_t_n(NR,NX), Atilde_p_n(NR,NX), A_z_n(NR,NX)
      real*8 B_t_n(NR,NX), Btilde_p_n(NR,NX), B_z_n(NR,NX)

      real*8 res
      integer i, j

      real*8  u(N1,N2)
      real*8  dx1, dx2
      real*8  x1, x1l, x1u, x2, x2l, x2u
      integer x1lN, x2lN, x1uN, x2uN
      real*8  x1min,x1max,x2min,x2max
      real*8, dimension(4) :: y, y1, y2, y12
      real*8  ansy,ansy1,ansy2

      real*8 f_x1, f_x2, f_x12
      real*8 f_x1_f, f_x2_f, f_x12_ff
      real*8 f_x1_b, f_x2_b, f_x12_bb
      real*8 f_x12_fb, f_x12_bf

      dx1 = (x1max-x1min)/dble(N1-1)
      dx2 = (x2max-x2min)/dble(N2-1)

      call rff(u,N1,N2)  ! N1, N2 needs to match actual size used in V.bin

      do i=1,NR
      do j=1,NX

        x1 = R(i)
        x2 = X(j)
        x1l = dx1*floor(x1/dx1)
        x1u = x1l+dx1
        x2l = dx2*floor(x2/dx2)
        x2u = x2l+dx2

        x1lN = floor(abs(x1min-x1)/dx1)+1
        x1uN = x1lN+1
        x2lN = floor(abs(x2min-x2)/dx2)+1
        x2uN = x2lN+1

        ! hacky fix for when interpolating rectangle extends beyond
        ! the outer boundary when x=xmax
        if (x1 .GE. x1max) then
          x1l = x1l-dx1
          x1u = x1u-dx1
          x1lN = x1lN-1
          x1uN = x1uN-1
        end if
        if (x2 .GE. x2max) then
          x2l = x2l-dx2
          x2u = x2u-dx2
          x2lN = x2lN-1
          x2uN = x2uN-1
        end if

        y(1)   = u(x1lN,x2lN)
        y(2)   = u(x1uN,x2lN)
        y(3)   = u(x1uN,x2uN)
        y(4)   = u(x1lN,x2uN)

        if ((x1lN .LE. 2) .AND. (x2lN .LE. 2)) then
          y1(1)  = f_x1_f(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_f(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_f(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_f(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_f(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_f(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_f(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_f(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_ff(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_ff(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_ff(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_ff(u,x1lN,x2uN,dx1,dx2)
        else if ((x1lN .LE. 2) .AND. (x2uN .GE. (N2-2))) then
          y1(1)  = f_x1_f(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_f(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_f(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_f(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_b(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_b(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_b(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_b(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_fb(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_fb(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_fb(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_fb(u,x1lN,x2uN,dx1,dx2)
        else if ((x1uN .GE. (N1-2)) .AND. (x2uN .LE. 2)) then
          y1(1)  = f_x1_b(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_b(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_b(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_b(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_f(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_f(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_f(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_f(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_bf(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_bf(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_bf(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_bf(u,x1lN,x2uN,dx1,dx2)
        else if ((x1uN .GE. (N1-2)) .AND. (x2uN .GE. (N2-2))) then
          y1(1)  = f_x1_b(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_b(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_b(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_b(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_b(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_b(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_b(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_b(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_bb(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_bb(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_bb(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_bb(u,x1lN,x2uN,dx1,dx2)
        else if (x1lN .LE. 2) then
          y1(1)  = f_x1_f(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_f(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_f(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_f(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_ff(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_ff(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_ff(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_ff(u,x1lN,x2uN,dx1,dx2)
        else if (x1uN .GE. (N1-2)) then
          y1(1)  = f_x1_b(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1_b(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1_b(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1_b(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_bb(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_bb(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_bb(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_bb(u,x1lN,x2uN,dx1,dx2)
        else if (x2lN .LE. 2) then
          y1(1)  = f_x1(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_f(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_f(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_f(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_f(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_ff(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_ff(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_ff(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_ff(u,x1lN,x2uN,dx1,dx2)
        else if (x2uN .GE. (N2-2)) then
          y1(1)  = f_x1(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2_b(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2_b(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2_b(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2_b(u,x1lN,x2uN,dx2)
          y12(1) = f_x12_bb(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12_bb(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12_bb(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12_bb(u,x1lN,x2uN,dx1,dx2)
        else
          y1(1)  = f_x1(u,x1lN,x2lN,dx1)
          y1(2)  = f_x1(u,x1uN,x2lN,dx1)
          y1(3)  = f_x1(u,x1uN,x2uN,dx1)
          y1(4)  = f_x1(u,x1lN,x2uN,dx1)
          y2(1)  = f_x2(u,x1lN,x2lN,dx2)
          y2(2)  = f_x2(u,x1uN,x2lN,dx2)
          y2(3)  = f_x2(u,x1uN,x2uN,dx2)
          y2(4)  = f_x2(u,x1lN,x2uN,dx2)
          y12(1) = f_x12(u,x1lN,x2lN,dx1,dx2)
          y12(2) = f_x12(u,x1uN,x2lN,dx1,dx2)
          y12(3) = f_x12(u,x1uN,x2uN,dx1,dx2)
          y12(4) = f_x12(u,x1lN,x2uN,dx1,dx2)
        end if

        call bicubic(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
        V_n(i,j) = ansy

      end do
      end do

      ! compute relevant quantities using O(h^2) differencing
      do i=2,NR-1
      do j=2,NX-1

        B_t_n(i,j) =
     # dble(2 * Atilde_p_n(i, j)) + pR_n(i, j) * (Atilde_p_n(i + 1,
     # j) - Atilde_p_n(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 + (A_z_n(i
     #, j + 1) - A_z_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1

        Btilde_p_n(i,j) =
     # ((A_t_n(i + 1, j) - A_t_n(i - 1, j)) / dR * dRdp_n(i, j) / 0
     #.2D1 - (V_n(i + 1, j) - V_n(i - 1, j)) / dR * dRdp_n(i, j)/0.2D1)/ 
     #pR_n(i, j)

        B_z_n(i,j) =
     # (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (V_n(i, j + 1) - V_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1

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
     #_n(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2*V_n(i,
     # j)-5*V_n(i + 1, j) + 4 * V_n(i + 2, j) - V_n(i + 3, j)) / dR** 2*
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(V_n(i, j)) +dble(2 *V
     #_n(i+1, j)) - dble(V_n(i + 2, j)) / 0.2D1) / dble(dR) *dRdp2_n(i,j
     #)
      Btilde_p_n(i,j)=res

      res = (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (V_n(i, j + 1) - V_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
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

      res = ((A_t_n(i + 1, j) - A_t_n(i - 1, j)) / dR * dRdp_n(i,j)/ 0
     #.2D1 - (V_n(i + 1, j) - V_n(i - 1, j)) / dR * dRdp_n(i, j)/0.2D1)/
     #pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1)) + A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1 / 0.2D1*V_n(i
     #, j) - dble(2 * V_n(i, j - 1)) + V_n(i, j - 2) / 0.2D1)/dX*dXdz_n(
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
     #A_t_n(i - 2, j)/0.2D1)/ dR * dRdp_n(i, j) - (0.3D1 / 0.2D1 * V_n(
     #i, j) - dble(2 * V_n(i - 1, j)) + V_n(i - 2, j)/0.2D1)/dR*dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (A_t_n(i, j + 1) - A_t_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.
     #2D1 - (V_n(i, j + 1) - V_n(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1
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
     #.2D1 - (V_n(i + 1, j) - V_n(i - 1, j)) / dR * dRdp_n(i, j)/0.2D1)/
     #pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1)) - 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1 * V
     #_n(i, j)+dble(2 * V_n(i, j + 1)) - V_n(i, j + 2) / 0.2D1)/dX*dXdz_
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
     #/ 0.2D1 * dble(A_t_n(i, j)) + dble(2 * A_t_n(i + 1, j)) -dble(A_t
     #_n(i + 2,j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2 *V_n(i,
     # j)-5 * V_n(i + 1, j) + 4 * V_n(i + 2, j) - V_n(i + 3, j))/dR**2*
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(V_n(i, j)) +dble(2*V
     #_n(i+1,j)) - dble(V_n(i + 2, j)) / 0.2D1)/dble(dR)*dRdp2_n(i,j
     #)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1))- 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1*V
     #_n(i,j)+ dble(2 * V_n(i, j + 1)) - V_n(i, j + 2)/0.2D1)/dX*dXdz_
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
     #_n(i + 2, j)) / 0.2D1) / dble(dR) * dRdp2_n(i, j) - dble((2*V_n(i,
     # j) - 5 * V_n(i + 1, j)+4*V_n(i+2,j)-V_n(i + 3, j)) / dR ** 2 *
     # dRdp_n(i, j) ** 2) - (-0.3D1 / 0.2D1 * dble(V_n(i, j))+ dble(2*V
     #_n(i+1, j)) - dble(V_n(i + 2, j)) / 0.2D1) / dble(dR) *dRdp2_n(i,j
     #)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1)) + A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1/0.2D1*V_n(i
     #, j) - dble(2 * V_n(i, j - 1)) + V_n(i, j - 2) / 0.2D1)/dX*dXdz_n(
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
     #A_t_n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) - (0.3D1/ 0.2D1*V_n(
     #i, j) - dble(2 * V_n(i - 1, j)) + V_n(i - 2, j) /0.2D1)/dR*dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (0.3D1 / 0.2D1 * A_t_n(i, j) - dble(2 * A_t_n(i, j - 1))+A
     #_t_n(i, j - 2) / 0.2D1) / dX * dXdz_n(i, j) - (0.3D1/0.2D1*V_n(i
     #, j) - dble(2 * V_n(i, j - 1)) + V_n(i,j-2)/0.2D1)/dX*dXdz_n(
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
     #A_t_n(i - 2, j) / 0.2D1) / dR * dRdp_n(i, j) - (0.3D1/0.2D1*V_n(
     #i, j) - dble(2 * V_n(i - 1, j)) + V_n(i - 2, j) / 0.2D1)/dR*dRdp_n
     #(i, j)) / pR_n(i, j)
      Btilde_p_n(i,j)=res

      res = (-0.3D1 / 0.2D1 * A_t_n(i, j) + dble(2 * A_t_n(i, j + 1)) - 
     #A_t_n(i, j + 2) / 0.2D1) / dX * dXdz_n(i, j) - (-0.3D1 / 0.2D1 *V
     #_n(i, j)+dble(2 * V_n(i, j + 1)) - V_n(i, j + 2) / 0.2D1)/dX*dXdz_
     #n(i, j)
      B_z_n(i,j)=res

      return
      end subroutine init_qball2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     List of derivatives
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! centered derivative in x1 direction
      function f_x1(f,i,j,h)
      use globals
      real*8 f_x1
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x1 = (f(i+1,j)-f(i-1,j))/(2*h)
      return
      end function f_x1

      ! centered derivative in x2 direction
      function f_x2(f,i,j,h)
      use globals
      real*8 f_x2
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x2 = (f(i,j+1)-f(i,j-1))/(2*h)
      return
      end function f_x2

      ! centered mixed derivative
      function f_x12(f,i,j,h1,h2)
      use globals
      real*8 f_x12
      real*8 f(N1,N2)
      integer i, j
      real*8 h1, h2
      f_x12 = (f(i+1, j+1)-f(i+1, j-1)-f(i-1,
     &         j+1)+f(i-1, j-1))/(4*h1*h2)
      return
      end function f_x12

      ! forward derivative in x1 direction
      function f_x1_f(f,i,j,h)
      use globals
      real*8 f_x1_f
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x1_f = ((-3d0)*f(i,j)+4d0*f(i+1,j)-f(i+1,j))/(2*h)
      return
      end function f_x1_f

      ! forward derivative in x2 direction
      function f_x2_f(f,i,j,h)
      use globals
      real*8 f_x2_f
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x2_f = ((-3d0)*f(i,j)+4d0*f(i,j+1)-f(i,j+1))/(2*h)
      return
      end function f_x2_f

      ! forward mixed derivative
      function f_x12_ff(f,i,j,h1,h2)
      use globals
      real*8 f_x12_ff
      real*8 f(N1,N2)
      integer i, j
      real*8 h1, h2
      f_x12_ff = (9*f(i, j) - 12*f(i, j + 1) + 3*f(i, j + 2) -
     & 12*f(i + 1, j) + 16*f(i + 1, j + 1) - 4*f(i + 1, j + 2) +
     & 3*f(i + 2, j) - 4*f(i + 2, j + 1) + f(i + 2, j + 2))/(4*h2*h1)
      return
      end function f_x12_ff

      ! backward derivative in x1 direction
      function f_x1_b(f,i,j,h)
      use globals
      real*8 f_x1_b
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x1_b = (3*f(i, j) - 4*f(i - 1, j) + f(i - 2, j))/(2*h)
      return
      end function f_x1_b

      ! backward derivative in x2 direction
      function f_x2_b(f,i,j,h)
      use globals
      real*8 f_x2_b
      real*8 f(N1,N2)
      integer i, j
      real*8 h
      f_x2_b = (3*f(i, j) - 4*f(i, j-1) + f(i, j-2))/(2*h)
      return
      end function f_x2_b

      ! backward mixed derivative
      function f_x12_bb(f,i,j,h1,h2)
      use globals
      real*8 f_x12_bb
      real*8 f(N1,N2)
      integer i, j
      real*8 h1, h2
      f_x12_bb = (9*f(i, j) - 12*f(i, j - 1) + 3*f(i, j - 2) -
     &12*f(i - 1, j) + 16*f(i - 1, j - 1) - 4*f(i - 1, j - 2) +
     &3*f(i - 2, j) - 4*f(i - 2, j - 1) + f(i - 2, j - 2))/(4*h2*h1)
      return
      end function f_x12_bb

      ! mixed derivative: backward x1, forward x2
      function f_x12_bf(f,i,j,h1,h2)
      use globals
      real*8 f_x12_bf
      real*8 f(N1,N2)
      integer i, j
      real*8 h1, h2
      f_x12_bf = (-9*f(i, j) + 12*f(i, j + 1) - 3*f(i, j + 2) +
     &12*f(i - 1, j) - 16*f(i - 1, j + 1) + 4*f(i - 1, j + 2) -
     &3*f(i - 2, j) + 4*f(i - 2, j + 1) - f(i - 2, j + 2))/(4*h2*h1)
      return
      end function f_x12_bf

      ! mixed derivative: forward x1, backward x2
      function f_x12_fb(f,i,j,h1,h2)
      use globals
      real*8 f_x12_fb
      real*8 f(N1,N2)
      integer i, j
      real*8 h1, h2
      f_x12_fb = (-9*f(i, j) + 12*f(i, j - 1) - 3*f(i, j - 2) +
     &12*f(i + 1, j) - 16*f(i + 1, j - 1) + 4*f(i + 1, j - 2) -
     &3*f(i + 2, j) + 4*f(i + 2, j - 1) - f(i + 2, j - 2))/(4*h2*h1)
      return
      end function f_x12_fb

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     lop: differential operator L(V)=0 computed where cmask=CMASK_ON
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lop(LV,V,cmask,R,X,NR,NX,dR,dX,pR,zX,dRdp,dRdp2,dXdz,
     &dXdz2,e,phi1,phi2,pi1,pi2,modphi,A_t)
      use globals

      implicit none

      integer NR, NX
      real*8 V(NR,NX)
      real*8 cmask(NR,NX), LV(NR,NX)
      real*8 R(NR), X(NX)
      real*8 dR, dX
      real*8 pR(NR,NX), zX(NR,NX)
      real*8 dRdp(NR,NX), dRdp2(NR,NX)
      real*8 dXdz(NR,NX), dXdz2(NR,NX)

      integer i,j
      real*8 f(NR,NX)

      real*8 e
      real*8 phi1(NR,NX), phi2(NR,NX), pi1(NR,NX), pi2(NR,NX)
      real*8 modphi(NR,NX), A_t(NR,NX)

      include 'cmask.inc'

      ! source term
      do i=1,NR
      do j=1,NX
c        f(i,j)=Q/(sigma*sqrt(2d0*pi))**3d0
c     &         *exp(-sqrt(pR(i,j)**2d0+zX(i,j)**2d0)**2d0
c     &         /(2d0*sigma**2d0))
        f(i,j)=2d0*e*(phi2(i,j)*pi1(i,j)-phi1(i,j)*pi2(i,j))
     &         +2d0*e**2d0*A_t(i,j)*modphi(i,j)**2d0
      end do
      end do

      ! interior points
      do i=2,NR-1
      do j=2,NX-1
        if (cmask(i,j) .eq. CMASK_ON) then
          LV(i,j)=
     &dble((V(i - 1,j) - 2 * V(i,j) + V(i + 1,j)) / dR ** 2 * dRdp(i,
     &j) ** 2) + dble((-V(i - 1,j) + V(i + 1,j)) / dR * dRdp2(i, j))
     &/ 0.2D1 + dble((-V(i - 1,j) + V(i + 1,j)) / dR * dRdp(i, j) /
     &pR(i, j)) / 0.2D1 + dble((V(i,j - 1) - 2 * V(i,j) + V(i,j + 1))
     &/ dX ** 2 * dXdz(i, j) ** 2) + dble((-V(i,j - 1) + V(i,j + 1))
     &/ dX * dXdz2(i, j)) / 0.2D1 - f(i, j)
        end if
      end do
      end do

      ! Rmin boundary
      i=1
      do j=1,NX
        if (cmask(i,j).eq.CMASK_ON) then
          LV(i,j)=
     #(-0.3D1 / 0.2D1 * V(i,j) + 0.2D1 * V(i + 1,j) - V(i + 2,j) / 0.2D1
     #) / dR * dRdp(i,j)
        end if
      end do

      ! Rmax boundary
      i=NR
      do j=1,NX
        if (cmask(i,j).eq.CMASK_ON) then
          LV(i,j)=
     #(0.3D1 / 0.2D1 * V(i,j) - 0.2D1 * V(i - 1,j) + V(i - 2,j) / 0.2D1)
     #/ dR * dRdp(i,j) + pR(i,j) * V(i,j) / (pR(i,j) ** 2 + zX(i,j) **2)
        end if
      end do

      ! Xmin boundary
      j=1
      do i=1,NR
        if (cmask(i,j).eq.CMASK_ON) then
          LV(i,j)=
     #(-0.3D1 / 0.2D1 * V(i,j) + 0.2D1 * V(i,j + 1) - V(i,j + 2) / 0.2D1
     #) / dX * dXdz(i,j) + zX(i,j) * V(i,j) / (pR(i,j) ** 2 +zX(i,j)**2)
        end if
      end do

      ! Xmax boundary
      j=NX
      do i=1,NR
        if (cmask(i,j).eq.CMASK_ON) then
          LV(i,j)=
     #(0.3D1 / 0.2D1 * V(i,j) - 0.2D1 * V(i,j - 1) + V(i,j - 2) / 0.2D1)
     #/ dX * dXdz(i,j) + zX(i,j) * V(i,j) / (pR(i,j) ** 2 + zX(i,j) **2)
        end if
      end do

c      write(*,*) counter, "lop:     ", NR, NX

      return
      end subroutine lop

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     residual: residual L[V]-rhs computed where cmask=CMASK_ON
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine residual(res,rhs,V,cmask,R,X,norm,NR,NX,dR,dX,pR,zX,
     &dRdp,dRdp2,dXdz,dXdz2,e,phi1,phi2,pi1,pi2,modphi,A_t)
      use globals

      implicit none

      integer NR, NX
      real*8 V(NR,NX)
      real*8 cmask(NR,NX), res(NR,NX), rhs(NR,NX)
      real*8 R(NR), X(NX), norm
      real*8 dR, dX
      real*8 pR(NR,NX), zX(NR,NX)
      real*8 dRdp(NR,NX), dRdp2(NR,NX)
      real*8 dXdz(NR,NX), dXdz2(NR,NX)
      integer i, j, sum

      real*8 e
      real*8 phi1(NR,NX), phi2(NR,NX), pi1(NR,NX), pi2(NR,NX)
      real*8 modphi(NR,NX), A_t(NR,NX)

      include 'cmask.inc'

      call lop(res,V,cmask,R,X,NR,NX,dR,dX,pR,zX,dRdp,dRdp2,dXdz,dXdz2,
     &e,phi1,phi2,pi1,pi2,modphi,A_t)

      norm = 0d0
      sum  = 0

      do i=1,NR
      do j=1,NX
        if (cmask(i,j).eq.CMASK_ON) then
          res(i,j) = res(i,j)-rhs(i,j)
          norm     = norm+res(i,j)**2
          sum      = sum+1
        end if
      end do
      end do

      norm = sqrt(norm/sum)

c      write(*,*) counter, "residual:", NR, NX, norm

      return
      end subroutine residual

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     relax: applies alternating-direction line relaxation using LAPACK
c            'dgtsv' solver
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine relax(V,V_rhs,cmask,phys_bdy,norm,R,X,NR,NX,dR,dX,
     &pR,zX,dRdp,dRdp2,dXdz,dXdz2,e,phi1,phi2,pi1,pi2,modphi,A_t)
      use globals

      implicit none

      integer NR, NX
      real*8 V(NR,NX), V_rhs(NR,NX), cmask(NR,NX)
      integer phys_bdy(4)
      real*8 norm
      real*8 dR, dX
      real*8 R(NR), X(NX)
      real*8 pR(NR,NX), zX(NR,NX)
      real*8 dRdp(NR,NX), dRdp2(NR,NX)
      real*8 dXdz(NR,NX), dXdz2(NR,NX)

      integer i, j, k, l
      real*8, dimension (NX) :: d, du, dl, rhs
      integer nrhs, info
      real*8, dimension (NR, NX) :: f
      real*8 res
      integer sum

      real*8 e
      real*8 phi1(NR,NX), phi2(NR,NX), pi1(NR,NX), pi2(NR,NX)
      real*8 modphi(NR,NX), A_t(NR,NX)

      include 'cmask.inc'

      norm = 0d0
      sum  = 0

      ! source term
      do i=1,NR
      do j=1,NX
        f(i,j)=2d0*e*(phi2(i,j)*pi1(i,j)-phi1(i,j)*pi2(i,j))
     &         +2d0*e**2d0*A_t(i,j)*modphi(i,j)**2d0
      end do
      end do

c     X-DIRECTION LINE RELAXATION
c     set up tridiagonal system; note that indexing on
c     lower diagonal is always (j-1) when implementing the
c     jth equation
c      do l=0,1  ! zebra sweep
c      do i=2+l,NR-1-l,2
      do i=2,NR-1
        d(1)   = 1.0d0
        du(1)  = 0.0d0
        rhs(1) = V(i,1)
        do j=2,NX-1  ! using 5-point stencil
          dl(j-1) =
     #0.1D1 / dX ** 2 * dXdz(i, j) ** 2 - 0.1D1 / dX * dXdz2(i,
     # j) / 0.2D1
          d(j) =
     #-2 / dR ** 2 * dRdp(i, j) ** 2 - 2 / dX ** 2 * dXdz(i, j
     #) ** 2
          du(j) =
     #0.1D1 / dX ** 2 * dXdz(i, j) ** 2 + 0.1D1 / dX * dXdz2(i
     #, j) / 0.2D1
          rhs(j) =
     #f(i, j) - (V(i - 1,int(j)) + V(i + 1,int(j))) / dR**
     # 2 * dRdp(i, j) ** 2 - (-V(i-1,int(j))+V(i+1,int(j))) /
     #dR * dRdp2(i, j) / 0.2D1-(-V(i-1,int(j))+V(i+1,int(j))
     #) / dR * dRdp(i, j) / pR(i, j) / 0.2D1 + V_rhs(i,j)
        end do
        dl(NX-1) = 0.0d0
        d(NX)    = 1.0d0
        rhs(NX)  = V(i,NX)

        nrhs = 1
        call dgtsv( NX, nrhs, dl, d, du, rhs, NX, info )

        do j=1,NX
          V(i,j)=rhs(j)
        end do

      end do

c     Y-DIRECTION LINE RELAXATION
c     set up tridiagonal system; note that indexing on
c     lower diagonal is always (j-1) when implementing the
c     jth equation
      do j=2,NX-1
        d(1)   = 1.0d0
        du(1)  = 0.0d0
        rhs(1) = V(1,j)
        do i=2,NR-1  ! using 5-point stencil
          dl(i-1) =
     #0.1D1 / dR ** 2 * dRdp(i, j) ** 2 - 0.1D1 / dR * dRdp2(i
     #, j) / 0.2D1 - 0.1D1 / dR * dRdp(i, j) / pR(i, j) / 0.2D1
          d(i) =
     # -2 / dR ** 2 * dRdp(i, j) ** 2 - 2 / dX ** 2 * dXdz(i, j
     #) ** 2
          du(i) =
     #0.1D1 / dR ** 2 * dRdp(i, j) ** 2 + 0.1D1 / dR * dRdp2(i
     #, j) / 0.2D1 + 0.1D1 / dR * dRdp(i, j) / pR(i, j) / 0.2D1
          rhs(i) =
     #f(i, j) - (V(int(i),j - 1) + V(int(i),j + 1)) / dX *
     #* 2 * dXdz(i, j)**2 - (-V(int(i),j-1) + V(int(i),j+1))/
     # dX * dXdz2(i, j) / 0.2D1 + V_rhs(i,j)
        end do
        dl(NR-1) = 0.0d0
        d(NR)    = 1.0d0
        rhs(NR)  = V(NR,j)

        nrhs = 1
        call dgtsv( NR, nrhs, dl, d, du, rhs, NR, info )

        do i=1,NR
          V(i,j)=rhs(i)
        end do

      end do

      ! Rmin boundary condition
      i=1
      do j=1,NX
        V(i,j) =
     #-dble((2 * dR * V_rhs(i,j) - 4 * V(i + 1,j) * dRdp(i,j) + V(i + 2,
     #j) * dRdp(i,j)) / dRdp(i,j)) / 0.3D1
      end do

      ! Rmax boundary condition
      i=NR
      do j=1,NX
        V(i,j) =
     #(2 * dR * V_rhs(i,j) - V(i - 2,j) * dRdp(i,j) + 4 * V(i - 1,j) *
     #dRdp(i,j)) / (3 * dRdp(i,j) * pR(i,j) ** 2 + 3 * dRdp(i,j) * zX(i,
     #j) ** 2 + 2 * dR * pR(i,j)) * (pR(i,j) ** 2 + zX(i,j) ** 2)
      end do

      ! Xmin boundary condition
      j=1
      do i=1,NR
        V(i,j) =
     #(2 * dX * V_rhs(i,j) - 4 * V(i,j + 1) * dXdz(i,j) + V(i,j + 2) *
     #dXdz(i,j)) / (-3 * dXdz(i,j) * pR(i,j) ** 2 - 3 * dXdz(i,j) *
     #zX(i,j) ** 2 + 2 * dX * zX(i,j)) * (pR(i,j) ** 2 + zX(i,j) ** 2)
      end do

      ! Xmax boundary condition
      j=NX
      do i=1,NR
        V(i,j) =
     #(2 * dX * V_rhs(i,j) - V(i,j - 2) * dXdz(i,j) + 4 * V(i,j - 1) *
     #dXdz(i,j)) / (3 * dXdz(i,j) * pR(i,j) ** 2 + 3 * dXdz(i,j) *
     #zX(i,j) ** 2 + 2 * dX * zX(i,j)) * (pR(i,j) ** 2 + zX(i,j) ** 2)
      end do

      ! compute residual
      do i=2,NR-1
      do j=2,NX-1
        if (cmask(i,j).eq.CMASK_ON) then
          res =
     &dble((V(i - 1,j) - 2 * V(i,j) + V(i + 1,j)) / dR ** 2 * dRdp(i,
     &j) ** 2) + dble((-V(i - 1,j) + V(i + 1,j)) / dR * dRdp2(i, j))
     &/ 0.2D1 + dble((-V(i - 1,j) + V(i + 1,j)) / dR * dRdp(i, j) /
     &pR(i, j)) / 0.2D1 + dble((V(i,j - 1) - 2 * V(i,j) + V(i,j + 1))
     &/ dX ** 2 * dXdz(i, j) ** 2) + dble((-V(i,j - 1) + V(i,j + 1))
     &/ dX * dXdz2(i, j)) / 0.2D1 - f(i, j) - V_rhs(i,j)
          norm = norm + res**2
          sum  = sum + 1
        end if
      end do
      end do

      norm = sqrt(norm/sum)

      counter = counter + 1

c      write(*,*) counter, "relax:   ", NR, NX, norm

      return
      end subroutine relax
