!----------------------------------------------------------------------
!  A Fortran function designed to be called by a PAMR hook function.
!  Integrates a gridfunction over the entire domain using an O(h^2)
!  trapezoidal rule.
!----------------------------------------------------------------------
      real*8 function qtotcalcpamr(qdensity_np1,qdensity_n,g1_Nx,g1_Ny,
     &                             dx,dy,gw) 
        implicit none

        include 'globals.inc'

        integer g1_Nx,g1_Ny
        real*8  qdensity_np1(1:g1_Nx,1:g1_Ny)
        real*8  qdensity_n(1:g1_Nx,1:g1_Ny)
        real*8  dx  
        real*8  dy  
        integer gw(4)

c=======================================================================
c     Code to manually update (integrate) the total charge Q
c=======================================================================

        integer i, j
        real*8  qtotsum    

        qtotsum = 0.0d0

        do i=1+gw(1), g1_Nx-1-gw(2)
        do j=1+gw(3), g1_Ny-1-gw(4)
          qtotsum=qtotsum+dx*dy/4.0d0*
     &                    (  qdensity_n(i,j)+qdensity_n(i+1,j)
     &                     + qdensity_n(i,j+1)+qdensity_n(i+1,j+1) )
        end do
        end do

        ! The function returns the total charge integrated over
        ! the entire domain
        qtotcalcpamr = qtotsum
        return
      end 
