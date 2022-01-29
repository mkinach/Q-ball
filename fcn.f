c===========================================================
c     Driver routine which integrates ODEs defining global
c     Q-ball shape.
c===========================================================
      subroutine fcn(neq,r,y,yprime)
         implicit    none

         include    'fcn.inc'

         integer     neq
         real*8      r,     y(neq),    yprime(neq)

         real*8      s,     u

         s = y(1)  ! initial conditions defined in init_qball.inc
         u = y(2)

         yprime(1) = u

         if ( r .eq. 0.0d0 ) then ! L'Hopital's rule for r = 0
            yprime(2) = 0.5d0*(-(w**2.0d0-1.0d0)*s
     &                -A*s**2.0d0+B*s**3.0d0)
         else  
            yprime(2) = -(u/r)-(w**2.0d0-1.0d0)*s
     &                -A*s**2.0d0+B*s**3.0d0
         endif
 
         return
      end 

c===========================================================
c     Dummy Jacobian routine.
c===========================================================
      subroutine jac
         implicit    none

         include    'fcn.inc'

         return
      end

c===========================================================
c     Interpolation routine and dependencies from lvutil
c     and dveclib.f in rnpletal
c===========================================================
 
c-----------------------------------------------------------------------
c
c     (NINTRP-1)th degree interpolation from one set of values to another
c     (Both sets must be in ascending order.)
c
c-----------------------------------------------------------------------
c
      subroutine dvinqn(v,x,vbar,xbar,n,nbar,vs,vf,nintrp)
c
         implicit     logical (a-z)
c
         real*8       flipn
         integer      llogic
         logical      leven
c
         integer      n, nbar, nintrp
         real*8       v(n), vbar(nbar), x(n), xbar(nbar)
c
         real*8       vs, vf
         integer      i, j, jbar

         logical      ltrace
         parameter  ( ltrace = .false. )
c
         if( ltrace ) then
            WRITE(*,*) '>>> dvinqn:: n',n
            CALL DVPDMP(X,V,N,'V(X)',6)
         end if
c
         jbar = 1
         j = 1
 100     continue
         if( jbar .gt. nbar ) go to 200
            if( xbar(jbar) .lt. x(j) ) then
               if( j .eq. 1 ) then
                  vbar(jbar) = vs
               else
                  if( leven(nintrp) ) then
                     i = min(n-nintrp+1,max(1,j-(nintrp/2)))
                  else
                     i = min(n-nintrp+1,max(1,j-(nintrp+1)/2 +
     *                       llogic(xbar(jbar) .gt.
     *                             0.5d0 * (x(j) + x(j-1)))))
                  end if
                  vbar(jbar) = flipn(xbar(jbar),x(i),v(i),nintrp)
               end if
               jbar = jbar + 1
            else if( xbar(jbar) .eq. x(j) ) then
               vbar(jbar) = v(j)
               jbar = jbar + 1
            else
               if( j .eq. n ) then
                  vbar(jbar) = vf
                  jbar = jbar + 1
               else
                  j = j + 1
               end if
            end if
         go to 100
c
 200     continue
c
         if( ltrace ) then
            WRITE(*,*) '>>> dvinqn:: nbar',nbar
            CALL DVPDMP(XBAR,VBAR,NBAR,'VBAR(XBAR)',6)
         end if

         return
c
      end

c-----------------------------------------------------------------------
c
c     Local version of EVEN to stop Loader messages.
c
c-----------------------------------------------------------------------
c
      logical function leven(n)
c
         integer      n
c
         leven = mod(n,2) .eq. 0
c
         return
c
      end

c-----------------------------------------------------------------------
c
c     "Converts" FORTRAN Boolean to {0,1} Boolean.
c
c     Local version to get rid of loader messages.
c
c-----------------------------------------------------------------------
c
      integer function llogic(tvalue)
         logical     tvalue
         if( tvalue ) then
            llogic = 1
         else
            llogic = 0
         end if
         return
      end

c-----------------------------------------------------------------------
c
c     Low level routine for (N-1)th (N > 1) order polynomial interp-
C     olation. Straightforward implementation of Nevilles algortihm.
c     Should only be used for reasonably well-conditioned problems.
c
c-----------------------------------------------------------------------
c
      double precision function flipn(xbar,x,y,n)
c
         implicit     logical (a-z)
c
         integer      nmax
         parameter    ( nmax = 20 )
         real*8       p(nmax)
c
         integer      n
         real*8       x(n), y(n)
         real*8       xbar
c
         integer      j, l
c
         do 10 l = 1 , n
            p(l) = y(l)
 10      continue
         do 30 l = 1 , n-1
            do 20 j = 1 , n-l
               p(j) = ((xbar - x(j+l)) * p(j) +
     *                 (x(j) - xbar)   * p(j+1)) / (x(j) - x(j+l))
 20         continue
 30      continue
         flipn = p(1)

        return

      end
