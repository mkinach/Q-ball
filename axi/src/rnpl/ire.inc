cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ire: manually update independent residuals using compactified
c          coordinates. Note that RNPL supplies the routine header,
c          declaration of some routine arguments, and RETURN/END
c          statements
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer i, j
      real*8  res

      ! excludes the boundaries
      do i=2, g1_NR-2
      do j=2, g1_NX-2

      ! logarithmic potential
      if ( Vtype == 1 ) then

      ! phi1 independent residual
      if ( (phi1_np1(i, j)**2 + phi2_np1(i, j)**2) .LT. 1.0d-50 ) then

      res = (phi1_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 - phi1_np1(i, j) * A_t_np1(i, j) **
     #2 + phi1_np1(i, j) * A_z_np1(i, j) ** 2 + phi1_np1(i, j) * (d * ex
     #p(c * R(i)) - d * exp(-c * R(i))) ** 2 * Atilde_phi_np1(i, j) ** 2
     #) * e ** 2 + (-0.2D1 * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #Atilde_p_np1(i, j) * (phi2_np1(i + 1, j) - phi2_np1(i, j)) / dR /
     #c / d / (exp(c * R(i)) + exp(-c * R(i))) + 0.2D1 * A_t_np1(i, j) *
     # pi2_np1(i, j) - 0.2D1 * A_z_np1(i, j) * (phi2_np1(i, j + 1) - phi
     #2_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j)))) * e
     #- ((phi1_np1(i + 2, j) - dble(2 * phi1_np1(i + 1, j)) + phi1_np1(i
     #, j)) / dR ** 2 / c / d / (exp(c * R(i)) + exp(-c * R(i))) - (dble
     #(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c * R(i
     #)) + exp(-c * R(i))) ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)
     #))) / c / d / (exp(c * R(i)) + exp(-c * R(i))) + (pi1_np1(i, j) -
     #pi1_n(i, j)) / dt - ((phi1_np1(i, j + 2) - dble(2 * phi1_np1(i, j
     #+ 1)) + phi1_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-
     #c * X(j))) - (dble(phi1_np1(i, j + 1)) - phi1_np1(i, j)) / dX / c
     #/ d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) -
     # c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) + exp(-c * X(j))) -
     # (dble(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) / (d * exp(c * R(i)) - d * exp(-c * R(i
     #)))

      else

      res = (phi1_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 - phi1_np1(i, j) * A_t_np1(i, j) **
     #2 + phi1_np1(i, j) * A_z_np1(i, j) ** 2 + phi1_np1(i, j) * (d * ex
     #p(c * R(i)) - d * exp(-c * R(i))) ** 2 * Atilde_phi_np1(i, j) ** 2
     #) * e ** 2 + (-0.2D1 * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #Atilde_p_np1(i, j) * (phi2_np1(i + 1, j) - phi2_np1(i, j)) / dR /
     #c / d / (exp(c * R(i)) + exp(-c * R(i))) + 0.2D1 * A_t_np1(i, j) *
     # pi2_np1(i, j) - 0.2D1 * A_z_np1(i, j) * (phi2_np1(i, j + 1) - phi
     #2_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j)))) * e
     #- ((phi1_np1(i + 2, j) - dble(2 * phi1_np1(i + 1, j)) + phi1_np1(i
     #, j)) / dR ** 2 / c / d / (exp(c * R(i)) + exp(-c * R(i))) - (dble
     #(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c * R(i
     #)) + exp(-c * R(i))) ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)
     #))) / c / d / (exp(c * R(i)) + exp(-c * R(i))) + (pi1_np1(i, j) -
     #pi1_n(i, j)) / dt - ((phi1_np1(i, j + 2) - dble(2 * phi1_np1(i, j
     #+ 1)) + phi1_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-
     #c * X(j))) - (dble(phi1_np1(i, j + 1)) - phi1_np1(i, j)) / dX / c
     #/ d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) -
     # c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) + exp(-c * X(j))) -
     # (dble(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) / (d * exp(c * R(i)) - d * exp(-c * R(i
     #))) - phi1_np1(i, j) * mu ** 2 * (log(beta ** 2 * (phi1_np1(i, j)
     #** 2 + phi2_np1(i, j) ** 2)) + 0.1D1)

      end if

      phi1ire_np1(i,j) = res

      ! phi2 independent residual
      if ( (phi1_np1(i, j)**2 + phi2_np1(i, j)**2) .LT. 1.0d-50 ) then

      res = (phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 + phi2_np1(i, j) * A_z_np1(i, j) ** 
     #2 + phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) ** 2
     # * Atilde_phi_np1(i, j) ** 2 - phi2_np1(i, j) * A_t_np1(i, j) ** 2
     #) * e ** 2 + (0.2D1 * A_z_np1(i, j) * (phi1_np1(i, j + 1) - phi1_n
     #p1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) - 0.2D1 
     #* A_t_np1(i, j) * pi1_np1(i, j) + 0.2D1 * (d * exp(c * R(i)) - d *
     # exp(-c * R(i))) * Atilde_p_np1(i, j) * (phi1_np1(i + 1, j) - phi1
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))) * e -
     # ((phi2_np1(i, j + 2) - dble(2 * phi2_np1(i, j + 1)) + phi2_np1(i,
     # j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-c * X(j))) - (dble(
     #phi2_np1(i, j + 1)) - phi2_np1(i, j)) / dX / c / d / (exp(c * X(j)
     #) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) - c * exp(-c * X(j))
     #)) / c / d / (exp(c * X(j)) + exp(-c * X(j))) - ((phi2_np1(i + 2, 
     #j) - dble(2 * phi2_np1(i + 1, j)) + phi2_np1(i, j)) / dR ** 2 / c 
     #/ d / (exp(c * R(i)) + exp(-c * R(i))) - (dble(phi2_np1(i + 1, j))
     # - phi2_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))
     # ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) + (pi
     #2_np1(i, j) - pi2_n(i, j)) / dt - (dble(phi2_np1(i + 1, j)) - phi2
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) / (d *
     # exp(c * R(i)) - d * exp(-c * R(i)))

      else

      res = (phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 + phi2_np1(i, j) * A_z_np1(i, j) ** 
     #2 + phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) ** 2
     # * Atilde_phi_np1(i, j) ** 2 - phi2_np1(i, j) * A_t_np1(i, j) ** 2
     #) * e ** 2 + (0.2D1 * A_z_np1(i, j) * (phi1_np1(i, j + 1) - phi1_n
     #p1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) - 0.2D1 
     #* A_t_np1(i, j) * pi1_np1(i, j) + 0.2D1 * (d * exp(c * R(i)) - d *
     # exp(-c * R(i))) * Atilde_p_np1(i, j) * (phi1_np1(i + 1, j) - phi1
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))) * e -
     # ((phi2_np1(i, j + 2) - dble(2 * phi2_np1(i, j + 1)) + phi2_np1(i,
     # j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-c * X(j))) - (dble(
     #phi2_np1(i, j + 1)) - phi2_np1(i, j)) / dX / c / d / (exp(c * X(j)
     #) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) - c * exp(-c * X(j))
     #)) / c / d / (exp(c * X(j)) + exp(-c * X(j))) - ((phi2_np1(i + 2, 
     #j) - dble(2 * phi2_np1(i + 1, j)) + phi2_np1(i, j)) / dR ** 2 / c 
     #/ d / (exp(c * R(i)) + exp(-c * R(i))) - (dble(phi2_np1(i + 1, j))
     # - phi2_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))
     # ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) - phi2_np1(i, j) * mu ** 2 * (log(beta 
     #** 2 * (phi1_np1(i, j) ** 2 + phi2_np1(i, j) ** 2)) + 0.1D1) + (pi
     #2_np1(i, j) - pi2_n(i, j)) / dt - (dble(phi2_np1(i + 1, j)) - phi2
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) / (d *
     # exp(c * R(i)) - d * exp(-c * R(i)))

      end if

      phi2ire_np1(i,j) = res

      ! A_t independent residual
      res = dble((2 * A_t_np1(i, j) * phi1_np1(i, j) ** 2 + 2 * A_t_np1(
     #i, j) * phi2_np1(i, j) ** 2) * e ** 2) + dble((-2 * phi1_np1(i, j)
     # * pi2_np1(i, j) + 2 * phi2_np1(i, j) * pi1_np1(i, j)) * e) + (B_t
     #_np1(i, j) - B_t_n(i, j)) / dt - (dble(A_t_np1(i, j + 2) - 2 * A_t
     #_np1(i, j + 1) + A_t_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(j))
     # + exp(-c * X(j))) - dble(A_t_np1(i, j + 1) - A_t_np1(i, j)) / dX
     #/ c / d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j
     #)) - c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) + exp(-c * X(j)
     #)) - (dble(A_t_np1(i + 2, j) - 2 * A_t_np1(i + 1, j) + A_t_np1(i,
     #j)) / dR ** 2 / c / d / (exp(c * R(i)) + exp(-c * R(i))) - dble(A_
     #t_np1(i + 1, j) - A_t_np1(i, j)) / dR / c / d / (exp(c * R(i)) + e
     #xp(-c * R(i))) ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) /
     #c / d / (exp(c * R(i)) + exp(-c * R(i))) - dble(A_t_np1(i + 1, j)
     #- A_t_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) /
     # (d * exp(c * R(i)) - d * exp(-c * R(i)))

      atire_np1(i,j) = res

      ! A_p independent residual
      res = (0.2D1 * (d * exp(c * R(i)) - d * exp(-c * R(i))) * Atilde_p
     #_np1(i, j) * phi1_np1(i, j) ** 2 + 0.2D1 * (d * exp(c * R(i)) - d 
     #* exp(-c * R(i))) * Atilde_p_np1(i, j) * phi2_np1(i, j) ** 2) * e 
     #** 2 + (-0.2D1 * phi1_np1(i, j) * (phi2_np1(i + 1, j) - phi2_np1(i
     #, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) + 0.2D1 * ph
     #i2_np1(i, j) * (phi1_np1(i + 1, j) - phi1_np1(i, j)) / dR / c / d 
     #/ (exp(c * R(i)) + exp(-c * R(i)))) * e - (d * exp(c * R(i)) - d *
     # exp(-c * R(i))) * ((Atilde_p_np1(i, j + 2) - dble(2 * Atilde_p_np
     #1(i, j + 1)) + Atilde_p_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(
     #j)) + exp(-c * X(j))) - (dble(Atilde_p_np1(i, j + 1)) - Atilde_p_n
     #p1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (
     #c * exp(c * X(j)) - c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) 
     #+ exp(-c * X(j))) - 0.2D1 * (Atilde_p_np1(i + 1, j) - Atilde_p_np1
     #(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) - (d * exp
     #(c * R(i)) - d * exp(-c * R(i))) * ((Atilde_p_np1(i + 2, j) - 0.2D
     #1 * Atilde_p_np1(i + 1, j) + Atilde_p_np1(i, j)) / dR ** 2 / c / d
     # / (exp(c * R(i)) + exp(-c * R(i))) - (Atilde_p_np1(i + 1, j) - At
     #ilde_p_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) 
     #** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c 
     #* R(i)) + exp(-c * R(i))) + (d * exp(c * R(i)) - d * exp(-c * R(i)
     #)) * (Btilde_p_np1(i, j) - Btilde_p_n(i, j)) / dt - (Atilde_p_np1(
     #i, j) + (d * exp(c * R(i)) - d * exp(-c * R(i))) * (Atilde_p_np1(i
     # + 1, j) - Atilde_p_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp
     #(-c * R(i)))) / (d * exp(c * R(i)) - d * exp(-c * R(i))) + 0.1D1 /
     # (d * exp(c * R(i)) - d * exp(-c * R(i))) * Atilde_p_np1(i, j)

      apire_np1(i,j) = res

      ! A_z independent residual
      res = dble((2 * A_z_np1(i, j) * phi1_np1(i, j) ** 2 + 2 * A_z_np1(
     #i, j) * phi2_np1(i, j) ** 2) * e ** 2) + (0.2D1 * dble(phi1_np1(i,
     # j + 1) - phi1_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c *
     # X(j))) * dble(phi2_np1(i, j)) - 0.2D1 * dble(phi2_np1(i, j + 1) -
     # phi2_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) *
     # dble(phi1_np1(i, j))) * dble(e) + (B_z_np1(i, j) - B_z_n(i, j)) /
     # dt - (dble(A_z_np1(i, j + 2) - 2 * A_z_np1(i, j + 1) + A_z_np1(i,
     # j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-c * X(j))) - dble(A
     #_z_np1(i, j + 1) - A_z_np1(i, j)) / dX / c / d / (exp(c * X(j)) + 
     #exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) - c * exp(-c * X(j)))) /
     # c / d / (exp(c * X(j)) + exp(-c * X(j))) - (dble(A_z_np1(i + 2, j
     #) - 2 * A_z_np1(i + 1, j) + A_z_np1(i, j)) / dR ** 2 / c / d / (ex
     #p(c * R(i)) + exp(-c * R(i))) - dble(A_z_np1(i + 1, j) - A_z_np1(i
     #, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) ** 2 * (c * 
     #exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c * R(i)) + ex
     #p(-c * R(i))) - dble(A_z_np1(i + 1, j) - A_z_np1(i, j)) / dR / c /
     # d / (exp(c * R(i)) + exp(-c * R(i))) / (d * exp(c * R(i)) - d * e
     #xp(-c * R(i)))

      azire_np1(i,j) = res

      ! polynomial potential
      else if ( Vtype == 2 ) then

      ! phi1 independent residual
      res = (phi1_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 - phi1_np1(i, j) * A_t_np1(i, j) ** 
     #2 + phi1_np1(i, j) * A_z_np1(i, j) ** 2 + phi1_np1(i, j) * (d * ex
     #p(c * R(i)) - d * exp(-c * R(i))) ** 2 * Atilde_phi_np1(i, j) ** 2
     #) * e ** 2 + (-0.2D1 * (d * exp(c * R(i)) - d * exp(-c * R(i))) * 
     #Atilde_p_np1(i, j) * (phi2_np1(i + 1, j) - phi2_np1(i, j)) / dR / 
     #c / d / (exp(c * R(i)) + exp(-c * R(i))) + 0.2D1 * A_t_np1(i, j) *
     # pi2_np1(i, j) - 0.2D1 * A_z_np1(i, j) * (phi2_np1(i, j + 1) - phi
     #2_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j)))) * e 
     #- ((phi1_np1(i + 2, j) - dble(2 * phi1_np1(i + 1, j)) + phi1_np1(i
     #, j)) / dR ** 2 / c / d / (exp(c * R(i)) + exp(-c * R(i))) - (dble
     #(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c * R(i
     #)) + exp(-c * R(i))) ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)
     #))) / c / d / (exp(c * R(i)) + exp(-c * R(i))) + (pi1_np1(i, j) - 
     #pi1_n(i, j)) / dt - ((phi1_np1(i, j + 2) - dble(2 * phi1_np1(i, j 
     #+ 1)) + phi1_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-
     #c * X(j))) - (dble(phi1_np1(i, j + 1)) - phi1_np1(i, j)) / dX / c 
     #/ d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) -
     # c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) + exp(-c * X(j))) -
     # (dble(phi1_np1(i + 1, j)) - phi1_np1(i, j)) / dR / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) / (d * exp(c * R(i)) - d * exp(-c * R(i
     #))) + phi1_np1(i, j) * (h * phi1_np1(i, j) ** 4 + (0.2D1 * h * phi
     #2_np1(i, j) ** 2 - gg) * phi1_np1(i, j) ** 2 + h * phi2_np1(i, j) 
     #** 4 - gg * phi2_np1(i, j) ** 2 + mm)

      phi1ire_np1(i,j) = res

      ! phi2 independent residual
      res = (phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) *
     #* 2 * Atilde_p_np1(i, j) ** 2 + phi2_np1(i, j) * A_z_np1(i, j) ** 
     #2 + phi2_np1(i, j) * (d * exp(c * R(i)) - d * exp(-c * R(i))) ** 2
     # * Atilde_phi_np1(i, j) ** 2 - phi2_np1(i, j) * A_t_np1(i, j) ** 2
     #) * e ** 2 + (0.2D1 * A_z_np1(i, j) * (phi1_np1(i, j + 1) - phi1_n
     #p1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) - 0.2D1 
     #* A_t_np1(i, j) * pi1_np1(i, j) + 0.2D1 * (d * exp(c * R(i)) - d *
     # exp(-c * R(i))) * Atilde_p_np1(i, j) * (phi1_np1(i + 1, j) - phi1
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))) * e -
     # ((phi2_np1(i, j + 2) - dble(2 * phi2_np1(i, j + 1)) + phi2_np1(i,
     # j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-c * X(j))) - (dble(
     #phi2_np1(i, j + 1)) - phi2_np1(i, j)) / dX / c / d / (exp(c * X(j)
     #) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) - c * exp(-c * X(j))
     #)) / c / d / (exp(c * X(j)) + exp(-c * X(j))) - ((phi2_np1(i + 2, 
     #j) - dble(2 * phi2_np1(i + 1, j)) + phi2_np1(i, j)) / dR ** 2 / c 
     #/ d / (exp(c * R(i)) + exp(-c * R(i))) - (dble(phi2_np1(i + 1, j))
     # - phi2_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))
     # ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c
     # * R(i)) + exp(-c * R(i))) + phi2_np1(i, j) * (h * phi2_np1(i, j) 
     #** 4 + (0.2D1 * h * phi1_np1(i, j) ** 2 - gg) * phi2_np1(i, j) ** 
     #2 + h * phi1_np1(i, j) ** 4 - gg * phi1_np1(i, j) ** 2 + mm) + (pi
     #2_np1(i, j) - pi2_n(i, j)) / dt - (dble(phi2_np1(i + 1, j)) - phi2
     #_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) / (d *
     # exp(c * R(i)) - d * exp(-c * R(i)))

      phi2ire_np1(i,j) = res

      ! A_t independent residual
      res = dble((2 * A_t_np1(i, j) * phi1_np1(i, j) ** 2 + 2 * A_t_np1(
     #i, j) * phi2_np1(i, j) ** 2) * e ** 2) + dble((-2 * phi1_np1(i, j)
     # * pi2_np1(i, j) + 2 * phi2_np1(i, j) * pi1_np1(i, j)) * e) + (B_t
     #_np1(i, j) - B_t_n(i, j)) / dt - (dble(A_t_np1(i, j + 2) - 2 * A_t
     #_np1(i, j + 1) + A_t_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(j))
     # + exp(-c * X(j))) - dble(A_t_np1(i, j + 1) - A_t_np1(i, j)) / dX 
     #/ c / d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (c * exp(c * X(j
     #)) - c * exp(-c * X(j)))) / c / d / (exp(c * X(j)) + exp(-c * X(j)
     #)) - (dble(A_t_np1(i + 2, j) - 2 * A_t_np1(i + 1, j) + A_t_np1(i, 
     #j)) / dR ** 2 / c / d / (exp(c * R(i)) + exp(-c * R(i))) - dble(A_
     #t_np1(i + 1, j) - A_t_np1(i, j)) / dR / c / d / (exp(c * R(i)) + e
     #xp(-c * R(i))) ** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / 
     #c / d / (exp(c * R(i)) + exp(-c * R(i))) - dble(A_t_np1(i + 1, j) 
     #- A_t_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) /
     # (d * exp(c * R(i)) - d * exp(-c * R(i)))

      atire_np1(i,j) = res

      ! A_p independent residual
      res = (0.2D1 * (d * exp(c * R(i)) - d * exp(-c * R(i))) * Atilde_p
     #_np1(i, j) * phi1_np1(i, j) ** 2 + 0.2D1 * (d * exp(c * R(i)) - d
     #* exp(-c * R(i))) * Atilde_p_np1(i, j) * phi2_np1(i, j) ** 2) * e
     #** 2 + (-0.2D1 * phi1_np1(i, j) * (phi2_np1(i + 1, j) - phi2_np1(i
     #, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) + 0.2D1 * ph
     #i2_np1(i, j) * (phi1_np1(i + 1, j) - phi1_np1(i, j)) / dR / c / d
     #/ (exp(c * R(i)) + exp(-c * R(i)))) * e - (d * exp(c * R(i)) - d *
     # exp(-c * R(i))) * ((Atilde_p_np1(i, j + 2) - dble(2 * Atilde_p_np
     #1(i, j + 1)) + Atilde_p_np1(i, j)) / dX ** 2 / c / d / (exp(c * X(
     #j)) + exp(-c * X(j))) - (dble(Atilde_p_np1(i, j + 1)) - Atilde_p_n
     #p1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) ** 2 * (
     #c * exp(c * X(j)) - c * exp(-c * X(j)))) / c / d / (exp(c * X(j))
     #+ exp(-c * X(j))) - 0.2D1 * (Atilde_p_np1(i + 1, j) - Atilde_p_np1
     #(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) - (d * exp
     #(c * R(i)) - d * exp(-c * R(i))) * ((Atilde_p_np1(i + 2, j) - 0.2D
     #1 * Atilde_p_np1(i + 1, j) + Atilde_p_np1(i, j)) / dR ** 2 / c / d
     # / (exp(c * R(i)) + exp(-c * R(i))) - (Atilde_p_np1(i + 1, j) - At
     #ilde_p_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i)))
     #** 2 * (c * exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c
     #* R(i)) + exp(-c * R(i))) + (d * exp(c * R(i)) - d * exp(-c * R(i)
     #)) * (Btilde_p_np1(i, j) - Btilde_p_n(i, j)) / dt - (Atilde_p_np1(
     #i, j) + (d * exp(c * R(i)) - d * exp(-c * R(i))) * (Atilde_p_np1(i
     # + 1, j) - Atilde_p_np1(i, j)) / dR / c / d / (exp(c * R(i)) + exp
     #(-c * R(i)))) / (d * exp(c * R(i)) - d * exp(-c * R(i))) + 0.1D1 /
     # (d * exp(c * R(i)) - d * exp(-c * R(i))) * Atilde_p_np1(i, j)

      apire_np1(i,j) = res

      ! A_z independent residual
      res = dble((2 * A_z_np1(i, j) * phi1_np1(i, j) ** 2 + 2 * A_z_np1(
     #i, j) * phi2_np1(i, j) ** 2) * e ** 2) + (0.2D1 * dble(phi1_np1(i,
     # j + 1) - phi1_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c *
     # X(j))) * dble(phi2_np1(i, j)) - 0.2D1 * dble(phi2_np1(i, j + 1) -
     # phi2_np1(i, j)) / dX / c / d / (exp(c * X(j)) + exp(-c * X(j))) *
     # dble(phi1_np1(i, j))) * dble(e) + (B_z_np1(i, j) - B_z_n(i, j)) /
     # dt - (dble(A_z_np1(i, j + 2) - 2 * A_z_np1(i, j + 1) + A_z_np1(i,
     # j)) / dX ** 2 / c / d / (exp(c * X(j)) + exp(-c * X(j))) - dble(A
     #_z_np1(i, j + 1) - A_z_np1(i, j)) / dX / c / d / (exp(c * X(j)) + 
     #exp(-c * X(j))) ** 2 * (c * exp(c * X(j)) - c * exp(-c * X(j)))) /
     # c / d / (exp(c * X(j)) + exp(-c * X(j))) - (dble(A_z_np1(i + 2, j
     #) - 2 * A_z_np1(i + 1, j) + A_z_np1(i, j)) / dR ** 2 / c / d / (ex
     #p(c * R(i)) + exp(-c * R(i))) - dble(A_z_np1(i + 1, j) - A_z_np1(i
     #, j)) / dR / c / d / (exp(c * R(i)) + exp(-c * R(i))) ** 2 * (c * 
     #exp(c * R(i)) - c * exp(-c * R(i)))) / c / d / (exp(c * R(i)) + ex
     #p(-c * R(i))) - dble(A_z_np1(i + 1, j) - A_z_np1(i, j)) / dR / c /
     # d / (exp(c * R(i)) + exp(-c * R(i))) / (d * exp(c * R(i)) - d * e
     #xp(-c * R(i)))

      azire_np1(i,j) = res

      end if

      end do
      end do
