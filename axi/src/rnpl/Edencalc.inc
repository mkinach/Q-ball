cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Edencalc: manually update the energy density. Note that RNPL
c               supplies the routine header, declaration of some routine
c               arguments, and RETURN/END statements
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      integer i, j
      real*8  res

      ! excludes the boundaries
      do i=2, g1_NR-1
      do j=2, g1_NX-1

      ! logarithmic potential
      if ( Vtype == 1 ) then

      if ( (phi1_np1(i,j)**2.0 + phi2_np1(i,j)**2.0) .lt. (1D-50) ) then

      res = pR_n(i, j) ** 2 * Btilde_p_np1(i, j) ** 2 / 0.2D1 + pR_n(i, 
     #j) ** 2 * Btilde_phi_np1(i, j) ** 2 / 0.2D1 + pi1_np1(i, j) ** 2 +
     # pi2_np1(i, j) ** 2 + dble(2 * Atilde_phi_np1(i, j) ** 2) + B_z_np
     #1(i, j) ** 2 / 0.2D1 + xi_np1(i, j) ** 2 + pR_n(i, j) * dble(Atild
     #e_phi_np1(i, j)) * (Atilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 
     #1, j)) / dR * dRdp_n(i, j) - pR_n(i, j) * Btilde_p_np1(i, j) * (A_
     #t_np1(i + 1, j) - A_t_np1(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 +
     # A_t_np1(i, j) ** 2 * e ** 2 * phi2_np1(i, j) ** 2 + A_t_np1(i, j)
     # ** 2 * e ** 2 * phi1_np1(i, j) ** 2 + A_z_np1(i, j) ** 2 * e ** 2
     # * phi1_np1(i, j) ** 2 + A_z_np1(i, j) ** 2 * e ** 2 * phi2_np1(i,
     # j) ** 2 + cc * chi_np1(i, j) ** 2 * phi1_np1(i, j) ** 2 + cc * ch
     #i_np1(i, j) ** 2 * phi2_np1(i, j) ** 2 + (phi1_np1(i + 1, j) - phi
     #1_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (A_z
     #_np1(i + 1, j) - A_z_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) 
     #** 2 / 0.8D1 + (phi2_np1(i + 1, j) - phi2_np1(i - 1, j)) ** 2 / dR
     # ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (chi_np1(i + 1, j) - chi_np1(i
     # - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (A_t_np1(i 
     #+ 1, j) - A_t_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 
     #0.8D1 + (chi_np1(i, j + 1) - chi_np1(i, j - 1)) ** 2 / dX ** 2 * d
     #Xdz_n(i, j) ** 2 / 0.4D1 + (phi1_np1(i, j + 1) - phi1_np1(i, j - 1
     #)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1 + (phi2_np1(i, j + 1
     #) - phi2_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1
     # + (A_t_np1(i, j + 1) - A_t_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n
     #(i, j) ** 2 / 0.8D1 + A_z_np1(i, j) * e * (phi1_np1(i, j + 1) - ph
     #i1_np1(i, j - 1)) / dX * dXdz_n(i, j) * phi2_np1(i, j) - A_z_np1(i
     #, j) * e * (phi2_np1(i, j + 1) - phi2_np1(i, j - 1)) / dX * dXdz_n
     #(i, j) * phi1_np1(i, j) - pR_n(i, j) * Atilde_p_np1(i, j) * e * (p
     #hi2_np1(i + 1, j) - phi2_np1(i - 1, j)) / dR * dRdp_n(i, j) * phi1
     #_np1(i, j) + pR_n(i, j) * Atilde_p_np1(i, j) * e * (phi1_np1(i + 1
     #, j) - phi1_np1(i - 1, j)) / dR * dRdp_n(i, j) * phi2_np1(i, j) - 
     #pR_n(i, j) * (A_z_np1(i + 1, j) - A_z_np1(i - 1, j)) / dR * dRdp_n
     #(i, j) * (Atilde_p_np1(i, j + 1) - Atilde_p_np1(i, j - 1)) / dX * 
     #dXdz_n(i, j) / 0.4D1 + 0.2D1 * A_t_np1(i, j) * e * pi1_np1(i, j) *
     # phi2_np1(i, j) - 0.2D1 * A_t_np1(i, j) * e * pi2_np1(i, j) * phi1
     #_np1(i, j) + pR_n(i, j) ** 2 * Atilde_p_np1(i, j) ** 2 * e ** 2 * 
     #phi1_np1(i, j) ** 2 + pR_n(i, j) ** 2 * Atilde_p_np1(i, j) ** 2 * 
     #e ** 2 * phi2_np1(i, j) ** 2 + pR_n(i, j) ** 2 * e ** 2 * dble(Ati
     #lde_phi_np1(i, j) ** 2) * phi1_np1(i, j) ** 2 + pR_n(i, j) ** 2 * 
     #e ** 2 * dble(Atilde_phi_np1(i, j) ** 2) * phi2_np1(i, j) ** 2 + p
     #R_n(i, j) ** 2 * (Atilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 1,
     # j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.8D1 - B_z_np1(i, j) * 
     #(A_t_np1(i, j + 1) - A_t_np1(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D
     #1 + pR_n(i, j) ** 2 * (Atilde_p_np1(i, j + 1) - Atilde_p_np1(i, j 
     #- 1)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.8D1 + pR_n(i, j) ** 2
     # * (Atilde_phi_np1(i, j + 1) - Atilde_phi_np1(i, j - 1)) ** 2 / dX
     # ** 2 * dXdz_n(i, j) ** 2 / 0.8D1

      else

      res = pR_n(i, j) ** 2 * Btilde_p_np1(i, j) ** 2 / 0.2D1 + pR_n(i, 
     #j) ** 2 * Btilde_phi_np1(i, j) ** 2 / 0.2D1 + pi1_np1(i, j) ** 2 +
     # pi2_np1(i, j) ** 2 + dble(2 * Atilde_phi_np1(i, j) ** 2) + B_z_np
     #1(i, j) ** 2 / 0.2D1 + xi_np1(i, j) ** 2 + pR_n(i, j) * dble(Atild
     #e_phi_np1(i, j)) * (Atilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 
     #1, j)) / dR * dRdp_n(i, j) - pR_n(i, j) * Btilde_p_np1(i, j) * (A_
     #t_np1(i + 1, j) - A_t_np1(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 -
     # mu ** 2 * log(beta ** 2 * (phi1_np1(i, j) ** 2 + phi2_np1(i, j) *
     #* 2)) * phi2_np1(i, j) ** 2 + A_t_np1(i, j) ** 2 * e ** 2 * phi2_n
     #p1(i, j) ** 2 - mu ** 2 * log(beta ** 2 * (phi1_np1(i, j) ** 2 + p
     #hi2_np1(i, j) ** 2)) * phi1_np1(i, j) ** 2 + A_t_np1(i, j) ** 2 * 
     #e ** 2 * phi1_np1(i, j) ** 2 + A_z_np1(i, j) ** 2 * e ** 2 * phi1_
     #np1(i, j) ** 2 + A_z_np1(i, j) ** 2 * e ** 2 * phi2_np1(i, j) ** 2
     # + cc * chi_np1(i, j) ** 2 * phi1_np1(i, j) ** 2 + cc * chi_np1(i,
     # j) ** 2 * phi2_np1(i, j) ** 2 + (phi1_np1(i + 1, j) - phi1_np1(i 
     #- 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (A_z_np1(i +
     # 1, j) - A_z_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0
     #.8D1 + (phi2_np1(i + 1, j) - phi2_np1(i - 1, j)) ** 2 / dR ** 2 * 
     #dRdp_n(i, j) ** 2 / 0.4D1 + (chi_np1(i + 1, j) - chi_np1(i - 1, j)
     #) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (A_t_np1(i + 1, j) 
     #- A_t_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.8D1 + 
     #(chi_np1(i, j + 1) - chi_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i,
     # j) ** 2 / 0.4D1 + (phi1_np1(i, j + 1) - phi1_np1(i, j - 1)) ** 2 
     #/ dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1 + (phi2_np1(i, j + 1) - phi2
     #_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1 + (A_t_
     #np1(i, j + 1) - A_t_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j) *
     #* 2 / 0.8D1 + A_z_np1(i, j) * e * (phi1_np1(i, j + 1) - phi1_np1(i
     #, j - 1)) / dX * dXdz_n(i, j) * phi2_np1(i, j) - A_z_np1(i, j) * e
     # * (phi2_np1(i, j + 1) - phi2_np1(i, j - 1)) / dX * dXdz_n(i, j) *
     # phi1_np1(i, j) - pR_n(i, j) * Atilde_p_np1(i, j) * e * (phi2_np1(
     #i + 1, j) - phi2_np1(i - 1, j)) / dR * dRdp_n(i, j) * phi1_np1(i, 
     #j) + pR_n(i, j) * Atilde_p_np1(i, j) * e * (phi1_np1(i + 1, j) - p
     #hi1_np1(i - 1, j)) / dR * dRdp_n(i, j) * phi2_np1(i, j) - pR_n(i, 
     #j) * (A_z_np1(i + 1, j) - A_z_np1(i - 1, j)) / dR * dRdp_n(i, j) *
     # (Atilde_p_np1(i, j + 1) - Atilde_p_np1(i, j - 1)) / dX * dXdz_n(i
     #, j) / 0.4D1 + 0.2D1 * A_t_np1(i, j) * e * pi1_np1(i, j) * phi2_np
     #1(i, j) - 0.2D1 * A_t_np1(i, j) * e * pi2_np1(i, j) * phi1_np1(i, 
     #j) + pR_n(i, j) ** 2 * Atilde_p_np1(i, j) ** 2 * e ** 2 * phi1_np1
     #(i, j) ** 2 + pR_n(i, j) ** 2 * Atilde_p_np1(i, j) ** 2 * e ** 2 *
     # phi2_np1(i, j) ** 2 + pR_n(i, j) ** 2 * e ** 2 * dble(Atilde_phi_
     #np1(i, j) ** 2) * phi1_np1(i, j) ** 2 + pR_n(i, j) ** 2 * e ** 2 *
     # dble(Atilde_phi_np1(i, j) ** 2) * phi2_np1(i, j) ** 2 + pR_n(i, j
     #) ** 2 * (Atilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 1, j)) ** 
     #2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.8D1 - B_z_np1(i, j) * (A_t_np1
     #(i, j + 1) - A_t_np1(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1 + pR_n
     #(i, j) ** 2 * (Atilde_p_np1(i, j + 1) - Atilde_p_np1(i, j - 1)) **
     # 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.8D1 + pR_n(i, j) ** 2 * (Atil
     #de_phi_np1(i, j + 1) - Atilde_phi_np1(i, j - 1)) ** 2 / dX ** 2 * 
     #dXdz_n(i, j) ** 2 / 0.8D1

      end if

      ! polynomial potential
      else if ( Vtype == 2 ) then

      res = pi1_np1(i, j) ** 2 + pi2_np1(i, j) ** 2 + dble(2 * Atilde_ph
     #i_np1(i, j) ** 2) + B_z_np1(i, j) ** 2 / 0.2D1 + xi_np1(i, j) ** 2
     # + mm * (phi1_np1(i, j) ** 2 + phi2_np1(i, j) ** 2) - gg * (phi1_n
     #p1(i, j) ** 2 + phi2_np1(i, j) ** 2) ** 2 / 0.2D1 + h * (phi1_np1(
     #i, j) ** 2 + phi2_np1(i, j) ** 2) ** 3 / 0.3D1 + pR_n(i, j) ** 2 *
     # Btilde_p_np1(i, j) ** 2 / 0.2D1 + pR_n(i, j) ** 2 * Btilde_phi_np
     #1(i, j) ** 2 / 0.2D1 + pR_n(i, j) * dble(Atilde_phi_np1(i, j)) * (
     #Atilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 1, j)) / dR * dRdp_n
     #(i, j) - pR_n(i, j) * Btilde_p_np1(i, j) * (A_t_np1(i + 1, j) - A_
     #t_np1(i - 1, j)) / dR * dRdp_n(i, j) / 0.2D1 + 0.2D1 * A_t_np1(i, 
     #j) * e * pi1_np1(i, j) * phi2_np1(i, j) - 0.2D1 * A_t_np1(i, j) * 
     #e * pi2_np1(i, j) * phi1_np1(i, j) + pR_n(i, j) * Atilde_p_np1(i, 
     #j) * e * (phi1_np1(i + 1, j) - phi1_np1(i - 1, j)) / dR * dRdp_n(i
     #, j) * phi2_np1(i, j) - pR_n(i, j) * (A_z_np1(i + 1, j) - A_z_np1(
     #i - 1, j)) / dR * dRdp_n(i, j) * (Atilde_p_np1(i, j + 1) - Atilde_
     #p_np1(i, j - 1)) / dX * dXdz_n(i, j) / 0.4D1 - pR_n(i, j) * Atilde
     #_p_np1(i, j) * e * (phi2_np1(i + 1, j) - phi2_np1(i - 1, j)) / dR 
     #* dRdp_n(i, j) * phi1_np1(i, j) + A_z_np1(i, j) * e * phi2_np1(i, 
     #j) * (phi1_np1(i, j + 1) - phi1_np1(i, j - 1)) / dX * dXdz_n(i, j)
     # - A_z_np1(i, j) * e * phi1_np1(i, j) * (phi2_np1(i, j + 1) - phi2
     #_np1(i, j - 1)) / dX * dXdz_n(i, j) + A_t_np1(i, j) ** 2 * e ** 2 
     #* phi1_np1(i, j) ** 2 + (chi_np1(i + 1, j) - chi_np1(i - 1, j)) **
     # 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.4D1 + (phi1_np1(i, j + 1) - p
     #hi1_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1 + (c
     #hi_np1(i, j + 1) - chi_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j
     #) ** 2 / 0.4D1 + (phi2_np1(i, j + 1) - phi2_np1(i, j - 1)) ** 2 / 
     #dX ** 2 * dXdz_n(i, j) ** 2 / 0.4D1 + (A_t_np1(i, j + 1) - A_t_np1
     #(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j) ** 2 / 0.8D1 + A_z_np1(i
     #, j) ** 2 * e ** 2 * phi1_np1(i, j) ** 2 + cc * chi_np1(i, j) ** 2
     # * phi1_np1(i, j) ** 2 + A_t_np1(i, j) ** 2 * e ** 2 * phi2_np1(i,
     # j) ** 2 + (phi1_np1(i + 1, j) - phi1_np1(i - 1, j)) ** 2 / dR ** 
     #2 * dRdp_n(i, j) ** 2 / 0.4D1 + (A_z_np1(i + 1, j) - A_z_np1(i - 1
     #, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.8D1 + (A_t_np1(i + 1,
     # j) - A_t_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2 / 0.8D
     #1 + A_z_np1(i, j) ** 2 * e ** 2 * phi2_np1(i, j) ** 2 + (phi2_np1(
     #i + 1, j) - phi2_np1(i - 1, j)) ** 2 / dR ** 2 * dRdp_n(i, j) ** 2
     # / 0.4D1 + cc * chi_np1(i, j) ** 2 * phi2_np1(i, j) ** 2 + pR_n(i,
     # j) ** 2 * Atilde_p_np1(i, j) ** 2 * e ** 2 * phi1_np1(i, j) ** 2 
     #+ pR_n(i, j) ** 2 * Atilde_p_np1(i, j) ** 2 * e ** 2 * phi2_np1(i,
     # j) ** 2 + pR_n(i, j) ** 2 * e ** 2 * dble(Atilde_phi_np1(i, j) **
     # 2) * phi1_np1(i, j) ** 2 + pR_n(i, j) ** 2 * e ** 2 * dble(Atilde
     #_phi_np1(i, j) ** 2) * phi2_np1(i, j) ** 2 + pR_n(i, j) ** 2 * (At
     #ilde_phi_np1(i + 1, j) - Atilde_phi_np1(i - 1, j)) ** 2 / dR ** 2 
     #* dRdp_n(i, j) ** 2 / 0.8D1 - B_z_np1(i, j) * (A_t_np1(i, j + 1) -
     # A_t_np1(i, j - 1)) / dX * dXdz_n(i, j) / 0.2D1 + pR_n(i, j) ** 2 
     #* (Atilde_phi_np1(i, j + 1) - Atilde_phi_np1(i, j - 1)) ** 2 / dX 
     #** 2 * dXdz_n(i, j) ** 2 / 0.8D1 + pR_n(i, j) ** 2 * (Atilde_p_np1
     #(i, j + 1) - Atilde_p_np1(i, j - 1)) ** 2 / dX ** 2 * dXdz_n(i, j)
     # ** 2 / 0.8D1

      end if

      Eden_np1(i,j) = res

      end do
      end do
