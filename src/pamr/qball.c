//======================================================================
// qball.c
//
// This program interfaces with the PAMR and AMRD libraries to solve the
// gauged Q-ball equations of motion using adaptive mesh refinement and
// multigrid in three spatial dimensions.
//
// For more information about built-in PAMR/AMRD functionality, look at
// the PAMR or AMRD reference manual, pamr.h, amrd.h, or pamr.c
//======================================================================

// #include <fenv.h> // uncomment when debugging FPEs
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <amrd.h>
#include <pamr.h>
#include "qball.h"

#include <sys/timeb.h>
void elapsed_time(void);

// model parameters (values set in input parameter files)
int anti, V_type, ini_type;
real e, m, g, h, beta, mu;
real w, w2;
real cx, cy, cz, cx1, cx2, cy1, cy2, cz1, cz2;
real v, phase;
real c, d;
real cc, amp, chiX0, chiY0, chiZ0, wwX, wwY, wwZ, r0, delta;
real myzero;

// pointers for physical coordinate values
real *x, *y, *z;

// pointers for grid function values
real *nm1_phi1, *n_phi1, *np1_phi1;
real *nm1_phi2, *n_phi2, *np1_phi2;
real *nm1_pi1, *n_pi1, *np1_pi1;
real *nm1_pi2, *n_pi2, *np1_pi2;
real *nm1_At, *n_At, *np1_At;
real *nm1_Ax, *n_Ax, *np1_Ax;
real *nm1_Ay, *n_Ay, *np1_Ay;
real *nm1_Az, *n_Az, *np1_Az;
real *nm1_Bt, *n_Bt, *np1_Bt;
real *nm1_Bx, *n_Bx, *np1_Bx;
real *nm1_By, *n_By, *np1_By;
real *nm1_Bz, *n_Bz, *np1_Bz;
real *nm1_chi, *n_chi, *np1_chi;
real *nm1_xi, *n_xi, *np1_xi;
real *nm1_modphi, *n_modphi, *np1_modphi;
real *nm1_qden, *n_qden, *np1_qden;
real *nm1_eden, *n_eden, *np1_eden;
real *nm1_jdenx, *n_jdenx, *np1_jdenx;
real *nm1_lorenz, *n_lorenz, *np1_lorenz;
real *nm1_elecx, *n_elecx, *np1_elecx;
real *nm1_elecy, *n_elecy, *np1_elecy;
real *nm1_elecz, *n_elecz, *np1_elecz;
real *nm1_magx, *n_magx, *np1_magx;
real *nm1_magy, *n_magy, *np1_magy;
real *nm1_magz, *n_magz, *np1_magz;
real *nm1_eem, *n_eem, *np1_eem;
real *nm1_gaussE, *n_gaussE, *np1_gaussE;
real *nm1_gaussB, *n_gaussB, *np1_gaussB;

real *nm1_Xx, *n_Xx, *np1_Xx;
real *nm1_dxdX, *n_dxdX, *np1_dxdX;
real *nm1_dxdX2, *n_dxdX2, *np1_dxdX2;
real *nm1_Yy, *n_Yy, *np1_Yy;
real *nm1_dydY, *n_dydY, *np1_dydY;
real *nm1_dydY2, *n_dydY2, *np1_dydY2;
real *nm1_Zz, *n_Zz, *np1_Zz;
real *nm1_dzdZ, *n_dzdZ, *np1_dzdZ;
real *nm1_dzdZ2, *n_dzdZ2, *np1_dzdZ2;
real *nm1_jac, *n_jac, *np1_jac;

real *cmask;

real *phi1_k1, *phi1_k2, *phi1_k3, *phi1_k4;
real *phi2_k1, *phi2_k2, *phi2_k3, *phi2_k4;
real *pi1_k1, *pi1_k2, *pi1_k3, *pi1_k4;
real *pi2_k1, *pi2_k2, *pi2_k3, *pi2_k4;
real *At_k1, *At_k2, *At_k3, *At_k4;
real *Ax_k1, *Ax_k2, *Ax_k3, *Ax_k4;
real *Ay_k1, *Ay_k2, *Ay_k3, *Ay_k4;
real *Az_k1, *Az_k2, *Az_k3, *Az_k4;
real *Bt_k1, *Bt_k2, *Bt_k3, *Bt_k4;
real *Bx_k1, *Bx_k2, *Bx_k3, *Bx_k4;
real *By_k1, *By_k2, *By_k3, *By_k4;
real *Bz_k1, *Bz_k2, *Bz_k3, *Bz_k4;
real *chi_k1, *chi_k2, *chi_k3, *chi_k4;
real *xi_k1, *xi_k2, *xi_k3, *xi_k4;

// pointers for multigrid functionality
real *V, *V_lop, *V_res, *V_rhs;
real *nm1_V, *n_V, *np1_V;
real *defect;
real *mask, *mask_mg;
real *Xx;
real *dxdX;
real *dxdX2;
real *Yy;
real *dydY;
real *dydY2;
real *Zz;
real *dzdZ;
real *dzdZ2;
real *phi1, *phi2, *pi1, *pi2, *modphi, *At;

// grid function numbers (gfn) for evolved and derived quantities
int nm1_phi1_gfn, n_phi1_gfn, np1_phi1_gfn;
int nm1_phi2_gfn, n_phi2_gfn, np1_phi2_gfn;
int nm1_pi1_gfn, n_pi1_gfn, np1_pi1_gfn;
int nm1_pi2_gfn, n_pi2_gfn, np1_pi2_gfn;
int nm1_At_gfn, n_At_gfn, np1_At_gfn;
int nm1_Ax_gfn, n_Ax_gfn, np1_Ax_gfn;
int nm1_Ay_gfn, n_Ay_gfn, np1_Ay_gfn;
int nm1_Az_gfn, n_Az_gfn, np1_Az_gfn;
int nm1_Bt_gfn, n_Bt_gfn, np1_Bt_gfn;
int nm1_Bx_gfn, n_Bx_gfn, np1_Bx_gfn;
int nm1_By_gfn, n_By_gfn, np1_By_gfn;
int nm1_Bz_gfn, n_Bz_gfn, np1_Bz_gfn;
int nm1_chi_gfn, n_chi_gfn, np1_chi_gfn;
int nm1_xi_gfn, n_xi_gfn, np1_xi_gfn;
int nm1_modphi_gfn, n_modphi_gfn, np1_modphi_gfn;
int nm1_qden_gfn, n_qden_gfn, np1_qden_gfn;
int nm1_eden_gfn, n_eden_gfn, np1_eden_gfn;
int nm1_jdenx_gfn, n_jdenx_gfn, np1_jdenx_gfn;
int nm1_lorenz_gfn, n_lorenz_gfn, np1_lorenz_gfn;
int nm1_elecx_gfn, n_elecx_gfn, np1_elecx_gfn;
int nm1_elecy_gfn, n_elecy_gfn, np1_elecy_gfn;
int nm1_elecz_gfn, n_elecz_gfn, np1_elecz_gfn;
int nm1_magx_gfn, n_magx_gfn, np1_magx_gfn;
int nm1_magy_gfn, n_magy_gfn, np1_magy_gfn;
int nm1_magz_gfn, n_magz_gfn, np1_magz_gfn;
int nm1_eem_gfn, n_eem_gfn, np1_eem_gfn;
int nm1_gaussE_gfn, n_gaussE_gfn, np1_gaussE_gfn;
int nm1_gaussB_gfn, n_gaussB_gfn, np1_gaussB_gfn;

// grid function numbers (gfn) for coordinate quantities
int nm1_Xx_gfn, n_Xx_gfn, np1_Xx_gfn;
int nm1_dxdX_gfn, n_dxdX_gfn, np1_dxdX_gfn;
int nm1_dxdX2_gfn, n_dxdX2_gfn, np1_dxdX2_gfn;
int nm1_Yy_gfn, n_Yy_gfn, np1_Yy_gfn;
int nm1_dydY_gfn, n_dydY_gfn, np1_dydY_gfn;
int nm1_dydY2_gfn, n_dydY2_gfn, np1_dydY2_gfn;
int nm1_Zz_gfn, n_Zz_gfn, np1_Zz_gfn;
int nm1_dzdZ_gfn, n_dzdZ_gfn, np1_dzdZ_gfn;
int nm1_dzdZ2_gfn, n_dzdZ2_gfn, np1_dzdZ2_gfn;
int nm1_jac_gfn, n_jac_gfn, np1_jac_gfn;

// grid function number (gfn) for AMR child mask
int cmask_gfn;

// grid function numbers (gfn) for RK4 scheme
int phi1_k1_gfn, phi1_k2_gfn, phi1_k3_gfn, phi1_k4_gfn;
int phi2_k1_gfn, phi2_k2_gfn, phi2_k3_gfn, phi2_k4_gfn;
int pi1_k1_gfn, pi1_k2_gfn, pi1_k3_gfn, pi1_k4_gfn;
int pi2_k1_gfn, pi2_k2_gfn, pi2_k3_gfn, pi2_k4_gfn;
int At_k1_gfn, At_k2_gfn, At_k3_gfn, At_k4_gfn;
int Ax_k1_gfn, Ax_k2_gfn, Ax_k3_gfn, Ax_k4_gfn;
int Ay_k1_gfn, Ay_k2_gfn, Ay_k3_gfn, Ay_k4_gfn;
int Az_k1_gfn, Az_k2_gfn, Az_k3_gfn, Az_k4_gfn;
int Bt_k1_gfn, Bt_k2_gfn, Bt_k3_gfn, Bt_k4_gfn;
int Bx_k1_gfn, Bx_k2_gfn, Bx_k3_gfn, Bx_k4_gfn;
int By_k1_gfn, By_k2_gfn, By_k3_gfn, By_k4_gfn;
int Bz_k1_gfn, Bz_k2_gfn, Bz_k3_gfn, Bz_k4_gfn;
int chi_k1_gfn, chi_k2_gfn, chi_k3_gfn, chi_k4_gfn;
int xi_k1_gfn, xi_k2_gfn, xi_k3_gfn, xi_k4_gfn;

// grid function numbers (gfn) for multigrid
int V_gfn, V_lop_gfn, V_res_gfn, V_rhs_gfn, mask_gfn, mask_mg_gfn;
int nm1_V_gfn, n_V_gfn, np1_V_gfn;
int defect_gfn;
int Xx_gfn, dxdX_gfn, dxdX2_gfn;
int Yy_gfn, dydY_gfn, dydY2_gfn;
int Zz_gfn, dzdZ_gfn, dzdZ2_gfn;
int phi1_gfn, phi2_gfn, pi1_gfn, pi2_gfn, modphi_gfn, At_gfn;

// toggles for multigrid functionality
int mg_toggle;
int MGiter = 0;
int MGshape = 0;
char mg_data[] = "V.bin";

int dim;             // grid dimensionality
int size;            // grid size
int shape[3];        // number of grid points per dimension
int ghost_width[6];  // width of AMR ghost region per dimension
int phys_bdy[6];     // toggle for AMR/physical boundaries
int g_rank;          // MPI rank of execution
int g_L;             // current level in AMR hierarchy
int Nx, Ny, Nz;      // number of grid points in each dimension
real dx, dy, dz, dt; // grid spacing in each dimension
real base_bbox[6];   // stores physical bounding box data
real bbox[6];        // stores AMR bounding box data
real *g_norms;       // stores L-infinity norm

/***********************************************************************
 * set_gfns
 *
 * Maps grid functions (e.g., phi1, phi2, ...) to their respective grid
 * function numbers (gfn) with error checking. Specifically,
 * PAMR_get_gfn() returns the gfn of a grid function in the hierarchy
 * defined by the PAMR_AMRH variable (from pamr.h) at the appropriate
 * time level.
 *
 **********************************************************************/
void set_gfns(void) {
  if ((nm1_phi1_gfn = PAMR_get_gfn("phi1", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_phi1_gfn = PAMR_get_gfn("phi1", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_phi1_gfn = PAMR_get_gfn("phi1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_phi2_gfn = PAMR_get_gfn("phi2", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_phi2_gfn = PAMR_get_gfn("phi2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_phi2_gfn = PAMR_get_gfn("phi2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_pi1_gfn = PAMR_get_gfn("pi1", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_pi1_gfn = PAMR_get_gfn("pi1", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_pi1_gfn = PAMR_get_gfn("pi1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_pi2_gfn = PAMR_get_gfn("pi2", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_pi2_gfn = PAMR_get_gfn("pi2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_pi2_gfn = PAMR_get_gfn("pi2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_At_gfn = PAMR_get_gfn("At", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_At_gfn = PAMR_get_gfn("At", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_At_gfn = PAMR_get_gfn("At", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Ax_gfn = PAMR_get_gfn("Ax", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Ax_gfn = PAMR_get_gfn("Ax", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Ax_gfn = PAMR_get_gfn("Ax", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Ay_gfn = PAMR_get_gfn("Ay", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Ay_gfn = PAMR_get_gfn("Ay", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Ay_gfn = PAMR_get_gfn("Ay", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Az_gfn = PAMR_get_gfn("Az", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Az_gfn = PAMR_get_gfn("Az", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Az_gfn = PAMR_get_gfn("Az", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Bt_gfn = PAMR_get_gfn("Bt", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Bt_gfn = PAMR_get_gfn("Bt", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Bt_gfn = PAMR_get_gfn("Bt", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Bx_gfn = PAMR_get_gfn("Bx", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Bx_gfn = PAMR_get_gfn("Bx", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Bx_gfn = PAMR_get_gfn("Bx", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_By_gfn = PAMR_get_gfn("By", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_By_gfn = PAMR_get_gfn("By", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_By_gfn = PAMR_get_gfn("By", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Bz_gfn = PAMR_get_gfn("Bz", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Bz_gfn = PAMR_get_gfn("Bz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Bz_gfn = PAMR_get_gfn("Bz", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_chi_gfn = PAMR_get_gfn("chi", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_chi_gfn = PAMR_get_gfn("chi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_chi_gfn = PAMR_get_gfn("chi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_xi_gfn = PAMR_get_gfn("xi", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_xi_gfn = PAMR_get_gfn("xi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_xi_gfn = PAMR_get_gfn("xi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_modphi_gfn = PAMR_get_gfn("modphi", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_modphi_gfn = PAMR_get_gfn("modphi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_modphi_gfn = PAMR_get_gfn("modphi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_qden_gfn = PAMR_get_gfn("qden", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_qden_gfn = PAMR_get_gfn("qden", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_qden_gfn = PAMR_get_gfn("qden", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_eden_gfn = PAMR_get_gfn("eden", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_eden_gfn = PAMR_get_gfn("eden", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_eden_gfn = PAMR_get_gfn("eden", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_jdenx_gfn = PAMR_get_gfn("jdenx", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_jdenx_gfn = PAMR_get_gfn("jdenx", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_jdenx_gfn = PAMR_get_gfn("jdenx", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_lorenz_gfn = PAMR_get_gfn("lorenz", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_lorenz_gfn = PAMR_get_gfn("lorenz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_lorenz_gfn = PAMR_get_gfn("lorenz", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_elecx_gfn = PAMR_get_gfn("elecx", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_elecx_gfn = PAMR_get_gfn("elecx", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_elecx_gfn = PAMR_get_gfn("elecx", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_elecy_gfn = PAMR_get_gfn("elecy", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_elecy_gfn = PAMR_get_gfn("elecy", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_elecy_gfn = PAMR_get_gfn("elecy", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_elecz_gfn = PAMR_get_gfn("elecz", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_elecz_gfn = PAMR_get_gfn("elecz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_elecz_gfn = PAMR_get_gfn("elecz", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_magx_gfn = PAMR_get_gfn("magx", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_magx_gfn = PAMR_get_gfn("magx", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_magx_gfn = PAMR_get_gfn("magx", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_magy_gfn = PAMR_get_gfn("magy", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_magy_gfn = PAMR_get_gfn("magy", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_magy_gfn = PAMR_get_gfn("magy", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_magz_gfn = PAMR_get_gfn("magz", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_magz_gfn = PAMR_get_gfn("magz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_magz_gfn = PAMR_get_gfn("magz", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_eem_gfn = PAMR_get_gfn("eem", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_eem_gfn = PAMR_get_gfn("eem", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_eem_gfn = PAMR_get_gfn("eem", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_gaussE_gfn = PAMR_get_gfn("gaussE", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_gaussE_gfn = PAMR_get_gfn("gaussE", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_gaussE_gfn = PAMR_get_gfn("gaussE", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_gaussB_gfn = PAMR_get_gfn("gaussB", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_gaussB_gfn = PAMR_get_gfn("gaussB", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_gaussB_gfn = PAMR_get_gfn("gaussB", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Xx_gfn = PAMR_get_gfn("Xx", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Xx_gfn = PAMR_get_gfn("Xx", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Xx_gfn = PAMR_get_gfn("Xx", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dxdX_gfn = PAMR_get_gfn("dxdX", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dxdX_gfn = PAMR_get_gfn("dxdX", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dxdX_gfn = PAMR_get_gfn("dxdX", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dxdX2_gfn = PAMR_get_gfn("dxdX2", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dxdX2_gfn = PAMR_get_gfn("dxdX2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dxdX2_gfn = PAMR_get_gfn("dxdX2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Yy_gfn = PAMR_get_gfn("Yy", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Yy_gfn = PAMR_get_gfn("Yy", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Yy_gfn = PAMR_get_gfn("Yy", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dydY_gfn = PAMR_get_gfn("dydY", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dydY_gfn = PAMR_get_gfn("dydY", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dydY_gfn = PAMR_get_gfn("dydY", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dydY2_gfn = PAMR_get_gfn("dydY2", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dydY2_gfn = PAMR_get_gfn("dydY2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dydY2_gfn = PAMR_get_gfn("dydY2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_Zz_gfn = PAMR_get_gfn("Zz", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_Zz_gfn = PAMR_get_gfn("Zz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_Zz_gfn = PAMR_get_gfn("Zz", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dzdZ_gfn = PAMR_get_gfn("dzdZ", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dzdZ_gfn = PAMR_get_gfn("dzdZ", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dzdZ_gfn = PAMR_get_gfn("dzdZ", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_dzdZ2_gfn = PAMR_get_gfn("dzdZ2", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_dzdZ2_gfn = PAMR_get_gfn("dzdZ2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_dzdZ2_gfn = PAMR_get_gfn("dzdZ2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_jac_gfn = PAMR_get_gfn("jac", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_jac_gfn = PAMR_get_gfn("jac", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_jac_gfn = PAMR_get_gfn("jac", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((cmask_gfn = PAMR_get_gfn("cmask", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((phi1_k1_gfn = PAMR_get_gfn("phi1_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1_k2_gfn = PAMR_get_gfn("phi1_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1_k3_gfn = PAMR_get_gfn("phi1_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1_k4_gfn = PAMR_get_gfn("phi1_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_k1_gfn = PAMR_get_gfn("phi2_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_k2_gfn = PAMR_get_gfn("phi2_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_k3_gfn = PAMR_get_gfn("phi2_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_k4_gfn = PAMR_get_gfn("phi2_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_k1_gfn = PAMR_get_gfn("pi1_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_k2_gfn = PAMR_get_gfn("pi1_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_k3_gfn = PAMR_get_gfn("pi1_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_k4_gfn = PAMR_get_gfn("pi1_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_k1_gfn = PAMR_get_gfn("pi2_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_k2_gfn = PAMR_get_gfn("pi2_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_k3_gfn = PAMR_get_gfn("pi2_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_k4_gfn = PAMR_get_gfn("pi2_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((At_k1_gfn = PAMR_get_gfn("At_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((At_k2_gfn = PAMR_get_gfn("At_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((At_k3_gfn = PAMR_get_gfn("At_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((At_k4_gfn = PAMR_get_gfn("At_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ax_k1_gfn = PAMR_get_gfn("Ax_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ax_k2_gfn = PAMR_get_gfn("Ax_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ax_k3_gfn = PAMR_get_gfn("Ax_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ax_k4_gfn = PAMR_get_gfn("Ax_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ay_k1_gfn = PAMR_get_gfn("Ay_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ay_k2_gfn = PAMR_get_gfn("Ay_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ay_k3_gfn = PAMR_get_gfn("Ay_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Ay_k4_gfn = PAMR_get_gfn("Ay_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Az_k1_gfn = PAMR_get_gfn("Az_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Az_k2_gfn = PAMR_get_gfn("Az_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Az_k3_gfn = PAMR_get_gfn("Az_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Az_k4_gfn = PAMR_get_gfn("Az_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bt_k1_gfn = PAMR_get_gfn("Bt_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bt_k2_gfn = PAMR_get_gfn("Bt_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bt_k3_gfn = PAMR_get_gfn("Bt_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bt_k4_gfn = PAMR_get_gfn("Bt_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bx_k1_gfn = PAMR_get_gfn("Bx_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bx_k2_gfn = PAMR_get_gfn("Bx_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bx_k3_gfn = PAMR_get_gfn("Bx_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bx_k4_gfn = PAMR_get_gfn("Bx_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((By_k1_gfn = PAMR_get_gfn("By_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((By_k2_gfn = PAMR_get_gfn("By_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((By_k3_gfn = PAMR_get_gfn("By_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((By_k4_gfn = PAMR_get_gfn("By_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bz_k1_gfn = PAMR_get_gfn("Bz_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bz_k2_gfn = PAMR_get_gfn("Bz_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bz_k3_gfn = PAMR_get_gfn("Bz_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Bz_k4_gfn = PAMR_get_gfn("Bz_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((chi_k1_gfn = PAMR_get_gfn("chi_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((chi_k2_gfn = PAMR_get_gfn("chi_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((chi_k3_gfn = PAMR_get_gfn("chi_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((chi_k4_gfn = PAMR_get_gfn("chi_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_k1_gfn = PAMR_get_gfn("xi_k1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_k2_gfn = PAMR_get_gfn("xi_k2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_k3_gfn = PAMR_get_gfn("xi_k3", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_k4_gfn = PAMR_get_gfn("xi_k4", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((nm1_V_gfn = PAMR_get_gfn("V", PAMR_AMRH, 3)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((n_V_gfn = PAMR_get_gfn("V", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((np1_V_gfn = PAMR_get_gfn("V", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if (mg_toggle == 1) {
    if ((V_gfn = PAMR_get_gfn("V", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((V_res_gfn = PAMR_get_gfn("V_res", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((V_lop_gfn = PAMR_get_gfn("V_lop", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((V_rhs_gfn = PAMR_get_gfn("V_rhs", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((mask_mg_gfn = PAMR_get_gfn("cmask", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((mask_gfn = PAMR_get_gfn("cmask", PAMR_AMRH, 1)) < 0)
      AMRD_stop("set_gnfs error", 0);

    if ((defect_gfn = PAMR_get_gfn("defect", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);

    if ((dxdX_gfn = PAMR_get_gfn("dxdX", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dxdX2_gfn = PAMR_get_gfn("dxdX2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dydY_gfn = PAMR_get_gfn("dydY", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dydY2_gfn = PAMR_get_gfn("dydY2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dzdZ_gfn = PAMR_get_gfn("dzdZ", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dzdZ2_gfn = PAMR_get_gfn("dzdZ2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((Xx_gfn = PAMR_get_gfn("Xx", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((Yy_gfn = PAMR_get_gfn("Yy", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((Zz_gfn = PAMR_get_gfn("Zz", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);

    if ((phi1_gfn = PAMR_get_gfn("phi1", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((phi2_gfn = PAMR_get_gfn("phi2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((pi1_gfn = PAMR_get_gfn("pi1", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((pi2_gfn = PAMR_get_gfn("pi2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((modphi_gfn = PAMR_get_gfn("modphi", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((At_gfn = PAMR_get_gfn("At", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
  }

  // compute L-infinity norm of each grid function
  g_norms = AMRD_get_global_norms();
}

/***********************************************************************
 * ldptr
 *
 * Initializes grid parameters and sets up pointers for grid functions.
 * The following operations are performed:
 * - Initializes the coordinate bounding box of the base grid
 * - Retrieves grid dimensions, rank, shape, and ghost width
 * - Calculates and stores the grid spacings dx, dy, dz
 * - Checks for physical boundaries and computes the grid size
 * - Initializes pointer arrays for coordinates and grid function data
 *
 **********************************************************************/
void ldptr(void) {
  real dx0[3];              // stores grid spacings dx, dy, dz
  real *x0[3];              // stores pointers to coordinates x, y, z
  real *gfs[PAMR_MAX_GFNS]; // stores pointers to grid function data

  static int first = 1; // only initialize grid functions once
  if (first) {
    first = 0;
    set_gfns();

    // initialize coordinate bounding box of the base grid
    // (e.g. [xmin, xmax, ymin, ymax, zmin, zmax])
    PAMR_get_global_bbox(base_bbox);
  }

  // retrieve properties of the grid
  PAMR_get_g_dim(&dim);
  PAMR_get_g_rank(&g_rank);
  PAMR_get_g_shape(shape);
  PAMR_get_g_bbox(bbox);
  PAMR_get_g_ghost_width(ghost_width);
  PAMR_get_g_level(&g_L);
  PAMR_get_dxdt(g_L, dx0, &dt); // fills dx0 and dt

  // compute grid spacings
  dx = dx0[0];
  dy = dx0[1];
  dz = dx0[2];

  // check for AMR and physical boundaries and compute grid size
  if ((bbox[0] - base_bbox[0]) < dx / 2)
    phys_bdy[0] = 1;
  else
    phys_bdy[0] = 0;
  if ((base_bbox[1] - bbox[1]) < dx / 2)
    phys_bdy[1] = 1;
  else
    phys_bdy[1] = 0;
  Nx = shape[0];
  size = Nx;
  Ny = 1;
  if ((bbox[2] - base_bbox[2]) < dy / 2)
    phys_bdy[2] = 1;
  else
    phys_bdy[2] = 0;
  if ((base_bbox[3] - bbox[3]) < dy / 2)
    phys_bdy[3] = 1;
  else
    phys_bdy[3] = 0;
  Ny = shape[1];
  size *= Ny;
  if ((bbox[4] - base_bbox[4]) < dz / 4)
    phys_bdy[4] = 1;
  else
    phys_bdy[4] = 0;
  if ((base_bbox[5] - bbox[5]) < dz / 4)
    phys_bdy[5] = 1;
  else
    phys_bdy[5] = 0;
  Nz = shape[2];
  size *= Nz;

  PAMR_get_g_x(x0); // retrieve pointers to arrays containing coordinates
  x = x0[0];
  y = x0[1];
  z = x0[2];

  // retrieve pointers to arrays containing grid function data
  PAMR_get_g_gfs(gfs);
  nm1_phi1 = gfs[nm1_phi1_gfn - 1];
  n_phi1 = gfs[n_phi1_gfn - 1];
  np1_phi1 = gfs[np1_phi1_gfn - 1];
  nm1_phi2 = gfs[nm1_phi2_gfn - 1];
  n_phi2 = gfs[n_phi2_gfn - 1];
  np1_phi2 = gfs[np1_phi2_gfn - 1];
  nm1_pi1 = gfs[nm1_pi1_gfn - 1];
  n_pi1 = gfs[n_pi1_gfn - 1];
  np1_pi1 = gfs[np1_pi1_gfn - 1];
  nm1_pi2 = gfs[nm1_pi2_gfn - 1];
  n_pi2 = gfs[n_pi2_gfn - 1];
  np1_pi2 = gfs[np1_pi2_gfn - 1];
  nm1_At = gfs[nm1_At_gfn - 1];
  n_At = gfs[n_At_gfn - 1];
  np1_At = gfs[np1_At_gfn - 1];
  nm1_Ax = gfs[nm1_Ax_gfn - 1];
  n_Ax = gfs[n_Ax_gfn - 1];
  np1_Ax = gfs[np1_Ax_gfn - 1];
  nm1_Ay = gfs[nm1_Ay_gfn - 1];
  n_Ay = gfs[n_Ay_gfn - 1];
  np1_Ay = gfs[np1_Ay_gfn - 1];
  nm1_Az = gfs[nm1_Az_gfn - 1];
  n_Az = gfs[n_Az_gfn - 1];
  np1_Az = gfs[np1_Az_gfn - 1];
  nm1_Bt = gfs[nm1_Bt_gfn - 1];
  n_Bt = gfs[n_Bt_gfn - 1];
  np1_Bt = gfs[np1_Bt_gfn - 1];
  nm1_Bx = gfs[nm1_Bx_gfn - 1];
  n_Bx = gfs[n_Bx_gfn - 1];
  np1_Bx = gfs[np1_Bx_gfn - 1];
  nm1_By = gfs[nm1_By_gfn - 1];
  n_By = gfs[n_By_gfn - 1];
  np1_By = gfs[np1_By_gfn - 1];
  nm1_Bz = gfs[nm1_Bz_gfn - 1];
  n_Bz = gfs[n_Bz_gfn - 1];
  np1_Bz = gfs[np1_Bz_gfn - 1];
  nm1_chi = gfs[nm1_chi_gfn - 1];
  n_chi = gfs[n_chi_gfn - 1];
  np1_chi = gfs[np1_chi_gfn - 1];
  nm1_xi = gfs[nm1_xi_gfn - 1];
  n_xi = gfs[n_xi_gfn - 1];
  np1_xi = gfs[np1_xi_gfn - 1];
  nm1_modphi = gfs[nm1_modphi_gfn - 1];
  n_modphi = gfs[n_modphi_gfn - 1];
  np1_modphi = gfs[np1_modphi_gfn - 1];
  nm1_qden = gfs[nm1_qden_gfn - 1];
  n_qden = gfs[n_qden_gfn - 1];
  np1_qden = gfs[np1_qden_gfn - 1];
  nm1_eden = gfs[nm1_eden_gfn - 1];
  n_eden = gfs[n_eden_gfn - 1];
  np1_eden = gfs[np1_eden_gfn - 1];
  nm1_jdenx = gfs[nm1_jdenx_gfn - 1];
  n_jdenx = gfs[n_jdenx_gfn - 1];
  np1_jdenx = gfs[np1_jdenx_gfn - 1];
  nm1_lorenz = gfs[nm1_lorenz_gfn - 1];
  n_lorenz = gfs[n_lorenz_gfn - 1];
  np1_lorenz = gfs[np1_lorenz_gfn - 1];
  nm1_elecx = gfs[nm1_elecx_gfn - 1];
  n_elecx = gfs[n_elecx_gfn - 1];
  np1_elecx = gfs[np1_elecx_gfn - 1];
  nm1_elecy = gfs[nm1_elecy_gfn - 1];
  n_elecy = gfs[n_elecy_gfn - 1];
  np1_elecy = gfs[np1_elecy_gfn - 1];
  nm1_elecz = gfs[nm1_elecz_gfn - 1];
  n_elecz = gfs[n_elecz_gfn - 1];
  np1_elecz = gfs[np1_elecz_gfn - 1];
  nm1_magx = gfs[nm1_magx_gfn - 1];
  n_magx = gfs[n_magx_gfn - 1];
  np1_magx = gfs[np1_magx_gfn - 1];
  nm1_magy = gfs[nm1_magy_gfn - 1];
  n_magy = gfs[n_magy_gfn - 1];
  np1_magy = gfs[np1_magy_gfn - 1];
  nm1_magz = gfs[nm1_magz_gfn - 1];
  n_magz = gfs[n_magz_gfn - 1];
  np1_magz = gfs[np1_magz_gfn - 1];
  nm1_eem = gfs[nm1_eem_gfn - 1];
  n_eem = gfs[n_eem_gfn - 1];
  np1_eem = gfs[np1_eem_gfn - 1];
  nm1_gaussE = gfs[nm1_gaussE_gfn - 1];
  n_gaussE = gfs[n_gaussE_gfn - 1];
  np1_gaussE = gfs[np1_gaussE_gfn - 1];
  nm1_gaussB = gfs[nm1_gaussB_gfn - 1];
  n_gaussB = gfs[n_gaussB_gfn - 1];
  np1_gaussB = gfs[np1_gaussB_gfn - 1];

  nm1_Xx = gfs[nm1_Xx_gfn - 1];
  n_Xx = gfs[n_Xx_gfn - 1];
  np1_Xx = gfs[np1_Xx_gfn - 1];
  nm1_dxdX = gfs[nm1_dxdX_gfn - 1];
  n_dxdX = gfs[n_dxdX_gfn - 1];
  np1_dxdX = gfs[np1_dxdX_gfn - 1];
  nm1_dxdX2 = gfs[nm1_dxdX2_gfn - 1];
  n_dxdX2 = gfs[n_dxdX2_gfn - 1];
  np1_dxdX2 = gfs[np1_dxdX2_gfn - 1];
  nm1_Yy = gfs[nm1_Yy_gfn - 1];
  n_Yy = gfs[n_Yy_gfn - 1];
  np1_Yy = gfs[np1_Yy_gfn - 1];
  nm1_dydY = gfs[nm1_dydY_gfn - 1];
  n_dydY = gfs[n_dydY_gfn - 1];
  np1_dydY = gfs[np1_dydY_gfn - 1];
  nm1_dydY2 = gfs[nm1_dydY2_gfn - 1];
  n_dydY2 = gfs[n_dydY2_gfn - 1];
  np1_dydY2 = gfs[np1_dydY2_gfn - 1];
  nm1_Zz = gfs[nm1_Zz_gfn - 1];
  n_Zz = gfs[n_Zz_gfn - 1];
  np1_Zz = gfs[np1_Zz_gfn - 1];
  nm1_dzdZ = gfs[nm1_dzdZ_gfn - 1];
  n_dzdZ = gfs[n_dzdZ_gfn - 1];
  np1_dzdZ = gfs[np1_dzdZ_gfn - 1];
  nm1_dzdZ2 = gfs[nm1_dzdZ2_gfn - 1];
  n_dzdZ2 = gfs[n_dzdZ2_gfn - 1];
  np1_dzdZ2 = gfs[np1_dzdZ2_gfn - 1];
  nm1_jac = gfs[nm1_jac_gfn - 1];
  n_jac = gfs[n_jac_gfn - 1];
  np1_jac = gfs[np1_jac_gfn - 1];

  cmask = gfs[cmask_gfn - 1];

  phi1_k1 = gfs[phi1_k1_gfn - 1];
  phi1_k2 = gfs[phi1_k2_gfn - 1];
  phi1_k3 = gfs[phi1_k3_gfn - 1];
  phi1_k4 = gfs[phi1_k4_gfn - 1];
  phi2_k1 = gfs[phi2_k1_gfn - 1];
  phi2_k2 = gfs[phi2_k2_gfn - 1];
  phi2_k3 = gfs[phi2_k3_gfn - 1];
  phi2_k4 = gfs[phi2_k4_gfn - 1];
  pi1_k1 = gfs[pi1_k1_gfn - 1];
  pi1_k2 = gfs[pi1_k2_gfn - 1];
  pi1_k3 = gfs[pi1_k3_gfn - 1];
  pi1_k4 = gfs[pi1_k4_gfn - 1];
  pi2_k1 = gfs[pi2_k1_gfn - 1];
  pi2_k2 = gfs[pi2_k2_gfn - 1];
  pi2_k3 = gfs[pi2_k3_gfn - 1];
  pi2_k4 = gfs[pi2_k4_gfn - 1];
  At_k1 = gfs[At_k1_gfn - 1];
  At_k2 = gfs[At_k2_gfn - 1];
  At_k3 = gfs[At_k3_gfn - 1];
  At_k4 = gfs[At_k4_gfn - 1];
  Ax_k1 = gfs[Ax_k1_gfn - 1];
  Ax_k2 = gfs[Ax_k2_gfn - 1];
  Ax_k3 = gfs[Ax_k3_gfn - 1];
  Ax_k4 = gfs[Ax_k4_gfn - 1];
  Ay_k1 = gfs[Ay_k1_gfn - 1];
  Ay_k2 = gfs[Ay_k2_gfn - 1];
  Ay_k3 = gfs[Ay_k3_gfn - 1];
  Ay_k4 = gfs[Ay_k4_gfn - 1];
  Az_k1 = gfs[Az_k1_gfn - 1];
  Az_k2 = gfs[Az_k2_gfn - 1];
  Az_k3 = gfs[Az_k3_gfn - 1];
  Az_k4 = gfs[Az_k4_gfn - 1];
  Bt_k1 = gfs[Bt_k1_gfn - 1];
  Bt_k2 = gfs[Bt_k2_gfn - 1];
  Bt_k3 = gfs[Bt_k3_gfn - 1];
  Bt_k4 = gfs[Bt_k4_gfn - 1];
  Bx_k1 = gfs[Bx_k1_gfn - 1];
  Bx_k2 = gfs[Bx_k2_gfn - 1];
  Bx_k3 = gfs[Bx_k3_gfn - 1];
  Bx_k4 = gfs[Bx_k4_gfn - 1];
  By_k1 = gfs[By_k1_gfn - 1];
  By_k2 = gfs[By_k2_gfn - 1];
  By_k3 = gfs[By_k3_gfn - 1];
  By_k4 = gfs[By_k4_gfn - 1];
  Bz_k1 = gfs[Bz_k1_gfn - 1];
  Bz_k2 = gfs[Bz_k2_gfn - 1];
  Bz_k3 = gfs[Bz_k3_gfn - 1];
  Bz_k4 = gfs[Bz_k4_gfn - 1];
  chi_k1 = gfs[chi_k1_gfn - 1];
  chi_k2 = gfs[chi_k2_gfn - 1];
  chi_k3 = gfs[chi_k3_gfn - 1];
  chi_k4 = gfs[chi_k4_gfn - 1];
  xi_k1 = gfs[xi_k1_gfn - 1];
  xi_k2 = gfs[xi_k2_gfn - 1];
  xi_k3 = gfs[xi_k3_gfn - 1];
  xi_k4 = gfs[xi_k4_gfn - 1];

  nm1_V = gfs[nm1_V_gfn - 1];
  n_V = gfs[n_V_gfn - 1];
  np1_V = gfs[np1_V_gfn - 1];

  if (mg_toggle == 1) {
    V = gfs[V_gfn - 1];

    mask_mg = gfs[mask_mg_gfn - 1];

    V_res = gfs[V_res_gfn - 1];
    V_rhs = gfs[V_rhs_gfn - 1];
    V_lop = gfs[V_lop_gfn - 1];

    defect = gfs[defect_gfn - 1];

    dxdX = gfs[dxdX_gfn - 1];
    dxdX2 = gfs[dxdX2_gfn - 1];
    dydY = gfs[dydY_gfn - 1];
    dydY2 = gfs[dydY2_gfn - 1];
    dzdZ = gfs[dzdZ_gfn - 1];
    dzdZ2 = gfs[dzdZ2_gfn - 1];
    Xx = gfs[Xx_gfn - 1];
    Yy = gfs[Yy_gfn - 1];
    Zz = gfs[Zz_gfn - 1];

    phi1 = gfs[phi1_gfn - 1];
    phi2 = gfs[phi2_gfn - 1];
    pi1 = gfs[pi1_gfn - 1];
    pi2 = gfs[pi2_gfn - 1];
    modphi = gfs[modphi_gfn - 1];
    At = gfs[At_gfn - 1];
  }

  // recompute grid spacing
  x = x0[0];
  dx = x[1] - x[0];
  y = x0[1];
  dy = y[1] - y[0];
  z = x0[2];
  dz = z[1] - z[0];
}

/***********************************************************************
 * const_f
 *
 * Sets a grid function to a constant value.
 *
 ***********************************************************************/
void const_f(real *f, real c) {
  int i;

  for (i = 0; i < Nx * Ny * Nz; i++)
    f[i] = c;
}

/***********************************************************************
 * zero
 *
 * Sets a grid function to zero.
 *
 ***********************************************************************/
void zero(real *f) { const_f(f, 0); }

/***********************************************************************
 * l2norm_calc
 *
 * Computes the L2-norm of a grid function.
 *
 ***********************************************************************/
real l2norm_calc(real *f) {
  int i;
  int sum = 0;
  real norm = 0;

  for (i = 0; i < size; i++) {
    sum++;
    norm += f[i] * f[i];
  }

  // if sum=0 somehow, set sum=1 to avoid NaN
  if (!sum)
    sum = 1;

  // if norm/sum is very small, return 0
  if (norm / sum < 1.0E-50)
    return 0;

  return (sqrt(norm / sum)); // L2-norm
}

/***********************************************************************
 * elapsed_time
 *
 * Maintains and reports elapsed wall-clock time.
 *
 ***********************************************************************/
void elapsed_time(void) {
  static int first = 1;
  struct timeb t;
  static real msinit;
  real mscurr, mselapsed;

  ftime(&t);
  mscurr = 1000.0 * t.time + t.millitm;
  if (first) {
    msinit = mscurr;
    first = 0;
  }
  mselapsed = mscurr - msinit;
  printf("\nelapsed_time: Seconds since initial call: %12.3f\n",
         mselapsed / 1000.0);
  printf("              Minutes since initial call: %12.3f\n",
         mselapsed / 60000.0);
  printf("              Hours   since initial call: %12.3f\n",
         mselapsed / 3600000.0);
}

/***********************************************************************
 * qball_id
 *
 * Determines the initialization behaviour for the hierarchy. If the
 * function is called by the process with rank 0, it also calculates and
 * logs the elapsed time of the calculation.
 *
 ***********************************************************************/
int qball_id(void) {
  if (my_rank == 0)
    elapsed_time();

  return 0;
}

/***********************************************************************
 * qball_var_pre_init
 *
 * Configures custom behaviour before parameters are read from the input
 * parameter files and the base hierarchy is initialized. In this case,
 * we print a "ready" message on each host to the log file and
 * initialize multigrid functionality.
 *
 ***********************************************************************/
void qball_var_pre_init(char *pfile) {
  int *base_shape = NULL;
  int dim;
  char hostname[256]; // name of the processor

  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready\n", getpid(), hostname);
  fflush(stdout);
  sleep(1);

  // uncomment if you want to attach debugger at the start of execution
  // volatile int i = 0;
  // while (0 == i) {sleep(5);} // sleeping won't use 100% CPU

  // initialize multigrid functionality
  mg_toggle = 0;
  AMRD_int_param(pfile, "mg_toggle", &mg_toggle, 1);
  PAMR_get_g_dim(&dim);
  AMRD_int_param_v(pfile, "base_shape", &base_shape, &dim);
  MGshape = base_shape[0]; // assumes equal points in all dimensions
}

/***********************************************************************
 * qball_var_post_init
 *
 * Configures custom behaviour after parameters are read from the input
 * parameter files and the base hierarchy is initialized. In this case,
 * we read values from the input parameter file and simultaneously print
 * them to a log file.
 *
 ***********************************************************************/
void qball_var_post_init(char *pfile) {
  if (my_rank == 0) {
    printf("==================================================================="
           "\n");
    printf("Reading parameters:\n\n");
  }

  // read custom parameters from parameter file (RNPL format)
  AMRD_real_param(pfile, "e", &e, 1);
  AMRD_real_param(pfile, "m", &m, 1);
  AMRD_real_param(pfile, "g", &g, 1);
  AMRD_real_param(pfile, "h", &h, 1);
  AMRD_real_param(pfile, "beta", &beta, 1);
  AMRD_real_param(pfile, "mu", &mu, 1);
  AMRD_real_param(pfile, "w", &w, 1);
  AMRD_real_param(pfile, "w2", &w2, 1);
  AMRD_real_param(pfile, "cx", &cx, 1);
  AMRD_real_param(pfile, "cy", &cy, 1);
  AMRD_real_param(pfile, "cz", &cz, 1);
  AMRD_real_param(pfile, "cx1", &cx1, 1);
  AMRD_real_param(pfile, "cx2", &cx2, 1);
  AMRD_real_param(pfile, "cy1", &cy1, 1);
  AMRD_real_param(pfile, "cy2", &cy2, 1);
  AMRD_real_param(pfile, "cz1", &cz1, 1);
  AMRD_real_param(pfile, "cz2", &cz2, 1);
  AMRD_real_param(pfile, "v", &v, 1);
  AMRD_real_param(pfile, "phase", &phase, 1);
  AMRD_int_param(pfile, "anti", &anti, 1);
  AMRD_int_param(pfile, "V_type", &V_type, 1);
  AMRD_int_param(pfile, "ini_type", &ini_type, 1);
  AMRD_real_param(pfile, "c", &c, 1);
  AMRD_real_param(pfile, "d", &d, 1);
  AMRD_real_param(pfile, "cc", &cc, 1);
  AMRD_real_param(pfile, "amp", &amp, 1);
  AMRD_real_param(pfile, "chiX0", &chiX0, 1);
  AMRD_real_param(pfile, "chiY0", &chiY0, 1);
  AMRD_real_param(pfile, "chiZ0", &chiY0, 1);
  AMRD_real_param(pfile, "wwX", &wwX, 1);
  AMRD_real_param(pfile, "wwY", &wwY, 1);
  AMRD_real_param(pfile, "wwZ", &wwZ, 1);
  AMRD_real_param(pfile, "r0", &r0, 1);
  AMRD_real_param(pfile, "delta", &delta, 1);
  AMRD_real_param(pfile, "myzero", &myzero, 1);

  if (my_rank == 0)
    printf("==================================================================="
           "\n");
}

/***********************************************************************
 * qball_AMRH_var_clear
 *
 * Sets all variables in the AMR hierarchy to their "zero" values.
 *
 ***********************************************************************/
void qball_AMRH_var_clear(void) {
  ldptr();

  zero(nm1_phi1);
  zero(n_phi1);
  zero(np1_phi1);
  zero(nm1_phi2);
  zero(n_phi2);
  zero(np1_phi2);
  zero(nm1_pi1);
  zero(n_pi1);
  zero(np1_pi1);
  zero(nm1_pi2);
  zero(n_pi2);
  zero(np1_pi2);
  zero(nm1_At);
  zero(n_At);
  zero(np1_At);
  zero(nm1_Ax);
  zero(n_Ax);
  zero(np1_Ax);
  zero(nm1_Ay);
  zero(n_Ay);
  zero(np1_Ay);
  zero(nm1_Az);
  zero(n_Az);
  zero(np1_Az);
  zero(nm1_Bt);
  zero(n_Bt);
  zero(np1_Bt);
  zero(nm1_Bx);
  zero(n_Bx);
  zero(np1_Bx);
  zero(nm1_By);
  zero(n_By);
  zero(np1_By);
  zero(nm1_Bz);
  zero(n_Bz);
  zero(np1_Bz);
  zero(nm1_chi);
  zero(n_chi);
  zero(np1_chi);
  zero(nm1_xi);
  zero(n_xi);
  zero(np1_xi);
  zero(nm1_modphi);
  zero(n_modphi);
  zero(np1_modphi);
  zero(nm1_qden);
  zero(n_qden);
  zero(np1_qden);
  zero(nm1_eden);
  zero(n_eden);
  zero(np1_eden);
  zero(nm1_jdenx);
  zero(n_jdenx);
  zero(np1_jdenx);
  zero(nm1_lorenz);
  zero(n_lorenz);
  zero(np1_lorenz);
  zero(nm1_elecx);
  zero(n_elecx);
  zero(np1_elecx);
  zero(nm1_elecy);
  zero(n_elecy);
  zero(np1_elecy);
  zero(nm1_elecz);
  zero(n_elecz);
  zero(np1_elecz);
  zero(nm1_magx);
  zero(n_magx);
  zero(np1_magx);
  zero(nm1_magy);
  zero(n_magy);
  zero(np1_magy);
  zero(nm1_magz);
  zero(n_magz);
  zero(np1_magz);
  zero(nm1_eem);
  zero(n_eem);
  zero(np1_eem);
  zero(nm1_gaussE);
  zero(n_gaussE);
  zero(np1_gaussE);
  zero(nm1_gaussB);
  zero(n_gaussB);
  zero(np1_gaussB);

  zero(nm1_Xx);
  zero(n_Xx);
  zero(np1_Xx);
  zero(nm1_dxdX);
  zero(n_dxdX);
  zero(np1_dxdX);
  zero(nm1_dxdX2);
  zero(n_dxdX2);
  zero(np1_dxdX2);
  zero(nm1_Yy);
  zero(n_Yy);
  zero(np1_Yy);
  zero(nm1_dydY);
  zero(n_dydY);
  zero(np1_dydY);
  zero(nm1_dydY2);
  zero(n_dydY2);
  zero(np1_dydY2);
  zero(nm1_Zz);
  zero(n_Zz);
  zero(np1_Zz);
  zero(nm1_dzdZ);
  zero(n_dzdZ);
  zero(np1_dzdZ);
  zero(nm1_dzdZ2);
  zero(n_dzdZ2);
  zero(np1_dzdZ2);
  zero(nm1_jac);
  zero(n_jac);
  zero(np1_jac);

  zero(cmask);

  zero(nm1_V);
  zero(n_V);
  zero(np1_V);
}

/***********************************************************************
 * qball_free_data
 *
 * Generates initial data for all relevant grid functions.
 *
 ***********************************************************************/
void qball_free_data(void) {
  ldptr();

  init_xx_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Xx);
  init_dxdx_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dxdX);
  init_dxdx2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dxdX2);
  init_yy_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Yy);
  init_dydy_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dydY);
  init_dydy2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dydY2);
  init_zz_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Zz);
  init_dzdz_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dzdZ);
  init_dzdz2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dzdZ2);
  init_jac_(x, y, z, &Nx, &Ny, &Nz, &c, &d, n_jac);

  init_qball_(n_Xx, n_dxdX, n_dxdX2, n_Yy, n_dydY, n_dydY2, n_Zz, n_dzdZ,
              n_dzdZ2, np1_phi1, n_phi1, np1_phi2, n_phi2, np1_pi1, n_pi1,
              np1_pi2, n_pi2, np1_At, n_At, np1_Ax, n_Ax, np1_Ay, n_Ay, np1_Az,
              n_Az, np1_Bt, n_Bt, np1_Bx, n_Bx, np1_By, n_By, np1_Bz, n_Bz, &Nx,
              &Ny, &Nz, x, y, z, &dx, &dy, &dz, &bbox[1], &bbox[0], &bbox[3],
              &bbox[2], &bbox[5], &bbox[4], &cx, &cy, &cz, &cx1, &cx2, &cy1,
              &cy2, &cz1, &cz2, &w, &v, &phase, &anti, &ini_type, &e);

  init_modphi_(n_phi1, n_phi2, &Nx, &Ny, &Nz, n_modphi);
  init_qden_(n_At, n_modphi, n_phi1, n_phi2, n_pi1, n_pi2, &Nx, &Ny, &Nz, &e, n_qden);

  // initialize energy differently depending on the scalar potential
  if (V_type == 0) { 
    init_eden_poly_(n_At, n_Ax, n_Ay, n_Az, n_Bx, n_By, n_Bz, n_dxdX, n_dydY,
                    n_dzdZ, n_phi1, n_phi2, n_pi1, n_pi2, &Nx, &Ny, &Nz, &e, &g,
                    &h, &dx, &dy, &dz, &m, n_eden);
  } else {
    init_eden_log_(n_At, n_Ax, n_Ay, n_Az, n_Bx, n_By, n_Bz, n_dxdX, n_dydY,
                   n_dzdZ, n_phi1, n_phi2, n_pi1, n_pi2, &Nx, &Ny, &Nz, &beta,
                   &e, &dx, &dy, &dz, &mu, n_eden);
  }
  init_jdenx_(n_At, n_Ax, n_Ay, n_Az, n_Bx, n_By, n_Bz, n_dxdX, n_dydY, n_dzdZ,
              n_phi1, n_phi2, n_pi1, n_pi2, y, z, &Nx, &Ny, &Nz, &c, &d, &e,
              &dx, &dy, &dz, n_jdenx);

  init_chi_(x, y, z, &Nx, &Ny, &Nz, &amp, &c, &chiX0, &chiY0, &chiZ0, &d,
            &delta, &r0, &wwX, &wwY, &wwZ, n_chi);
  init_xi_(n_chi, n_dxdX, n_dydY, n_dzdZ, x, y, z, &Nx, &Ny, &Nz, &c, &d, &dx,
           &dy, &dz, n_xi);

  // if multigrid data file exists, read it and evolve (overwriting previous data)
  // if multigrid data file doesn't exist, create it
  if (access(mg_data, F_OK) == 0) {
    init_qball2_(n_V, x, y, z, &Nx, &Ny, &Nz, &dx, &dy, &dz, n_Xx, n_Yy, n_Zz,
                 n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_At, n_Ax,
                 n_Ay, n_Az, n_Bt, n_Bx, n_By, n_Bz, &base_bbox[0],
                 &base_bbox[1], &base_bbox[2], &base_bbox[3], &base_bbox[4],
                 &base_bbox[5]);
  } else if (mg_toggle == 1)
  {
    initguess_(n_V, x, y, z, &Nx, &Ny, &Nz, n_Xx, n_Yy, n_Zz, n_At);
  }

  init_lorenz_(n_Ax, n_Ay, n_Az, n_Bt, n_dxdX, n_dydY, n_dzdZ, &Nx, &Ny, &Nz,
               &dx, &dy, &dz, n_lorenz);
  init_elecx_(n_At, n_Bx, n_dxdX, &Nx, &Ny, &Nz, &dx, n_elecx);
  init_elecy_(n_At, n_By, n_dydY, &Nx, &Ny, &Nz, &dy, n_elecy);
  init_elecz_(n_At, n_Bz, n_dzdZ, &Nx, &Ny, &Nz, &dz, n_elecz);
  init_magx_(n_Ay, n_Az, n_dydY, n_dzdZ, &Nx, &Ny, &Nz, &dy, &dz, n_magx);
  init_magy_(n_Ax, n_Az, n_dxdX, n_dzdZ, &Nx, &Ny, &Nz, &dx, &dz, n_magy);
  init_magz_(n_Ax, n_Ay, n_dxdX, n_dydY, &Nx, &Ny, &Nz, &dx, &dy, n_magz);
  init_eem_(n_elecx, n_elecy, n_elecz, n_magx, n_magy, n_magz, &Nx, &Ny, &Nz, n_eem);
  init_gausse_(n_At, n_Bx, n_By, n_Bz, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ,
               n_dzdZ2, n_qden, &Nx, &Ny, &Nz, &e, &dx, &dy, &dz, n_gaussE);
  init_gaussb_(n_dxdX, n_dydY, n_dzdZ, n_magx, n_magy, n_magz, &Nx, &Ny, &Nz,
               &dx, &dy, &dz, n_gaussB);
}

/***********************************************************************
 * qball_t0_cnst_data
 *
 * Initializes constraint data after each multigrid iteration.
 * Currently unused.
 *
 **********************************************************************/
void qball_t0_cnst_data(void) { return; }

/***********************************************************************
 * qball_pre_io_calc
 *
 * Performs calculations prior to saving information to disk.
 * Currently unused.
 *
 **********************************************************************/
void qball_pre_io_calc(void) { return; }

/***********************************************************************
 * qball_evo_residual
 *
 * Computes and returns a norm of the evolution residual.
 * Currently unused.
 *
 **********************************************************************/
real qball_evo_residual(void) { return 0.0; }

/***********************************************************************
 * qball_evolve
 *
 * Performs 1 iteration of the evolution equations.
 *
 **********************************************************************/
void qball_evolve(int iter, int *ifc_mask) {

  // if multigrid data file doesn't exist but mg_toggle turned on, write it and exit
  if ((access(mg_data, F_OK) != 0) && (mg_toggle == 1)) {
    w2f_(np1_V, &Nx, &Ny, &Nz); // write to binary file
    AMRD_stop("Multigrid processing complete.", 0);
  }

  ldptr();

  // update grid functions using classic RK4
  // this requires evo_max_iter = evo_min_iter = 4 
  //               t_interp_substeps  := [0 0.5 0.5 1] 
  //               diss_use_6th_order := 1
  switch (iter) {
  case 1:
    // re-initialize the coordinate variables on any AMR level required
    init_xx_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Xx);
    init_dxdx_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dxdX);
    init_dxdx2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dxdX2);
    init_yy_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Yy);
    init_dydy_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dydY);
    init_dydy2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dydY2);
    init_zz_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_Zz);
    init_dzdz_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dzdZ);
    init_dzdz2_(x, y, z, &Nx, &Ny, &Nz, &c, &d, &myzero, n_dzdZ2);
    init_jac_(x, y, z, &Nx, &Ny, &Nz, &c, &d, n_jac);

    // all grid functions should be at time level "n" to determine k1
    u_phi1_rk_(n_pi1, n_phi1, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi1_k1);
    u_phi2_rk_(n_pi2, n_phi2, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi2_k1);
    if (V_type == 0) {
      u_pi1_poly_rk_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                     n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi2, n_pi1, x,
                     y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h, &dx, &dy, &dz, &m,
                     &myzero, phys_bdy, pi1_k1);
      u_pi2_poly_rk_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                     n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi1, n_pi2, x,
                     y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h, &dx, &dy, &dz, &m,
                     &myzero, phys_bdy, pi2_k1);
    } else {
      u_pi1_log_rk_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                    n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi2, n_pi1, x,
                    y, z, &Nx, &Ny, &Nz, &beta, &cc, &e, &dx, &dy, &dz, &mu,
                    &myzero, phys_bdy, pi1_k1);
      u_pi2_log_rk_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                    n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi1, n_pi2, x,
                    y, z, &Nx, &Ny, &Nz, &beta, &cc, &e, &dx, &dy, &dz, &mu,
                    &myzero, phys_bdy, pi2_k1);
    }
    u_at_rk_(n_Bt, n_At, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, At_k1);
    u_ax_rk_(n_Bx, n_Ax, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ax_k1);
    u_ay_rk_(n_By, n_Ay, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ay_k1);
    u_az_rk_(n_Bz, n_Az, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Az_k1);
    u_bt_rk_(n_At, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
             n_phi2, n_pi1, n_pi2, n_Bt, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bt_k1);
    u_bx_rk_(n_Ax, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
             n_phi2, n_Bx, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy, &dz, &myzero,
             phys_bdy, Bx_k1);
    u_by_rk_(n_Ay, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
             n_phi2, n_By, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy, &dz, &myzero,
             phys_bdy, By_k1);
    u_bz_rk_(n_Az, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
             n_phi2, n_Bz, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy, &dz, &myzero,
             phys_bdy, Bz_k1);
    u_chi_rk_(n_xi, n_chi, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, chi_k1);
    u_xi_rk_(n_chi, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
             n_phi2, n_xi, x, y, z, &Nx, &Ny, &Nz, &cc, &dx, &dy, &dz, &myzero,
             phys_bdy, xi_k1);
    break;

  case 2: // all grid functions (except coordinates) should be at time level "np1" to determine k2
    u_phi1_rk_(np1_pi1, np1_phi1, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi1_k2);
    u_phi2_rk_(np1_pi2, np1_phi2, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi2_k2);
    if (V_type == 0) {
      u_pi1_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi1_k2);
      u_pi2_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi2_k2);
    } else {
      u_pi1_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi1_k2);
      u_pi2_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi2_k2);
    }
    u_at_rk_(np1_Bt, np1_At, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, At_k2);
    u_ax_rk_(np1_Bx, np1_Ax, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ax_k2);
    u_ay_rk_(np1_By, np1_Ay, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ay_k2);
    u_az_rk_(np1_Bz, np1_Az, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Az_k2);
    u_bt_rk_(np1_At, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_pi1, np1_pi2, np1_Bt, x, y, z, &Nx, &Ny,
             &Nz, &e, &dx, &dy, &dz, &myzero, phys_bdy, Bt_k2);
    u_bx_rk_(np1_Ax, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bx, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bx_k2);
    u_by_rk_(np1_Ay, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_By, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, By_k2);
    u_bz_rk_(np1_Az, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bz, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bz_k2);
    u_chi_rk_(np1_xi, np1_chi, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
              chi_k2);
    u_xi_rk_(np1_chi, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_xi, x, y, z, &Nx, &Ny, &Nz, &cc, &dx, &dy,
             &dz, &myzero, phys_bdy, xi_k2);
    break;

  case 3: // all grid functions (except coordinates) should be at time level "np1" to determine k3
    u_phi1_rk_(np1_pi1, np1_phi1, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi1_k3);
    u_phi2_rk_(np1_pi2, np1_phi2, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi2_k3);
    if (V_type == 0) {
      u_pi1_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi1_k3);
      u_pi2_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi2_k3);
    } else {
      u_pi1_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi1_k3);
      u_pi2_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi2_k3);
    }
    u_at_rk_(np1_Bt, np1_At, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, At_k3);
    u_ax_rk_(np1_Bx, np1_Ax, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ax_k3);
    u_ay_rk_(np1_By, np1_Ay, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ay_k3);
    u_az_rk_(np1_Bz, np1_Az, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Az_k3);
    u_bt_rk_(np1_At, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_pi1, np1_pi2, np1_Bt, x, y, z, &Nx, &Ny,
             &Nz, &e, &dx, &dy, &dz, &myzero, phys_bdy, Bt_k3);
    u_bx_rk_(np1_Ax, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bx, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bx_k3);
    u_by_rk_(np1_Ay, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_By, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, By_k3);
    u_bz_rk_(np1_Az, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bz, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bz_k3);
    u_chi_rk_(np1_xi, np1_chi, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
              chi_k3);
    u_xi_rk_(np1_chi, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_xi, x, y, z, &Nx, &Ny, &Nz, &cc, &dx, &dy,
             &dz, &myzero, phys_bdy, xi_k3);
    break;

  case 4: // all grid functions (except coordinates) should be at time level "np1" to determine k4
    u_phi1_rk_(np1_pi1, np1_phi1, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi1_k4);
    u_phi2_rk_(np1_pi2, np1_phi2, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
               phi2_k4);
    if (V_type == 0) {
      u_pi1_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi1_k4);
      u_pi2_poly_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                     n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                     np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &cc, &e, &g, &h,
                     &dx, &dy, &dz, &m, &myzero, phys_bdy, pi2_k4);
    } else {
      u_pi1_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi2, np1_pi1, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi1_k4);
      u_pi2_log_rk_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_chi, n_dxdX, n_dxdX2,
                    n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, np1_phi1, np1_phi2,
                    np1_pi1, np1_pi2, x, y, z, &Nx, &Ny, &Nz, &beta, &cc, &e,
                    &dx, &dy, &dz, &mu, &myzero, phys_bdy, pi2_k4);
    }
    u_at_rk_(np1_Bt, np1_At, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, At_k4);
    u_ax_rk_(np1_Bx, np1_Ax, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ax_k4);
    u_ay_rk_(np1_By, np1_Ay, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Ay_k4);
    u_az_rk_(np1_Bz, np1_Az, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy, Az_k4);
    u_bt_rk_(np1_At, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_pi1, np1_pi2, np1_Bt, x, y, z, &Nx, &Ny,
             &Nz, &e, &dx, &dy, &dz, &myzero, phys_bdy, Bt_k4);
    u_bx_rk_(np1_Ax, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bx, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bx_k4);
    u_by_rk_(np1_Ay, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_By, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, By_k4);
    u_bz_rk_(np1_Az, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_Bz, x, y, z, &Nx, &Ny, &Nz, &e, &dx, &dy,
             &dz, &myzero, phys_bdy, Bz_k4);
    u_chi_rk_(np1_xi, np1_chi, x, y, z, &Nx, &Ny, &Nz, &myzero, phys_bdy,
              chi_k4);
    u_xi_rk_(np1_chi, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2,
             np1_phi1, np1_phi2, np1_xi, x, y, z, &Nx, &Ny, &Nz, &cc, &dx, &dy,
             &dz, &myzero, phys_bdy, xi_k4);
    break;
  }

  switch (iter) {
  case 1: // we use time level "np1" as temporary storage for intermediate steps
    for (int k = 1; k < Nz - 1; k++) {
      for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
          int idx = (i) + (j * Nx) + (k * Nx * Ny);
          np1_phi1[idx] = n_phi1[idx] + phi1_k1[idx] * dt / 2.0;
          np1_phi2[idx] = n_phi2[idx] + phi2_k1[idx] * dt / 2.0;
          np1_pi1[idx] = n_pi1[idx] + pi1_k1[idx] * dt / 2.0;
          np1_pi2[idx] = n_pi2[idx] + pi2_k1[idx] * dt / 2.0;
          np1_At[idx] = n_At[idx] + At_k1[idx] * dt / 2.0;
          np1_Ax[idx] = n_Ax[idx] + Ax_k1[idx] * dt / 2.0;
          np1_Ay[idx] = n_Ay[idx] + Ay_k1[idx] * dt / 2.0;
          np1_Az[idx] = n_Az[idx] + Az_k1[idx] * dt / 2.0;
          np1_Bt[idx] = n_Bt[idx] + Bt_k1[idx] * dt / 2.0;
          np1_Bx[idx] = n_Bx[idx] + Bx_k1[idx] * dt / 2.0;
          np1_By[idx] = n_By[idx] + By_k1[idx] * dt / 2.0;
          np1_Bz[idx] = n_Bz[idx] + Bz_k1[idx] * dt / 2.0;
          np1_chi[idx] = n_chi[idx] + chi_k1[idx] * dt / 2.0;
          np1_xi[idx] = n_xi[idx] + xi_k1[idx] * dt / 2.0;
        }
      }
    }
    break;

  case 2:
    for (int k = 1; k < Nz - 1; k++) {
      for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
          int idx = (i) + (j * Nx) + (k * Nx * Ny);
          np1_phi1[idx] = n_phi1[idx] + phi1_k2[idx] * dt / 2.0;
          np1_phi2[idx] = n_phi2[idx] + phi2_k2[idx] * dt / 2.0;
          np1_pi1[idx] = n_pi1[idx] + pi1_k2[idx] * dt / 2.0;
          np1_pi2[idx] = n_pi2[idx] + pi2_k2[idx] * dt / 2.0;
          np1_At[idx] = n_At[idx] + At_k2[idx] * dt / 2.0;
          np1_Ax[idx] = n_Ax[idx] + Ax_k2[idx] * dt / 2.0;
          np1_Ay[idx] = n_Ay[idx] + Ay_k2[idx] * dt / 2.0;
          np1_Az[idx] = n_Az[idx] + Az_k2[idx] * dt / 2.0;
          np1_Bt[idx] = n_Bt[idx] + Bt_k2[idx] * dt / 2.0;
          np1_Bx[idx] = n_Bx[idx] + Bx_k2[idx] * dt / 2.0;
          np1_By[idx] = n_By[idx] + By_k2[idx] * dt / 2.0;
          np1_Bz[idx] = n_Bz[idx] + Bz_k2[idx] * dt / 2.0;
          np1_chi[idx] = n_chi[idx] + chi_k2[idx] * dt / 2.0;
          np1_xi[idx] = n_xi[idx] + xi_k2[idx] * dt / 2.0;
        }
      }
    }
    break;

  case 3:
    for (int k = 1; k < Nz - 1; k++) {
      for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
          int idx = (i) + (j * Nx) + (k * Nx * Ny);
          np1_phi1[idx] = n_phi1[idx] + phi1_k3[idx] * dt;
          np1_phi2[idx] = n_phi2[idx] + phi2_k3[idx] * dt;
          np1_pi1[idx] = n_pi1[idx] + pi1_k3[idx] * dt;
          np1_pi2[idx] = n_pi2[idx] + pi2_k3[idx] * dt;
          np1_At[idx] = n_At[idx] + At_k3[idx] * dt;
          np1_Ax[idx] = n_Ax[idx] + Ax_k3[idx] * dt;
          np1_Ay[idx] = n_Ay[idx] + Ay_k3[idx] * dt;
          np1_Az[idx] = n_Az[idx] + Az_k3[idx] * dt;
          np1_Bt[idx] = n_Bt[idx] + Bt_k3[idx] * dt;
          np1_Bx[idx] = n_Bx[idx] + Bx_k3[idx] * dt;
          np1_By[idx] = n_By[idx] + By_k3[idx] * dt;
          np1_Bz[idx] = n_Bz[idx] + Bz_k3[idx] * dt;
          np1_chi[idx] = n_chi[idx] + chi_k3[idx] * dt;
          np1_xi[idx] = n_xi[idx] + xi_k3[idx] * dt;
        }
      }
    }
    break;

  case 4: // perform the full RK4 update
    for (int k = 1; k < Nz - 1; k++) {
      for (int j = 1; j < Ny - 1; j++) {
        for (int i = 1; i < Nx - 1; i++) {
          int idx = (i) + (j * Nx) + (k * Nx * Ny);
          np1_phi1[idx] = n_phi1[idx] +
                          (phi1_k1[idx] + 2.0 * (phi1_k2[idx] + phi1_k3[idx]) +
                           phi1_k4[idx]) *
                              dt / 6.0;
          np1_phi2[idx] = n_phi2[idx] +
                          (phi2_k1[idx] + 2.0 * (phi2_k2[idx] + phi2_k3[idx]) +
                           phi2_k4[idx]) *
                              dt / 6.0;
          np1_pi1[idx] =
              n_pi1[idx] +
              (pi1_k1[idx] + 2.0 * (pi1_k2[idx] + pi1_k3[idx]) + pi1_k4[idx]) *
                  dt / 6.0;
          np1_pi2[idx] =
              n_pi2[idx] +
              (pi2_k1[idx] + 2.0 * (pi2_k2[idx] + pi2_k3[idx]) + pi2_k4[idx]) *
                  dt / 6.0;
          np1_At[idx] =
              n_At[idx] +
              (At_k1[idx] + 2.0 * (At_k2[idx] + At_k3[idx]) + At_k4[idx]) * dt /
                  6.0;
          np1_Ax[idx] =
              n_Ax[idx] +
              (Ax_k1[idx] + 2.0 * (Ax_k2[idx] + Ax_k3[idx]) + Ax_k4[idx]) * dt /
                  6.0;
          np1_Ay[idx] =
              n_Ay[idx] +
              (Ay_k1[idx] + 2.0 * (Ay_k2[idx] + Ay_k3[idx]) + Ay_k4[idx]) * dt /
                  6.0;
          np1_Az[idx] =
              n_Az[idx] +
              (Az_k1[idx] + 2.0 * (Az_k2[idx] + Az_k3[idx]) + Az_k4[idx]) * dt /
                  6.0;
          np1_Bt[idx] =
              n_Bt[idx] +
              (Bt_k1[idx] + 2.0 * (Bt_k2[idx] + Bt_k3[idx]) + Bt_k4[idx]) * dt /
                  6.0;
          np1_Bx[idx] =
              n_Bx[idx] +
              (Bx_k1[idx] + 2.0 * (Bx_k2[idx] + Bx_k3[idx]) + Bx_k4[idx]) * dt /
                  6.0;
          np1_By[idx] =
              n_By[idx] +
              (By_k1[idx] + 2.0 * (By_k2[idx] + By_k3[idx]) + By_k4[idx]) * dt /
                  6.0;
          np1_Bz[idx] =
              n_Bz[idx] +
              (Bz_k1[idx] + 2.0 * (Bz_k2[idx] + Bz_k3[idx]) + Bz_k4[idx]) * dt /
                  6.0;
          np1_chi[idx] =
              n_chi[idx] +
              (chi_k1[idx] + 2.0 * (chi_k2[idx] + chi_k3[idx]) + chi_k4[idx]) *
                  dt / 6.0;
          np1_xi[idx] =
              n_xi[idx] +
              (xi_k1[idx] + 2.0 * (xi_k2[idx] + xi_k3[idx]) + xi_k4[idx]) * dt /
                  6.0;
        }
      }
    }

    // update auxiliary grid functions
    u_modphi_(np1_phi1, np1_phi2, &Nx, &Ny, &Nz, phys_bdy, np1_modphi);
    u_qden_(np1_At, np1_modphi, np1_phi1, np1_phi2, np1_pi1, np1_pi2, &Nx, &Ny,
            &Nz, &e, phys_bdy, np1_qden);
    if (V_type == 0) {
      u_eden_poly_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_Bx, np1_By, np1_Bz,
                   np1_dxdX, np1_dydY, np1_dzdZ, np1_phi1, np1_phi2, np1_pi1,
                   np1_pi2, &Nx, &Ny, &Nz, &e, &g, &h, &dx, &dy, &dz, &m,
                   phys_bdy, np1_eden);
    } else {
      u_eden_log_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_Bx, np1_By, np1_Bz,
                  np1_dxdX, np1_dydY, np1_dzdZ, np1_phi1, np1_phi2, np1_pi1,
                  np1_pi2, &Nx, &Ny, &Nz, &beta, &e, &dx, &dy, &dz, &mu,
                  phys_bdy, np1_eden);
    }
    u_jdenx_(np1_At, np1_Ax, np1_Ay, np1_Az, np1_Bx, np1_By, np1_Bz, np1_dxdX,
             np1_dydY, np1_dzdZ, np1_phi1, np1_phi2, np1_pi1, np1_pi2, y, z,
             &Nx, &Ny, &Nz, &c, &d, &e, &dx, &dy, &dz, phys_bdy, np1_jdenx);
    u_lorenz_(np1_Ax, np1_Ay, np1_Az, np1_Bt, np1_dxdX, np1_dydY, np1_dzdZ, &Nx,
              &Ny, &Nz, &dx, &dy, &dz, phys_bdy, np1_lorenz);
    u_elecx_(np1_At, np1_Bx, np1_dxdX, &Nx, &Ny, &Nz, &dx, phys_bdy, np1_elecx);
    u_elecy_(np1_At, np1_By, np1_dydY, &Nx, &Ny, &Nz, &dy, phys_bdy, np1_elecy);
    u_elecz_(np1_At, np1_Bz, np1_dzdZ, &Nx, &Ny, &Nz, &dz, phys_bdy, np1_elecz);
    u_magx_(np1_Ay, np1_Az, np1_dydY, np1_dzdZ, &Nx, &Ny, &Nz, &dy, &dz,
            phys_bdy, np1_magx);
    u_magy_(np1_Ax, np1_Az, np1_dxdX, np1_dzdZ, &Nx, &Ny, &Nz, &dx, &dz,
            phys_bdy, np1_magy);
    u_magz_(np1_Ax, np1_Ay, np1_dxdX, np1_dydY, &Nx, &Ny, &Nz, &dx, &dy,
            phys_bdy, np1_magz);
    u_eem_(np1_elecx, np1_elecy, np1_elecz, np1_magx, np1_magy, np1_magz, &Nx,
           &Ny, &Nz, phys_bdy, np1_eem);
    u_gausse_(np1_At, np1_Bx, np1_By, np1_Bz, np1_dxdX, np1_dxdX2, np1_dydY,
              np1_dydY2, np1_dzdZ, np1_dzdZ2, np1_qden, &Nx, &Ny, &Nz, &e, &dx,
              &dy, &dz, phys_bdy, np1_gaussE);
    u_gaussb_(np1_dxdX, np1_dydY, np1_dzdZ, np1_magx, np1_magy, np1_magz, &Nx,
              &Ny, &Nz, &dx, &dy, &dz, phys_bdy, np1_gaussB);

    break;
  }
}

/***********************************************************************
 * qball_fill_ex_mask
 *
 * Sets excision mask.
 * Currently unused.
 *
 **********************************************************************/
void qball_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape,
                        int *shape_c, real *bbox, real excised) {
  return;
}

/***********************************************************************
 * qball_fill_bh_bboxes
 *
 * Fills excision arrays.
 * Currently unused.
 *
 **********************************************************************/
void qball_fill_bh_bboxes(real *bbox, int *num, int max_num) { return; }

/***********************************************************************
 * qball_post_tstep
 *
 * Performs calculations immediately after taking an evolution step.
 *
 **********************************************************************/
void qball_post_tstep(int L) {
  int i, valid, loopL, maxL, gridonlevelL = 0; // loop variables for grid iteration
  int base, half;            // integration region toggles
  int lshape[3];             // number of child grid points per dimension
  int lghost_width[6];       // width of child AMR ghost region per dimension
  real lbbox[6];             // stores child AMR bounding box data
  real ldx, ldy, ldz, ldt;   // child grid spacing in each dimension
  real ldx0[3];              // stores child grid spacings ldx, ldy, ldz
  real *lx0[3];              // stores pointers to coordinates x, y, z
  real *gfs0[PAMR_MAX_GFNS]; // stores pointers to child grid function data

  // AMR integrated quantities (charge, energy, angular momentum)
  real ltotalQ = 0.0, proctotalQ = 0.0, grandtotalQ = 0.0;
  real ltotalE = 0.0, proctotalE = 0.0, grandtotalE = 0.0;
  real ltotalJ = 0.0, proctotalJ = 0.0, grandtotalJ = 0.0;

  // AMR fractional integrated quantities (charge)
  real ltotalQhalf_upper = 0.0, proctotalQhalf_upper = 0.0;
  real ltotalQhalf_lower = 0.0, proctotalQhalf_lower = 0.0;
  real grandtotalQhalf_upper = 0.0, grandtotalQhalf_lower = 0.0;

  // unigrid integrated quantities (base level)
  real total_Qbase = 0.0, proctotal_Qbase = 0.0;
  real total_Ebase = 0.0, proctotal_Ebase = 0.0;
  real total_Jbase = 0.0, proctotal_Jbase = 0.0;
  real total_phi1IRE = 0.0, proctotal_phi1IRE = 0.0;
  real total_phi2IRE = 0.0, proctotal_phi2IRE = 0.0;
  real total_AtIRE = 0.0, proctotal_AtIRE = 0.0;
  real total_AxIRE = 0.0, proctotal_AxIRE = 0.0;
  real total_AyIRE = 0.0, proctotal_AyIRE = 0.0;
  real total_AzIRE = 0.0, proctotal_AzIRE = 0.0;
  real total_lorenz = 0.0, proctotal_lorenz = 0.0;
  real total_gaussE = 0.0, proctotal_gaussE = 0.0;
  real total_gaussB = 0.0, proctotal_gaussB = 0.0;

  // pointers for data output
  FILE *file_Qtot, *file_Etot, *file_Jtot;
  FILE *file_Qtothalf_upper, *file_Qtothalf_lower;
  FILE *file_Qtotbase, *file_Etotbase, *file_Jtotbase;
  FILE *file_phi1IRE, *file_phi2IRE;
  FILE *file_AtIRE, *file_AxIRE, *file_AyIRE, *file_AzIRE;
  FILE *file_lorenz, *file_gaussE, *file_gaussB;

  // pointers for receiving MPI data
  int num_procs;
  real *recvQ = NULL, *recvE = NULL, *recvJ = NULL;
  real *recvQhalf_upper = NULL, *recvQhalf_lower = NULL;
  real *recvQbase = NULL, *recvEbase = NULL, *recvJbase = NULL;
  real *recvIRE_phi1 = NULL, *recvIRE_phi2 = NULL;
  real *recvIRE_At = NULL, *recvIRE_Ax = NULL, *recvIRE_Ay = NULL, *recvIRE_Az = NULL;
  real *recvlorenz = NULL, *recvgaussE = NULL, *recvgaussB = NULL;

  // retrieve number of processors and allocate storage for receive buffer
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  if (my_rank == 0) {
    recvQ = (real *)calloc(num_procs, sizeof(real));
    recvE = (real *)calloc(num_procs, sizeof(real));
    recvJ = (real *)calloc(num_procs, sizeof(real));
    recvQhalf_upper = (real *)calloc(num_procs, sizeof(real));
    recvQhalf_lower = (real *)calloc(num_procs, sizeof(real));
    recvQbase = (real *)calloc(num_procs, sizeof(real));
    recvEbase = (real *)calloc(num_procs, sizeof(real));
    recvJbase = (real *)calloc(num_procs, sizeof(real));
    recvIRE_phi1 = (real *)calloc(num_procs, sizeof(real));
    recvIRE_phi2 = (real *)calloc(num_procs, sizeof(real));
    recvIRE_At = (real *)calloc(num_procs, sizeof(real));
    recvIRE_Ax = (real *)calloc(num_procs, sizeof(real));
    recvIRE_Ay = (real *)calloc(num_procs, sizeof(real));
    recvIRE_Az = (real *)calloc(num_procs, sizeof(real));
    recvlorenz = (real *)calloc(num_procs, sizeof(real));
    recvgaussE = (real *)calloc(num_procs, sizeof(real));
    recvgaussB = (real *)calloc(num_procs, sizeof(real));
  }

  // retrieve max AMR level in the hierarchy
  maxL = PAMR_get_max_lev(PAMR_AMRH);

  // if using AMR, compute base quantities on level 2 (always fully resolved)
  if (maxL == 1)
    loopL = 1;
  else
    loopL = 2;

  // to save resources, only do calculations on coarsest time steps
  if (L == 1) {

    // initialize iterator and loop through all grids on processor at level loopL
    valid = PAMR_init_s_iter(loopL, PAMR_AMRH, 0);
    while (valid) {

      // retrieve properties of the grid
      PAMR_get_g_ghost_width(lghost_width);
      PAMR_get_g_shape(lshape);
      PAMR_get_g_bbox(lbbox);
      PAMR_get_dxdt(loopL, ldx0, &ldt);
      ldx = ldx0[0];
      ldy = ldx0[1];
      ldz = ldx0[2];

      PAMR_get_g_x(lx0); // retrieve pointers to arrays containing coordinates

      // retrieve pointers to arrays containing grid function data
      PAMR_get_g_gfs(gfs0);

      nm1_modphi = gfs0[np1_modphi_gfn - 1];
      n_modphi = gfs0[nm1_modphi_gfn - 1];
      np1_modphi = gfs0[n_modphi_gfn - 1];

      nm1_phi1 = gfs0[np1_phi1_gfn - 1];
      n_phi1 = gfs0[nm1_phi1_gfn - 1];
      np1_phi1 = gfs0[n_phi1_gfn - 1];

      nm1_phi2 = gfs0[np1_phi2_gfn - 1];
      n_phi2 = gfs0[nm1_phi2_gfn - 1];
      np1_phi2 = gfs0[n_phi2_gfn - 1];

      nm1_pi1 = gfs0[np1_pi1_gfn - 1];
      n_pi1 = gfs0[nm1_pi1_gfn - 1];
      np1_pi1 = gfs0[n_pi1_gfn - 1];

      nm1_pi2 = gfs0[np1_pi2_gfn - 1];
      n_pi2 = gfs0[nm1_pi2_gfn - 1];
      np1_pi2 = gfs0[n_pi2_gfn - 1];

      nm1_At = gfs0[np1_At_gfn - 1];
      n_At = gfs0[nm1_At_gfn - 1];
      np1_At = gfs0[n_At_gfn - 1];

      nm1_Ax = gfs0[np1_Ax_gfn - 1];
      n_Ax = gfs0[nm1_Ax_gfn - 1];
      np1_Ax = gfs0[n_Ax_gfn - 1];

      nm1_Ay = gfs0[np1_Ay_gfn - 1];
      n_Ay = gfs0[nm1_Ay_gfn - 1];
      np1_Ay = gfs0[n_Ay_gfn - 1];

      nm1_Az = gfs0[np1_Az_gfn - 1];
      n_Az = gfs0[nm1_Az_gfn - 1];
      np1_Az = gfs0[n_Az_gfn - 1];

      nm1_Bt = gfs0[np1_Bt_gfn - 1];
      n_Bt = gfs0[nm1_Bt_gfn - 1];
      np1_Bt = gfs0[n_Bt_gfn - 1];

      nm1_Bx = gfs0[np1_Bx_gfn - 1];
      n_Bx = gfs0[nm1_Bx_gfn - 1];
      np1_Bx = gfs0[n_Bx_gfn - 1];

      nm1_By = gfs0[np1_By_gfn - 1];
      n_By = gfs0[nm1_By_gfn - 1];
      np1_By = gfs0[n_By_gfn - 1];

      nm1_Bz = gfs0[np1_Bz_gfn - 1];
      n_Bz = gfs0[nm1_Bz_gfn - 1];
      np1_Bz = gfs0[n_Bz_gfn - 1];

      nm1_chi = gfs0[np1_chi_gfn - 1];
      n_chi = gfs0[nm1_chi_gfn - 1];
      np1_chi = gfs0[n_chi_gfn - 1];

      nm1_xi = gfs0[np1_xi_gfn - 1];
      n_xi = gfs0[nm1_xi_gfn - 1];
      np1_xi = gfs0[n_xi_gfn - 1];

      // note change in cyclic switching for derived quantities
      nm1_qden = gfs0[nm1_qden_gfn - 1];
      n_qden = gfs0[n_qden_gfn - 1];
      np1_qden = gfs0[np1_qden_gfn - 1];

      nm1_eden = gfs0[nm1_eden_gfn - 1];
      n_eden = gfs0[n_eden_gfn - 1];
      np1_eden = gfs0[np1_eden_gfn - 1];

      nm1_jdenx = gfs0[nm1_jdenx_gfn - 1];
      n_jdenx = gfs0[n_jdenx_gfn - 1];
      np1_jdenx = gfs0[np1_jdenx_gfn - 1];

      nm1_lorenz = gfs0[nm1_lorenz_gfn - 1];
      n_lorenz = gfs0[n_lorenz_gfn - 1];
      np1_lorenz = gfs0[np1_lorenz_gfn - 1];

      nm1_elecx = gfs0[nm1_elecx_gfn - 1];
      n_elecx = gfs0[n_elecx_gfn - 1];
      np1_elecx = gfs0[np1_elecx_gfn - 1];

      nm1_elecy = gfs0[nm1_elecy_gfn - 1];
      n_elecy = gfs0[n_elecy_gfn - 1];
      np1_elecy = gfs0[np1_elecy_gfn - 1];

      nm1_elecz = gfs0[nm1_elecz_gfn - 1];
      n_elecz = gfs0[n_elecz_gfn - 1];
      np1_elecz = gfs0[np1_elecz_gfn - 1];

      nm1_magx = gfs0[nm1_magx_gfn - 1];
      n_magx = gfs0[n_magx_gfn - 1];
      np1_magx = gfs0[np1_magx_gfn - 1];

      nm1_magy = gfs0[nm1_magy_gfn - 1];
      n_magy = gfs0[n_magy_gfn - 1];
      np1_magy = gfs0[np1_magy_gfn - 1];

      nm1_magz = gfs0[nm1_magz_gfn - 1];
      n_magz = gfs0[n_magz_gfn - 1];
      np1_magz = gfs0[np1_magz_gfn - 1];

      nm1_eem = gfs0[nm1_eem_gfn - 1];
      n_eem = gfs0[n_eem_gfn - 1];
      np1_eem = gfs0[np1_eem_gfn - 1];

      nm1_gaussE = gfs0[nm1_gaussE_gfn - 1];
      n_gaussE = gfs0[n_gaussE_gfn - 1];
      np1_gaussE = gfs0[np1_gaussE_gfn - 1];

      nm1_gaussB = gfs0[nm1_gaussB_gfn - 1];
      n_gaussB = gfs0[n_gaussB_gfn - 1];
      np1_gaussB = gfs0[np1_gaussB_gfn - 1];

      nm1_Xx = gfs0[nm1_Xx_gfn - 1];
      n_Xx = gfs0[n_Xx_gfn - 1];
      np1_Xx = gfs0[np1_Xx_gfn - 1];

      nm1_dxdX = gfs0[nm1_dxdX_gfn - 1];
      n_dxdX = gfs0[n_dxdX_gfn - 1];
      np1_dxdX = gfs0[np1_dxdX_gfn - 1];

      nm1_dxdX2 = gfs0[nm1_dxdX2_gfn - 1];
      n_dxdX2 = gfs0[n_dxdX2_gfn - 1];
      np1_dxdX2 = gfs0[np1_dxdX2_gfn - 1];

      nm1_Yy = gfs0[nm1_Yy_gfn - 1];
      n_Yy = gfs0[n_Yy_gfn - 1];
      np1_Yy = gfs0[np1_Yy_gfn - 1];

      nm1_dydY = gfs0[nm1_dydY_gfn - 1];
      n_dydY = gfs0[n_dydY_gfn - 1];
      np1_dydY = gfs0[np1_dydY_gfn - 1];

      nm1_dydY2 = gfs0[nm1_dydY2_gfn - 1];
      n_dydY2 = gfs0[n_dydY2_gfn - 1];
      np1_dydY2 = gfs0[np1_dydY2_gfn - 1];

      nm1_Zz = gfs0[nm1_Zz_gfn - 1];
      n_Zz = gfs0[n_Zz_gfn - 1];
      np1_Zz = gfs0[np1_Zz_gfn - 1];

      nm1_dzdZ = gfs0[nm1_dzdZ_gfn - 1];
      n_dzdZ = gfs0[n_dzdZ_gfn - 1];
      np1_dzdZ = gfs0[np1_dzdZ_gfn - 1];

      nm1_dzdZ2 = gfs0[nm1_dzdZ2_gfn - 1];
      n_dzdZ2 = gfs0[n_dzdZ2_gfn - 1];
      np1_dzdZ2 = gfs0[np1_dzdZ2_gfn - 1];

      nm1_jac = gfs0[nm1_jac_gfn - 1];
      n_jac = gfs0[n_jac_gfn - 1];
      np1_jac = gfs0[np1_jac_gfn - 1];

      cmask = gfs0[cmask_gfn - 1];

      // set region toggles for fractional integration
      base = 1;
      half = 0;

      // compute desired quantities for each processor on the base level
      proctotal_Qbase =
          qtotcalc_(cmask, np1_qden, n_qden, n_jac, &lshape[0], &lshape[1],
                    &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                    lghost_width, &base, &half, &c, &d);

      proctotal_Ebase =
          etotcalc_(cmask, np1_eden, n_eden, n_jac, &lshape[0], &lshape[1],
                    &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                    lghost_width, &base, &c, &d);

      proctotal_Jbase =
          etotcalc_(cmask, np1_jdenx, n_jdenx, n_jac, &lshape[0], &lshape[1],
                    &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                    lghost_width, &base, &c, &d);

      // compute partial L2-norm of desired quantities for each processor on the base level
      // (this ignores ghost points, so will overcount)
      if (V_type == 0) {
        ire_phi1_poly_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                       n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi2,
                       nm1_phi1, np1_phi1, &lshape[0], &lshape[1], &lshape[2],
                       &cc, &e, &g, &h, &ldt, &ldx, &ldy, &ldz, &m,
                       &proctotal_phi1IRE);
        ire_phi2_poly_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                       n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi1,
                       nm1_phi2, np1_phi2, &lshape[0], &lshape[1], &lshape[2],
                       &cc, &e, &g, &h, &ldt, &ldx, &ldy, &ldz, &m,
                       &proctotal_phi2IRE);
      } else {
        ire_phi1_log_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                      n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi2, nm1_phi1,
                      np1_phi1, &lshape[0], &lshape[1], &lshape[2], &beta, &cc,
                      &e, &ldt, &ldx, &ldy, &ldz, &mu, &proctotal_phi1IRE);
        ire_phi2_log_(n_At, n_Ax, n_Ay, n_Az, n_chi, n_dxdX, n_dxdX2, n_dydY,
                      n_dydY2, n_dzdZ, n_dzdZ2, n_phi1, n_phi2, n_pi1, nm1_phi2,
                      np1_phi2, &lshape[0], &lshape[1], &lshape[2], &beta, &cc,
                      &e, &ldt, &ldx, &ldy, &ldz, &mu, &proctotal_phi2IRE);
      }
      ire_at_(n_At, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
              n_phi2, n_pi1, n_pi2, nm1_At, np1_At, &lshape[0], &lshape[1],
              &lshape[2], &e, &ldt, &ldx, &ldy, &ldz, &proctotal_AtIRE);
      ire_ax_(n_Ax, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
              n_phi2, nm1_Ax, np1_Ax, &lshape[0], &lshape[1], &lshape[2], &e,
              &ldt, &ldx, &ldy, &ldz, &proctotal_AxIRE);
      ire_ay_(n_Ay, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
              n_phi2, nm1_Ay, np1_Ay, &lshape[0], &lshape[1], &lshape[2], &e,
              &ldt, &ldx, &ldy, &ldz, &proctotal_AyIRE);
      ire_az_(n_Az, n_dxdX, n_dxdX2, n_dydY, n_dydY2, n_dzdZ, n_dzdZ2, n_phi1,
              n_phi2, nm1_Az, np1_Az, &lshape[0], &lshape[1], &lshape[2], &e,
              &ldt, &ldx, &ldy, &ldz, &proctotal_AzIRE);

      proctotal_lorenz = l2norm_calc(np1_lorenz);
      proctotal_gaussE = l2norm_calc(np1_gaussE);
      proctotal_gaussB = l2norm_calc(np1_gaussB);

      valid = PAMR_next_g(); // retrieve next grid on the specified level for this processor
    }

    // wait for all processors to reach this point
    MPI_Barrier(MPI_COMM_WORLD);

    // fill recv arrays with values from each processor
    MPI_Gather(&proctotal_phi1IRE, 1, MPI_DOUBLE, recvIRE_phi1, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_phi2IRE, 1, MPI_DOUBLE, recvIRE_phi2, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_AtIRE, 1, MPI_DOUBLE, recvIRE_At, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_AxIRE, 1, MPI_DOUBLE, recvIRE_Ax, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_AyIRE, 1, MPI_DOUBLE, recvIRE_Ay, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_AzIRE, 1, MPI_DOUBLE, recvIRE_Az, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_lorenz, 1, MPI_DOUBLE, recvlorenz, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_gaussE, 1, MPI_DOUBLE, recvgaussE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_gaussB, 1, MPI_DOUBLE, recvgaussB, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Qbase, 1, MPI_DOUBLE, recvQbase, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Ebase, 1, MPI_DOUBLE, recvEbase, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Jbase, 1, MPI_DOUBLE, recvJbase, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    // add up totals
    if (my_rank == 0) {
      total_phi1IRE = pow(proctotal_phi1IRE, 2);
      total_phi2IRE = pow(proctotal_phi2IRE, 2);
      total_AtIRE = pow(proctotal_AtIRE, 2);
      total_AxIRE = pow(proctotal_AxIRE, 2);
      total_AyIRE = pow(proctotal_AyIRE, 2);
      total_AzIRE = pow(proctotal_AzIRE, 2);
      total_lorenz = pow(proctotal_lorenz, 2);
      total_gaussE = pow(proctotal_gaussE, 2);
      total_gaussB = pow(proctotal_gaussB, 2);
      total_Qbase = proctotal_Qbase;
      total_Ebase = proctotal_Ebase;
      total_Jbase = proctotal_Jbase;

      for (i = 1; i < num_procs; i++) {
        total_phi1IRE = total_phi1IRE + pow(recvIRE_phi1[num_procs - i], 2);
        total_phi2IRE = total_phi2IRE + pow(recvIRE_phi2[num_procs - i], 2);
        total_AtIRE = total_AtIRE + pow(recvIRE_At[num_procs - i], 2);
        total_AxIRE = total_AxIRE + pow(recvIRE_Ax[num_procs - i], 2);
        total_AyIRE = total_AyIRE + pow(recvIRE_Ay[num_procs - i], 2);
        total_AzIRE = total_AzIRE + pow(recvIRE_Az[num_procs - i], 2);
        total_lorenz = total_lorenz + pow(recvlorenz[num_procs - i], 2);
        total_gaussE = total_gaussE + pow(recvgaussE[num_procs - i], 2);
        total_gaussB = total_gaussB + pow(recvgaussB[num_procs - i], 2);
        total_Qbase = total_Qbase + recvQbase[num_procs - i];
        total_Ebase = total_Ebase + recvEbase[num_procs - i];
        total_Jbase = total_Jbase + recvJbase[num_procs - i];
      }

      // stop evolution if NaN detected
      if (total_Qbase != total_Qbase) {
        AMRD_stop("NaN detected in app_post_tstep()... AMRD_stop", 0);
      }

      // write data to file
      file_phi1IRE = fopen("./IRE_phi1.dat", "a");
      fprintf(file_phi1IRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_phi1IRE));
      fclose(file_phi1IRE);

      file_phi2IRE = fopen("./IRE_phi2.dat", "a");
      fprintf(file_phi2IRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_phi2IRE));
      fclose(file_phi2IRE);

      file_AtIRE = fopen("./IRE_At.dat", "a");
      fprintf(file_AtIRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_AtIRE));
      fclose(file_AtIRE);

      file_AxIRE = fopen("./IRE_Ax.dat", "a");
      fprintf(file_AxIRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_AxIRE));
      fclose(file_AxIRE);

      file_AyIRE = fopen("./IRE_Ay.dat", "a");
      fprintf(file_AyIRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_AyIRE));
      fclose(file_AyIRE);

      file_AzIRE = fopen("./IRE_Az.dat", "a");
      fprintf(file_AzIRE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_AzIRE));
      fclose(file_AzIRE);

      file_lorenz = fopen("./lorenz.dat", "a");
      fprintf(file_lorenz, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_lorenz));
      fclose(file_lorenz);

      file_gaussE = fopen("./gaussE.dat", "a");
      fprintf(file_gaussE, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_gaussE));
      fclose(file_gaussE);

      file_gaussB = fopen("./gaussB.dat", "a");
      fprintf(file_gaussB, "%20.15f  %20.15f\n", PAMR_get_time(L), sqrt(total_gaussB));
      fclose(file_gaussB);

      file_Qtotbase = fopen("./Qtotbase.dat", "a");
      fprintf(file_Qtotbase, "%20.15f  %20.15f\n", PAMR_get_time(L), total_Qbase);
      fclose(file_Qtotbase);

      file_Etotbase = fopen("./Etotbase.dat", "a");
      fprintf(file_Etotbase, "%20.15f  %20.15f\n", PAMR_get_time(L), total_Ebase);
      fclose(file_Etotbase);

      file_Jtotbase = fopen("./Jtotbase.dat", "a");
      fprintf(file_Jtotbase, "%20.15f  %20.15f\n", PAMR_get_time(L), total_Jbase);
      fclose(file_Jtotbase);
    }

    // loop through all levels in the AMR hierarchy on each processor
    for (loopL = 1; loopL <= maxL; loopL = loopL + 1) {

      gridonlevelL = 1;

      // initialize iterator and loop through all grids on processor at level loopL
      PAMR_push_iter();
      valid = PAMR_init_s_iter(loopL, PAMR_AMRH, 0);

      // wherever child grid exists, set child mask
      while (valid) {
        PAMR_push_iter();
        set_cmask_child(loopL, PAMR_AMRH);
        PAMR_pop_iter();

        // retrieve properties of the child grid
        PAMR_get_g_ghost_width(lghost_width);
        PAMR_get_g_shape(lshape);
        PAMR_get_g_bbox(lbbox);
        PAMR_get_dxdt(loopL, ldx0, &ldt);
        ldx = ldx0[0];
        ldy = ldx0[1];
        ldz = ldx0[2];

        PAMR_get_g_x(lx0); // retrieve pointers to arrays containing coordinates

        // retrieve pointers to arrays containing grid function data
        PAMR_get_g_gfs(gfs0);

        nm1_modphi = gfs0[np1_modphi_gfn - 1];
        n_modphi = gfs0[nm1_modphi_gfn - 1];
        np1_modphi = gfs0[n_modphi_gfn - 1];

        nm1_phi1 = gfs0[np1_phi1_gfn - 1];
        n_phi1 = gfs0[nm1_phi1_gfn - 1];
        np1_phi1 = gfs0[n_phi1_gfn - 1];

        nm1_phi2 = gfs0[np1_phi2_gfn - 1];
        n_phi2 = gfs0[nm1_phi2_gfn - 1];
        np1_phi2 = gfs0[n_phi2_gfn - 1];

        nm1_pi1 = gfs0[np1_pi1_gfn - 1];
        n_pi1 = gfs0[nm1_pi1_gfn - 1];
        np1_pi1 = gfs0[n_pi1_gfn - 1];

        nm1_pi2 = gfs0[np1_pi2_gfn - 1];
        n_pi2 = gfs0[nm1_pi2_gfn - 1];
        np1_pi2 = gfs0[n_pi2_gfn - 1];

        nm1_At = gfs0[nm1_At_gfn - 1];
        n_At = gfs0[n_At_gfn - 1];
        np1_At = gfs0[np1_At_gfn - 1];

        nm1_Ax = gfs0[nm1_Ax_gfn - 1];
        n_Ax = gfs0[n_Ax_gfn - 1];
        np1_Ax = gfs0[np1_Ax_gfn - 1];

        nm1_Ay = gfs0[nm1_Ay_gfn - 1];
        n_Ay = gfs0[n_Ay_gfn - 1];
        np1_Ay = gfs0[np1_Ay_gfn - 1];

        nm1_Az = gfs0[nm1_Az_gfn - 1];
        n_Az = gfs0[n_Az_gfn - 1];
        np1_Az = gfs0[np1_Az_gfn - 1];

        nm1_Bt = gfs0[nm1_Bt_gfn - 1];
        n_Bt = gfs0[n_Bt_gfn - 1];
        np1_Bt = gfs0[np1_Bt_gfn - 1];

        nm1_Bx = gfs0[nm1_Bx_gfn - 1];
        n_Bx = gfs0[n_Bx_gfn - 1];
        np1_Bx = gfs0[np1_Bx_gfn - 1];

        nm1_By = gfs0[nm1_By_gfn - 1];
        n_By = gfs0[n_By_gfn - 1];
        np1_By = gfs0[np1_By_gfn - 1];

        nm1_Bz = gfs0[nm1_Bz_gfn - 1];
        n_Bz = gfs0[n_Bz_gfn - 1];
        np1_Bz = gfs0[np1_Bz_gfn - 1];

        nm1_chi = gfs0[nm1_chi_gfn - 1];
        n_chi = gfs0[n_chi_gfn - 1];
        np1_chi = gfs0[np1_chi_gfn - 1];

        nm1_xi = gfs0[nm1_xi_gfn - 1];
        n_xi = gfs0[n_xi_gfn - 1];
        np1_xi = gfs0[np1_xi_gfn - 1];

        // note change in cyclic switching for derived quantities
        nm1_jdenx = gfs0[nm1_jdenx_gfn - 1];
        n_jdenx = gfs0[n_jdenx_gfn - 1];
        np1_jdenx = gfs0[np1_jdenx_gfn - 1];

        nm1_eden = gfs0[nm1_eden_gfn - 1];
        n_eden = gfs0[n_eden_gfn - 1];
        np1_eden = gfs0[np1_eden_gfn - 1];

        nm1_qden = gfs0[nm1_qden_gfn - 1];
        n_qden = gfs0[n_qden_gfn - 1];
        np1_qden = gfs0[np1_qden_gfn - 1];

        nm1_Xx = gfs0[nm1_Xx_gfn - 1];
        n_Xx = gfs0[n_Xx_gfn - 1];
        np1_Xx = gfs0[np1_Xx_gfn - 1];

        nm1_dxdX = gfs0[nm1_dxdX_gfn - 1];
        n_dxdX = gfs0[n_dxdX_gfn - 1];
        np1_dxdX = gfs0[np1_dxdX_gfn - 1];

        nm1_dxdX2 = gfs0[nm1_dxdX2_gfn - 1];
        n_dxdX2 = gfs0[n_dxdX2_gfn - 1];
        np1_dxdX2 = gfs0[np1_dxdX2_gfn - 1];

        nm1_Yy = gfs0[nm1_Yy_gfn - 1];
        n_Yy = gfs0[n_Yy_gfn - 1];
        np1_Yy = gfs0[np1_Yy_gfn - 1];

        nm1_dydY = gfs0[nm1_dydY_gfn - 1];
        n_dydY = gfs0[n_dydY_gfn - 1];
        np1_dydY = gfs0[np1_dydY_gfn - 1];

        nm1_dydY2 = gfs0[nm1_dydY2_gfn - 1];
        n_dydY2 = gfs0[n_dydY2_gfn - 1];
        np1_dydY2 = gfs0[np1_dydY2_gfn - 1];

        nm1_Zz = gfs0[nm1_Zz_gfn - 1];
        n_Zz = gfs0[n_Zz_gfn - 1];
        np1_Zz = gfs0[np1_Zz_gfn - 1];

        nm1_dzdZ = gfs0[nm1_dzdZ_gfn - 1];
        n_dzdZ = gfs0[n_dzdZ_gfn - 1];
        np1_dzdZ = gfs0[np1_dzdZ_gfn - 1];

        nm1_dzdZ2 = gfs0[nm1_dzdZ2_gfn - 1];
        n_dzdZ2 = gfs0[n_dzdZ2_gfn - 1];
        np1_dzdZ2 = gfs0[np1_dzdZ2_gfn - 1];

        nm1_jac = gfs0[nm1_jac_gfn - 1];
        n_jac = gfs0[n_jac_gfn - 1];
        np1_jac = gfs0[np1_jac_gfn - 1];

        cmask = gfs0[cmask_gfn - 1];

        // set region toggles for fractional integration
        base = 0;
        half = 0;

        // compute desired quantities for each processor on the base level
        proctotalQ =
            qtotcalc_(cmask, np1_qden, n_qden, n_jac, &lshape[0], &lshape[1],
                      &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                      lghost_width, &base, &half, &c, &d);

        proctotalE = etotcalc_(cmask, np1_eden, n_eden, n_jac, &lshape[0],
                               &lshape[1], &lshape[2], lx0[0], lx0[1], lx0[2],
                               &ldx, &ldy, &ldz, lghost_width, &base, &c, &d);

        proctotalJ = etotcalc_(cmask, np1_jdenx, n_jdenx, n_jac, &lshape[0],
                               &lshape[1], &lshape[2], lx0[0], lx0[1], lx0[2],
                               &ldx, &ldy, &ldz, lghost_width, &base, &c, &d);

        // integrate in upper half-volume
        half = 1;
        ltotalQhalf_upper = qtotcalc_(cmask, np1_qden, n_qden, n_jac, &lshape[0], &lshape[1],
                                      &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                                      lghost_width, &base, &half, &c, &d);

        // integrate in lower half-volume
        half = -1;
        ltotalQhalf_lower = qtotcalc_(cmask, np1_qden, n_qden, n_jac, &lshape[0], &lshape[1],
                                      &lshape[2], lx0[0], lx0[1], lx0[2], &ldx, &ldy, &ldz,
                                      lghost_width, &base, &half, &c, &d);

        // add contributions from this level to processor total
        proctotalQ += ltotalQ;
        proctotalE += ltotalE;
        proctotalJ += ltotalJ;
        proctotalQhalf_upper += ltotalQhalf_upper;
        proctotalQhalf_lower += ltotalQhalf_lower;

        gridonlevelL += 1;
        valid = PAMR_next_g(); // retrieve next grid on the specified level for this processor
      }

      PAMR_pop_iter();
    }

    // wait for all processors to reach this point
    MPI_Barrier(MPI_COMM_WORLD);

    // fill recv arrays with values from each processor
    MPI_Gather(&proctotalQ, 1, MPI_DOUBLE, recvQ, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotalE, 1, MPI_DOUBLE, recvE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotalJ, 1, MPI_DOUBLE, recvJ, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotalQhalf_upper, 1, MPI_DOUBLE, recvQhalf_upper, 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&proctotalQhalf_lower, 1, MPI_DOUBLE, recvQhalf_lower, 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // add up totals
    if (my_rank == 0) {
      grandtotalQ += proctotalQ;
      grandtotalE += proctotalE;
      grandtotalJ += proctotalJ;
      grandtotalQhalf_upper += proctotalQhalf_upper;
      grandtotalQhalf_lower += proctotalQhalf_lower;

      for (i = 1; i < num_procs; i++) {
        grandtotalQ = grandtotalQ + recvQ[num_procs - i];
        grandtotalE = grandtotalE + recvE[num_procs - i];
        grandtotalJ = grandtotalJ + recvJ[num_procs - i];
        grandtotalQhalf_upper =
            grandtotalQhalf_upper + recvQhalf_upper[num_procs - i];
        grandtotalQhalf_lower =
            grandtotalQhalf_lower + recvQhalf_lower[num_procs - i];
      }

      // write data to file
      file_Qtot = fopen("./Qtot.dat", "a");
      fprintf(file_Qtot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalQ);
      fclose(file_Qtot);

      file_Etot = fopen("./Etot.dat", "a");
      fprintf(file_Etot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalE);
      fclose(file_Etot);

      file_Jtot = fopen("./Jtot.dat", "a");
      fprintf(file_Jtot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalJ);
      fclose(file_Jtot);

      file_Qtothalf_upper = fopen("./Qtothalf_upper.dat", "a");
      fprintf(file_Qtothalf_upper, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalQhalf_upper);
      fclose(file_Qtothalf_upper);

      file_Qtothalf_lower = fopen("./Qtothalf_lower.dat", "a");
      fprintf(file_Qtothalf_lower, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalQhalf_lower);
      fclose(file_Qtothalf_lower);
    }
  }

  // print elapsed time
  if (my_rank == 0 && L == 1)
    elapsed_time();

  if (my_rank == 0) {
    free(recvQ);
    free(recvE);
    free(recvJ);
    free(recvQhalf_upper);
    free(recvQhalf_lower);
    free(recvQbase);
    free(recvEbase);
    free(recvJbase);
    free(recvIRE_phi1);
    free(recvIRE_phi2);
  }
}

/***********************************************************************
 * qball_MG_residual
 *
 * Returns a norm of the residual and stores point-wise residuals for
 * multigrid variables.
 *
 ***********************************************************************/
real qball_MG_residual(void) {

  // if multigrid data file exists or mg_toggle turned off, skip calculations
  if ((access(mg_data, F_OK) == 0) || (mg_toggle == 0))
    return 0;

  real norm = 0.0;

  ldptr();

  // compute residuals
  residual_(V_res, V_rhs, V, defect, mask_mg, x, y, z, &norm, &Nx, &Ny, &Nz,
            &dx, &dy, &dz, Xx, Yy, Zz, dxdX, dxdX2, dydY, dydY2, dzdZ, dzdZ2,
            &e, phi1, phi2, pi1, pi2, modphi, At);

  // compute multigrid defect correction
  if (MGiter == 0 && Nx == MGshape) {
    MGiter = 1;
  } else if (MGiter == 1 && Nx == MGshape) {
    MGiter = 0;
    dc_(V, defect, x, y, z, &dx, &dy, &dz, &Nx, &Ny, &Nz, Xx, Yy, Zz, dxdX,
        dxdX2, dydY, dydY2, dzdZ, dzdZ2, &e, phi1, phi2, pi1, pi2, modphi, At);
    printf("\nApplied defect correction!\n");
  }

  return norm;
}

/***********************************************************************
 * qball_MG_relax
 *
 * Performs one relaxation sweep of the multigrid equations and returns
 * an estimate of the residual norm.
 *
 ***********************************************************************/
real qball_MG_relax(void) {

  // if multigrid data file exists or mg_toggle turned off, skip calculations
  if ((access(mg_data, F_OK) == 0) || (mg_toggle == 0))
    return 0;

  real norm = 0.0;

  ldptr();

  // perform relaxation sweep
  relax_(V, V_rhs, defect, mask_mg, phys_bdy, &norm, x, y, z, &Nx, &Ny, &Nz,
         &dx, &dy, &dz, Xx, Yy, Zz, dxdX, dxdX2, dydY, dydY2, dzdZ, dzdZ2, &e,
         phi1, phi2, pi1, pi2, modphi, At);

  return norm;
}

/***********************************************************************
 * qball_L_op
 *
 * Computes the elliptic differential operator acting on each multigrid
 * variable and stores the result.
 *
 ***********************************************************************/
void qball_L_op(void) {

  // if multigrid data file exists or mg_toggle turned off, skip calculations
  if ((access(mg_data, F_OK) == 0) || (mg_toggle == 0))
    return;

  ldptr();

  // compute differential operator
  lop_(V_lop, V, defect, mask_mg, x, y, z, &Nx, &Ny, &Nz, &dx, &dy, &dz, Xx, Yy,
       Zz, dxdX, dxdX2, dydY, dydY2, dzdZ, dzdZ2, &e, phi1, phi2, pi1, pi2,
       modphi, At);
}

/***********************************************************************
 * qball_scale_tre
 *
 * Modifies the truncation error estimates.
 * Currently unused.
 *
 ***********************************************************************/
void qball_scale_tre(void) { return; }

/***********************************************************************
 * qball_post_regrid
 *
 * Performs calculations immediately after an AMR regrid.
 * Currently unused.
 *
 ***********************************************************************/
void qball_post_regrid(void) { return; }

/***********************************************************************
 * main
 *
 * Calls the AMRD driver routine and reports the total execution time.
 ***********************************************************************/
int main(int argc, char **argv) {
//  // enables SIGFPE on floating-point exceptions (needs fenv.h)
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  amrd(argc, argv, &qball_id, &qball_var_pre_init, &qball_var_post_init,
       &qball_AMRH_var_clear, &qball_free_data, &qball_t0_cnst_data,
       &qball_evo_residual, &qball_MG_residual, &qball_evolve, &qball_MG_relax,
       &qball_L_op, &qball_pre_io_calc, &qball_scale_tre, &qball_post_regrid,
       &qball_post_tstep, &qball_fill_ex_mask, &qball_fill_bh_bboxes);
  if (my_rank == 0)
    elapsed_time();

  return 0;
}
