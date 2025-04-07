//======================================================================
// qball.c
//
// This program interfaces with the PAMR and AMRD libraries to solve the
// gauged Q-ball equations of motion using adaptive mesh refinement and
// multigrid in axisymmetry.
//
// For more information about built-in PAMR/AMRD functionality, look at
// the PAMR or AMRD reference manual, pamr.h, amrd.h, or pamr.c
//======================================================================

// #include <fenv.h>  // uncomment when debugging FPEs
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

// model parameters (values set in input parameter files or init_qball.inc)
real c, d;
real cR, cX, cR1, cR2, cX1, cX2;
real epsdis;
real e, w;
real beta, mu;
real h, mm, gg;
real v;
real amp, delta, R_center, Rw, r0, signum, X_center, Xw;
real cc;
int Vtype;

// pointers for physical coordinate values
real *R, *X;

// pointers for grid function values
real *modphi_n, *modphi_np1;
real *phi1_n, *phi1_np1;
real *phi2_n, *phi2_np1;
real *pi1_n, *pi1_np1;
real *pi2_n, *pi2_np1;
real *A_t_n, *A_t_np1;
real *Atilde_p_n, *Atilde_p_np1;
real *Atilde_phi_n, *Atilde_phi_np1;
real *A_z_n, *A_z_np1;
real *B_t_n, *B_t_np1;
real *Btilde_p_n, *Btilde_p_np1;
real *Btilde_phi_n, *Btilde_phi_np1;
real *B_z_n, *B_z_np1;
real *gaugecond_n, *gaugecond_np1;
real *divE_n, *divE_np1;
real *divB_n, *divB_np1;
real *chi_n, *chi_np1;
real *xi_n, *xi_np1;
real *Eden_n, *Eden_np1;
real *Eem_n, *Eem_np1;
real *Qden_n, *Qden_np1;
real *elec_p_n, *elec_p_np1;
real *elec_phi_n, *elec_phi_np1;
real *elec_z_n, *elec_z_np1;
real *mag_p_n, *mag_p_np1;
real *mag_phi_n, *mag_phi_np1;
real *mag_z_n, *mag_z_np1;

real *dRdp_n;
real *dRdp2_n;
real *dXdz_n;
real *dXdz2_n;
real *pR_n;
real *zX_n;

real *chi_res;
real *xi_res;
real *phi1_res;
real *phi2_res;
real *pi1_res;
real *pi2_res;
real *A_t_res;
real *Atilde_p_res;
real *Atilde_phi_res;
real *A_z_res;
real *B_t_res;
real *Btilde_p_res;
real *Btilde_phi_res;
real *B_z_res;

real *phi1ire_n, *phi1ire_np1;
real *phi2ire_n, *phi2ire_np1;
real *atire_n, *atire_np1;
real *apire_n, *apire_np1;
real *azire_n, *azire_np1;

real *cmask;

real *dummy_n;

// pointers for multigrid functionality
real *V, *V_lop, *V_res, *V_rhs;
real *V_n, *V_np1;
real *mask, *mask_mg;
real *dRdp;
real *dRdp2;
real *dXdz;
real *dXdz2;
real *pR;
real *zX;
real *phi1, *phi2, *pi1, *pi2, *modphi, *A_t;

// grid function numbers (gfn)
int modphi_n_gfn, modphi_np1_gfn;
int phi1_n_gfn, phi1_np1_gfn;
int phi2_n_gfn, phi2_np1_gfn;
int pi1_n_gfn, pi1_np1_gfn;
int pi2_n_gfn, pi2_np1_gfn;
int A_t_n_gfn, A_t_np1_gfn;
int Atilde_p_n_gfn, Atilde_p_np1_gfn;
int Atilde_phi_n_gfn, Atilde_phi_np1_gfn;
int A_z_n_gfn, A_z_np1_gfn;
int B_t_n_gfn, B_t_np1_gfn;
int Btilde_p_n_gfn, Btilde_p_np1_gfn;
int Btilde_phi_n_gfn, Btilde_phi_np1_gfn;
int B_z_n_gfn, B_z_np1_gfn;
int gaugecond_n_gfn, gaugecond_np1_gfn;
int divE_n_gfn, divE_np1_gfn;
int divB_n_gfn, divB_np1_gfn;
int chi_n_gfn, chi_np1_gfn;
int xi_n_gfn, xi_np1_gfn;
int Eden_n_gfn, Eden_np1_gfn;
int Eem_n_gfn, Eem_np1_gfn;
int Qden_n_gfn, Qden_np1_gfn;
int elec_p_n_gfn, elec_p_np1_gfn;
int elec_phi_n_gfn, elec_phi_np1_gfn;
int elec_z_n_gfn, elec_z_np1_gfn;
int mag_p_n_gfn, mag_p_np1_gfn;
int mag_phi_n_gfn, mag_phi_np1_gfn;
int mag_z_n_gfn, mag_z_np1_gfn;

int dRdp_n_gfn;
int dRdp2_n_gfn;
int dXdz_n_gfn;
int dXdz2_n_gfn;
int pR_n_gfn;
int zX_n_gfn;

int chi_res_gfn, xi_res_gfn;
int phi1_res_gfn, phi2_res_gfn, pi1_res_gfn, pi2_res_gfn;
int A_t_res_gfn, Atilde_p_res_gfn, Atilde_phi_res_gfn, A_z_res_gfn;
int B_t_res_gfn, Btilde_p_res_gfn, Btilde_phi_res_gfn, B_z_res_gfn;

int phi1ire_n_gfn, phi1ire_np1_gfn;
int phi2ire_n_gfn, phi2ire_np1_gfn;
int atire_n_gfn, atire_np1_gfn;
int apire_n_gfn, apire_np1_gfn;
int azire_n_gfn, azire_np1_gfn;

int cmask_gfn;

int dummy_n_gfn;

// grid function numbers (gfn) for multigrid
int V_gfn, V_lop_gfn, V_res_gfn, V_rhs_gfn, mask_gfn, mask_mg_gfn;
int V_n_gfn, V_np1_gfn;
int dRdp_gfn;
int dRdp2_gfn;
int dXdz_gfn;
int dXdz2_gfn;
int pR_gfn;
int zX_gfn;
int phi1_gfn, phi2_gfn, pi1_gfn, pi2_gfn, modphi_gfn, A_t_gfn;

// toggles for multigrid functionality
int mg_toggle;
char mg_data[] = "V.bin";

int dim;             // grid dimensionality
int size, *gridsize; // grid size
int shape[2];        // number of grid points per dimension
int ghost_width[4];  // width of AMR ghost region per dimension
int phys_bdy[4];     // toggle for AMR/physical boundaries
int g_rank;          // MPI rank of execution
int g_L;             // current level in AMR hierarchy
int NR, NX;          // number of grid points in each dimension
real dR, dX, dt;     // grid spacing in each dimension
real base_bbox[4];   // stores physical bounding box data
real bbox[4];        // stores AMR bounding box data
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
  if ((modphi_n_gfn = PAMR_get_gfn("modphi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((modphi_np1_gfn = PAMR_get_gfn("modphi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((phi1_n_gfn = PAMR_get_gfn("phi1", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1_np1_gfn = PAMR_get_gfn("phi1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((phi2_n_gfn = PAMR_get_gfn("phi2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_np1_gfn = PAMR_get_gfn("phi2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((pi1_n_gfn = PAMR_get_gfn("pi1", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_np1_gfn = PAMR_get_gfn("pi1", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((pi2_n_gfn = PAMR_get_gfn("pi2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_np1_gfn = PAMR_get_gfn("pi2", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((A_t_n_gfn = PAMR_get_gfn("A_t", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((A_t_np1_gfn = PAMR_get_gfn("A_t", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Atilde_p_n_gfn = PAMR_get_gfn("Atilde_p", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Atilde_p_np1_gfn = PAMR_get_gfn("Atilde_p", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Atilde_phi_n_gfn = PAMR_get_gfn("Atilde_phi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Atilde_phi_np1_gfn = PAMR_get_gfn("Atilde_phi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((A_z_n_gfn = PAMR_get_gfn("A_z", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((A_z_np1_gfn = PAMR_get_gfn("A_z", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((B_t_n_gfn = PAMR_get_gfn("B_t", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((B_t_np1_gfn = PAMR_get_gfn("B_t", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Btilde_p_n_gfn = PAMR_get_gfn("Btilde_p", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Btilde_p_np1_gfn = PAMR_get_gfn("Btilde_p", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Btilde_phi_n_gfn = PAMR_get_gfn("Btilde_phi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Btilde_phi_np1_gfn = PAMR_get_gfn("Btilde_phi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((B_z_n_gfn = PAMR_get_gfn("B_z", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((B_z_np1_gfn = PAMR_get_gfn("B_z", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((gaugecond_n_gfn = PAMR_get_gfn("gaugecond", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((gaugecond_np1_gfn = PAMR_get_gfn("gaugecond", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((divE_n_gfn = PAMR_get_gfn("divE", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((divE_np1_gfn = PAMR_get_gfn("divE", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((divB_n_gfn = PAMR_get_gfn("divB", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((divB_np1_gfn = PAMR_get_gfn("divB", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((chi_n_gfn = PAMR_get_gfn("chi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((chi_np1_gfn = PAMR_get_gfn("chi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((xi_n_gfn = PAMR_get_gfn("xi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_np1_gfn = PAMR_get_gfn("xi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Eden_n_gfn = PAMR_get_gfn("Eden", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Eden_np1_gfn = PAMR_get_gfn("Eden", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Eem_n_gfn = PAMR_get_gfn("Eem", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Eem_np1_gfn = PAMR_get_gfn("Eem", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((Qden_n_gfn = PAMR_get_gfn("Qden", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Qden_np1_gfn = PAMR_get_gfn("Qden", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((elec_p_n_gfn = PAMR_get_gfn("elec_p", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((elec_p_np1_gfn = PAMR_get_gfn("elec_p", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((elec_phi_n_gfn = PAMR_get_gfn("elec_phi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((elec_phi_np1_gfn = PAMR_get_gfn("elec_phi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((elec_z_n_gfn = PAMR_get_gfn("elec_z", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((elec_z_np1_gfn = PAMR_get_gfn("elec_z", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((mag_p_n_gfn = PAMR_get_gfn("mag_p", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((mag_p_np1_gfn = PAMR_get_gfn("mag_p", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((mag_phi_n_gfn = PAMR_get_gfn("mag_phi", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((mag_phi_np1_gfn = PAMR_get_gfn("mag_phi", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((mag_z_n_gfn = PAMR_get_gfn("mag_z", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((mag_z_np1_gfn = PAMR_get_gfn("mag_z", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((dRdp_n_gfn = PAMR_get_gfn("dRdp", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((dRdp2_n_gfn = PAMR_get_gfn("dRdp2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((dXdz_n_gfn = PAMR_get_gfn("dXdz", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((dXdz2_n_gfn = PAMR_get_gfn("dXdz2", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pR_n_gfn = PAMR_get_gfn("pR", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((zX_n_gfn = PAMR_get_gfn("zX", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((chi_res_gfn = PAMR_get_gfn("chi_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((xi_res_gfn = PAMR_get_gfn("xi_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1_res_gfn = PAMR_get_gfn("phi1_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2_res_gfn = PAMR_get_gfn("phi2_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi1_res_gfn = PAMR_get_gfn("pi1_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((pi2_res_gfn = PAMR_get_gfn("pi2_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((A_t_res_gfn = PAMR_get_gfn("A_t_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Atilde_p_res_gfn = PAMR_get_gfn("Atilde_p_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Atilde_phi_res_gfn = PAMR_get_gfn("Atilde_phi_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((A_z_res_gfn = PAMR_get_gfn("A_z_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((B_t_res_gfn = PAMR_get_gfn("B_t_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Btilde_p_res_gfn = PAMR_get_gfn("Btilde_p_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((Btilde_phi_res_gfn = PAMR_get_gfn("Btilde_phi_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((B_z_res_gfn = PAMR_get_gfn("B_z_res", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((phi1ire_n_gfn = PAMR_get_gfn("phi1ire", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi1ire_np1_gfn = PAMR_get_gfn("phi1ire", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2ire_n_gfn = PAMR_get_gfn("phi2ire", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((phi2ire_np1_gfn = PAMR_get_gfn("phi2ire", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((atire_n_gfn = PAMR_get_gfn("atire", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((atire_np1_gfn = PAMR_get_gfn("atire", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((apire_n_gfn = PAMR_get_gfn("apire", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((apire_np1_gfn = PAMR_get_gfn("apire", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((azire_n_gfn = PAMR_get_gfn("azire", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((azire_np1_gfn = PAMR_get_gfn("azire", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((cmask_gfn = PAMR_get_gfn("cmask", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((dummy_n_gfn = PAMR_get_gfn("dummy", PAMR_AMRH, 2)) < 0)
    AMRD_stop("set_gnfs error", 0);

  if ((V_np1_gfn = PAMR_get_gfn("V", PAMR_AMRH, 1)) < 0)
    AMRD_stop("set_gnfs error", 0);
  if ((V_n_gfn = PAMR_get_gfn("V", PAMR_AMRH, 2)) < 0)
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

    if ((dRdp_gfn = PAMR_get_gfn("dRdp", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dRdp2_gfn = PAMR_get_gfn("dRdp2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dXdz_gfn = PAMR_get_gfn("dXdz", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((dXdz2_gfn = PAMR_get_gfn("dXdz2", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((pR_gfn = PAMR_get_gfn("pR", PAMR_MGH, 0)) < 0)
      AMRD_stop("set_gnfs error", 0);
    if ((zX_gfn = PAMR_get_gfn("zX", PAMR_MGH, 0)) < 0)
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
    if ((A_t_gfn = PAMR_get_gfn("A_t", PAMR_MGH, 0)) < 0)
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
 * - Calculates and stores the grid spacings dR, dX
 * - Checks for physical boundaries and computes the grid size
 * - Initializes pointer arrays for coordinates and grid function data
 *
 **********************************************************************/
void ldptr(void) {
  real dx0[2];              // stores grid spacings dR, dX
  real *x0[2];              // stores pointers to coordinates R, X
  real *gfs[PAMR_MAX_GFNS]; // stores pointers to grid function data

  static int first = 1; // only initialize grid functions once
  if (first) {
    first = 0;
    set_gfns();

    // initialize coordinate bounding box of the base grid
    // (e.g. [Rmin, Rmax, Xmin, Xmax])
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
  dR = dx0[0];
  dX = dx0[1];

  // check for AMR and physical boundaries and compute grid size
  if ((bbox[0] - base_bbox[0]) < dR / 2)
    phys_bdy[0] = 1;
  else
    phys_bdy[0] = 0;
  if ((base_bbox[1] - bbox[1]) < dR / 2)
    phys_bdy[1] = 1;
  else
    phys_bdy[1] = 0;
  NR = shape[0];
  size = NR;
  NX = 1;
  if ((bbox[2] - base_bbox[2]) < dX / 2)
    phys_bdy[2] = 1;
  else
    phys_bdy[2] = 0;
  if ((base_bbox[3] - bbox[3]) < dX / 2)
    phys_bdy[3] = 1;
  else
    phys_bdy[3] = 0;
  NX = shape[1];
  size *= NX;

  PAMR_get_g_x(x0); // retrieve pointers to arrays containing coordinates
  R = x0[0];
  X = x0[1];

  // retrieve pointers to arrays containing grid function data
  PAMR_get_g_gfs(gfs);
  modphi_n = gfs[modphi_n_gfn - 1];
  modphi_np1 = gfs[modphi_np1_gfn - 1];
  phi1_n = gfs[phi1_n_gfn - 1];
  phi1_np1 = gfs[phi1_np1_gfn - 1];
  phi2_n = gfs[phi2_n_gfn - 1];
  phi2_np1 = gfs[phi2_np1_gfn - 1];
  pi1_n = gfs[pi1_n_gfn - 1];
  pi1_np1 = gfs[pi1_np1_gfn - 1];
  pi2_n = gfs[pi2_n_gfn - 1];
  pi2_np1 = gfs[pi2_np1_gfn - 1];
  A_t_n = gfs[A_t_n_gfn - 1];
  A_t_np1 = gfs[A_t_np1_gfn - 1];
  Atilde_p_n = gfs[Atilde_p_n_gfn - 1];
  Atilde_p_np1 = gfs[Atilde_p_np1_gfn - 1];
  Atilde_phi_n = gfs[Atilde_phi_n_gfn - 1];
  Atilde_phi_np1 = gfs[Atilde_phi_np1_gfn - 1];
  A_z_n = gfs[A_z_n_gfn - 1];
  A_z_np1 = gfs[A_z_np1_gfn - 1];
  B_t_n = gfs[B_t_n_gfn - 1];
  B_t_np1 = gfs[B_t_np1_gfn - 1];
  Btilde_p_n = gfs[Btilde_p_n_gfn - 1];
  Btilde_p_np1 = gfs[Btilde_p_np1_gfn - 1];
  Btilde_phi_n = gfs[Btilde_phi_n_gfn - 1];
  Btilde_phi_np1 = gfs[Btilde_phi_np1_gfn - 1];
  B_z_n = gfs[B_z_n_gfn - 1];
  B_z_np1 = gfs[B_z_np1_gfn - 1];
  gaugecond_n = gfs[gaugecond_n_gfn - 1];
  gaugecond_np1 = gfs[gaugecond_np1_gfn - 1];
  divE_n = gfs[divE_n_gfn - 1];
  divE_np1 = gfs[divE_np1_gfn - 1];
  divB_n = gfs[divB_n_gfn - 1];
  divB_np1 = gfs[divB_np1_gfn - 1];
  chi_n = gfs[chi_n_gfn - 1];
  chi_np1 = gfs[chi_np1_gfn - 1];
  xi_n = gfs[xi_n_gfn - 1];
  xi_np1 = gfs[xi_np1_gfn - 1];
  Eden_n = gfs[Eden_n_gfn - 1];
  Eden_np1 = gfs[Eden_np1_gfn - 1];
  Eem_n = gfs[Eem_n_gfn - 1];
  Eem_np1 = gfs[Eem_np1_gfn - 1];
  Qden_n = gfs[Qden_n_gfn - 1];
  Qden_np1 = gfs[Qden_np1_gfn - 1];
  elec_p_n = gfs[elec_p_n_gfn - 1];
  elec_p_np1 = gfs[elec_p_np1_gfn - 1];
  elec_phi_n = gfs[elec_phi_n_gfn - 1];
  elec_phi_np1 = gfs[elec_phi_np1_gfn - 1];
  elec_z_n = gfs[elec_z_n_gfn - 1];
  elec_z_np1 = gfs[elec_z_np1_gfn - 1];
  mag_p_n = gfs[mag_p_n_gfn - 1];
  mag_p_np1 = gfs[mag_p_np1_gfn - 1];
  mag_phi_n = gfs[mag_phi_n_gfn - 1];
  mag_phi_np1 = gfs[mag_phi_np1_gfn - 1];
  mag_z_n = gfs[mag_z_n_gfn - 1];
  mag_z_np1 = gfs[mag_z_np1_gfn - 1];

  dRdp_n = gfs[dRdp_n_gfn - 1];
  dRdp2_n = gfs[dRdp2_n_gfn - 1];
  dXdz_n = gfs[dXdz_n_gfn - 1];
  dXdz2_n = gfs[dXdz2_n_gfn - 1];
  pR_n = gfs[pR_n_gfn - 1];
  zX_n = gfs[zX_n_gfn - 1];

  chi_res = gfs[chi_res_gfn - 1];
  xi_res = gfs[xi_res_gfn - 1];
  phi1_res = gfs[phi1_res_gfn - 1];
  phi2_res = gfs[phi2_res_gfn - 1];
  pi1_res = gfs[pi1_res_gfn - 1];
  pi2_res = gfs[pi2_res_gfn - 1];
  A_t_res = gfs[A_t_res_gfn - 1];
  Atilde_p_res = gfs[Atilde_p_res_gfn - 1];
  Atilde_phi_res = gfs[Atilde_phi_res_gfn - 1];
  A_z_res = gfs[A_z_res_gfn - 1];
  B_t_res = gfs[B_t_res_gfn - 1];
  Btilde_p_res = gfs[Btilde_p_res_gfn - 1];
  Btilde_phi_res = gfs[Btilde_phi_res_gfn - 1];
  B_z_res = gfs[B_z_res_gfn - 1];

  phi1ire_n = gfs[phi1ire_n_gfn - 1];
  phi1ire_np1 = gfs[phi1ire_np1_gfn - 1];
  phi2ire_n = gfs[phi2ire_n_gfn - 1];
  phi2ire_np1 = gfs[phi2ire_np1_gfn - 1];
  atire_n = gfs[atire_n_gfn - 1];
  atire_np1 = gfs[atire_np1_gfn - 1];
  apire_n = gfs[apire_n_gfn - 1];
  apire_np1 = gfs[apire_np1_gfn - 1];
  azire_n = gfs[azire_n_gfn - 1];
  azire_np1 = gfs[azire_np1_gfn - 1];

  cmask = gfs[cmask_gfn - 1];

  dummy_n = gfs[dummy_n_gfn - 1];

  V_n = gfs[V_n_gfn - 1];
  V_np1 = gfs[V_np1_gfn - 1];

  if (mg_toggle == 1) {
    V = gfs[V_gfn - 1];

    mask_mg = gfs[mask_mg_gfn - 1];

    V_res = gfs[V_res_gfn - 1];
    V_rhs = gfs[V_rhs_gfn - 1];
    V_lop = gfs[V_lop_gfn - 1];

    dRdp = gfs[dRdp_gfn - 1];
    dRdp2 = gfs[dRdp2_gfn - 1];
    dXdz = gfs[dXdz_gfn - 1];
    dXdz2 = gfs[dXdz2_gfn - 1];
    pR = gfs[pR_gfn - 1];
    zX = gfs[zX_gfn - 1];

    phi1 = gfs[phi1_gfn - 1];
    phi2 = gfs[phi2_gfn - 1];
    pi1 = gfs[pi1_gfn - 1];
    pi2 = gfs[pi2_gfn - 1];
    modphi = gfs[modphi_gfn - 1];
    A_t = gfs[A_t_gfn - 1];
  }

  // recompute grid spacing
  R = x0[0];
  dR = R[1] - R[0];
  X = x0[1];
  dX = X[1] - X[0];
}

/***********************************************************************
 * const_f
 *
 * Sets a grid function to a constant value.
 *
 ***********************************************************************/
void const_f(real *f, real c) {
  int i;

  for (i = 0; i < NR * NX; i++)
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
  char hostname[256]; // name of the processor

  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready\n", getpid(), hostname);
  fflush(stdout);
  sleep(1);

  //// uncomment if you want to attach debugger at the start of execution
  //volatile int i = 0;
  //while (0 == i) {sleep(5);} // sleeping won't use 100% CPU

  // initialize multigrid functionality
  mg_toggle = 0;
  AMRD_int_param(pfile, "mg_toggle", &mg_toggle, 1);
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
  AMRD_real_param(pfile, "c", &c, 1);
  AMRD_real_param(pfile, "d", &d, 1);
  AMRD_real_param(pfile, "cR", &cR, 1);
  AMRD_real_param(pfile, "cR1", &cR1, 1);
  AMRD_real_param(pfile, "cR2", &cR2, 1);
  AMRD_real_param(pfile, "cX", &cX, 1);
  AMRD_real_param(pfile, "cX1", &cX1, 1);
  AMRD_real_param(pfile, "cX2", &cX2, 1);
  AMRD_real_param(pfile, "e", &e, 1);
  AMRD_real_param(pfile, "w", &w, 1);
  AMRD_real_param(pfile, "beta", &beta, 1);
  AMRD_real_param(pfile, "mu", &mu, 1);
  AMRD_real_param(pfile, "h", &h, 1);
  AMRD_real_param(pfile, "mm", &mm, 1);
  AMRD_real_param(pfile, "gg", &gg, 1);
  AMRD_real_param(pfile, "epsdis", &epsdis, 1);
  AMRD_real_param(pfile, "amp", &amp, 1);
  AMRD_real_param(pfile, "delta", &delta, 1);
  AMRD_real_param(pfile, "R_center", &R_center, 1);
  AMRD_real_param(pfile, "Rw", &Rw, 1);
  AMRD_real_param(pfile, "r0", &r0, 1);
  AMRD_real_param(pfile, "signum", &signum, 1);
  AMRD_real_param(pfile, "X_center", &X_center, 1);
  AMRD_real_param(pfile, "Xw", &Xw, 1);
  AMRD_real_param(pfile, "cc", &cc, 1);
  AMRD_int_param(pfile, "Vtype", &Vtype, 1);
  AMRD_int_param_v(pfile, "base_shape", &gridsize, &dim);

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

  zero(modphi_n);
  zero(modphi_np1);
  zero(phi1_n);
  zero(phi1_np1);
  zero(phi2_n);
  zero(phi2_np1);
  zero(pi1_n);
  zero(pi1_np1);
  zero(pi2_n);
  zero(pi2_np1);
  zero(A_t_n);
  zero(A_t_np1);
  zero(Atilde_p_n);
  zero(Atilde_p_np1);
  zero(Atilde_phi_n);
  zero(Atilde_phi_np1);
  zero(A_z_n);
  zero(A_z_np1);
  zero(B_t_n);
  zero(B_t_np1);
  zero(Btilde_p_n);
  zero(Btilde_p_np1);
  zero(Btilde_phi_n);
  zero(Btilde_phi_np1);
  zero(B_z_n);
  zero(B_z_np1);
  zero(gaugecond_n);
  zero(gaugecond_np1);
  zero(divE_n);
  zero(divE_np1);
  zero(divB_n);
  zero(divB_np1);
  zero(chi_n);
  zero(chi_np1);
  zero(xi_n);
  zero(xi_np1);
  zero(Eden_n);
  zero(Eden_np1);
  zero(Eem_n);
  zero(Eem_np1);
  zero(Qden_n);
  zero(Qden_np1);
  zero(elec_p_n);
  zero(elec_p_np1);
  zero(elec_phi_n);
  zero(elec_phi_np1);
  zero(elec_z_n);
  zero(elec_z_np1);
  zero(mag_p_n);
  zero(mag_p_np1);
  zero(mag_phi_n);
  zero(mag_phi_np1);
  zero(mag_z_n);
  zero(mag_z_np1);

  zero(dRdp_n);
  zero(dRdp2_n);
  zero(dXdz_n);
  zero(dXdz2_n);
  zero(pR_n);
  zero(zX_n);

  zero(chi_res);
  zero(xi_res);
  zero(phi1_res);
  zero(phi2_res);
  zero(pi1_res);
  zero(pi2_res);
  zero(A_t_res);
  zero(Atilde_p_res);
  zero(Atilde_phi_res);
  zero(A_z_res);
  zero(B_t_res);
  zero(Btilde_p_res);
  zero(Btilde_phi_res);
  zero(B_z_res);

  zero(phi1ire_n);
  zero(phi1ire_np1);
  zero(phi2ire_n);
  zero(phi2ire_np1);
  zero(atire_n);
  zero(atire_np1);
  zero(apire_n);
  zero(apire_np1);
  zero(azire_n);
  zero(azire_np1);

  zero(cmask);

  zero(dummy_n);

  zero(V_n);
  zero(V_np1);
}

/***********************************************************************
 * qball_free_data
 *
 * Generates initial data for all relevant grid functions.
 *
 ***********************************************************************/
void qball_free_data(void) {
  ldptr();

  initializer0_(dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, R, X,
                &c, &d);

  init_qball_(phi1_np1, phi1_n, phi2_np1, phi2_n, pi1_np1, pi1_n, pi2_np1,
              pi2_n, A_t_np1, A_t_n, Atilde_p_np1, Atilde_p_n, Atilde_phi_np1,
              Atilde_phi_n, A_z_np1, A_z_n, B_t_np1, B_t_n, Btilde_p_np1,
              Btilde_p_n, Btilde_phi_np1, Btilde_phi_n, B_z_np1, B_z_n, dRdp_n,
              dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, R, X, &dR, &dX,
              &c, &cR, &cR1, &cR2, &cX, &cX1, &cX2, &d, &e, &base_bbox[1],
              &base_bbox[0], &w, &base_bbox[3], &base_bbox[2]);

  initializer1_(phi1ire_n, phi2ire_n, atire_n, apire_n, azire_n, &NR, &NX);

  initializer3_(modphi_n, phi1_n, phi2_n, pi1_n, pi2_n, A_t_n, Eden_n, Eem_n,
                Qden_n, elec_p_n, elec_phi_n, elec_z_n, mag_p_n, mag_phi_n,
                mag_z_n, gaugecond_n, divE_n, divB_n, chi_n, xi_n, dummy_n,
                dRdp_n, dXdz_n, pR_n, zX_n, &NR, &NX, &dR, &dX, &amp, &delta,
                &e, &r0, &R_center, &Rw, &signum, &X_center, &Xw);

  // if multigrid data file exists, read it and evolve (overwriting previous data)
  // if multigrid data file doesn't exist, create it
  if (access(mg_data, F_OK) == 0) {
    init_qball2_(V_n, R, X, &NR, &NX, pR_n, zX_n, &dR, &dX, dRdp_n, dRdp2_n,
                 dXdz_n, dXdz2_n, A_t_n, Atilde_p_n, A_z_n, B_t_n, Btilde_p_n,
                 B_z_n, &base_bbox[0], &base_bbox[1], &base_bbox[2],
                 &base_bbox[3]);
  } else if (mg_toggle == 1) {
    initguess_(V_n, R, X, &NR, &NX, pR_n, zX_n, A_t_n);
  }
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
 *
 **********************************************************************/
real qball_evo_residual(void) {
  ldptr();

  real l2norm = 0.0;
  real l2norm_chi = 0.0;
  real l2norm_xi = 0.0;
  real l2norm_phi1 = 0.0;
  real l2norm_phi2 = 0.0;
  real l2norm_pi1 = 0.0;
  real l2norm_pi2 = 0.0;
  real l2norm_A_t = 0.0;
  real l2norm_Atilde_p = 0.0;
  real l2norm_Atilde_phi = 0.0;
  real l2norm_A_z = 0.0;
  real l2norm_B_t = 0.0;
  real l2norm_Btilde_p = 0.0;
  real l2norm_Btilde_phi = 0.0;
  real l2norm_B_z = 0.0;

  res_chi_(chi_res, chi_np1, chi_n, xi_np1, xi_n, &NR, &NX, &dt);

  res_xi_(xi_res, phi1_np1, phi1_n, phi2_np1, phi2_n, chi_np1, chi_n, xi_np1,
          xi_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, &dR,
          &dt, &dX, &cc, &phys_bdy[1], &phys_bdy[0], &phys_bdy[3],
          &phys_bdy[2]);

  res_phi1_(phi1_res, phi1_np1, phi1_n, pi1_np1, pi1_n, &NR, &NX, &dt);

  res_phi2_(phi2_res, phi2_np1, phi2_n, pi2_np1, pi2_n, &NR, &NX, &dt);

  res_pi1_(pi1_res, modphi_np1, modphi_n, phi1_np1, phi1_n, phi2_np1, phi2_n,
           pi1_np1, pi1_n, pi2_np1, pi2_n, A_t_np1, A_t_n, Atilde_p_np1,
           Atilde_p_n, Atilde_phi_np1, Atilde_phi_n, A_z_np1, A_z_n, chi_np1,
           chi_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, &dR,
           &dt, &dX, &beta, &cc, &e, &gg, &h, &mm, &mu, &phys_bdy[1],
           &phys_bdy[0], &Vtype, &phys_bdy[3], &phys_bdy[2]);

  res_pi2_(pi2_res, modphi_np1, modphi_n, phi1_np1, phi1_n, phi2_np1, phi2_n,
           pi1_np1, pi1_n, pi2_np1, pi2_n, A_t_np1, A_t_n, Atilde_p_np1,
           Atilde_p_n, Atilde_phi_np1, Atilde_phi_n, A_z_np1, A_z_n, chi_np1,
           chi_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, &dR,
           &dt, &dX, &beta, &cc, &e, &gg, &h, &mm, &mu, &phys_bdy[1],
           &phys_bdy[0], &Vtype, &phys_bdy[3], &phys_bdy[2]);

  res_a_t_(A_t_res, A_t_np1, A_t_n, B_t_np1, B_t_n, &NR, &NX, &dt);

  res_atilde_p_(Atilde_p_res, Atilde_p_np1, Atilde_p_n, Btilde_p_np1,
                Btilde_p_n, &NR, &NX, &dt);

  res_atilde_phi_(Atilde_phi_res, Atilde_phi_np1, Atilde_phi_n, Btilde_phi_np1,
                  Btilde_phi_n, &NR, &NX, &dt);

  res_a_z_(A_z_res, A_z_np1, A_z_n, B_z_np1, B_z_n, &NR, &NX, &dt);

  res_b_t_(B_t_res, phi1_np1, phi1_n, phi2_np1, phi2_n, pi1_np1, pi1_n, pi2_np1,
           pi2_n, A_t_np1, A_t_n, B_t_np1, B_t_n, dRdp_n, dRdp2_n, dXdz_n,
           dXdz2_n, pR_n, zX_n, &NR, &NX, &dR, &dt, &dX, &e, &phys_bdy[1],
           &phys_bdy[0], &phys_bdy[3], &phys_bdy[2]);

  res_btilde_p_(Btilde_p_res, phi1_np1, phi1_n, phi2_np1, phi2_n, Atilde_p_np1,
                Atilde_p_n, Btilde_p_np1, Btilde_p_n, dRdp_n, dRdp2_n, dXdz_n,
                dXdz2_n, pR_n, zX_n, &NR, &NX, &dR, &dt, &dX, &e, &phys_bdy[1],
                &phys_bdy[0], &phys_bdy[3], &phys_bdy[2]);

  res_btilde_phi_(Btilde_phi_res, phi1_np1, phi1_n, phi2_np1, phi2_n,
                  Atilde_phi_np1, Atilde_phi_n, Btilde_phi_np1, Btilde_phi_n,
                  dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, &dR,
                  &dt, &dX, &e, &phys_bdy[1], &phys_bdy[0], &phys_bdy[3],
                  &phys_bdy[2]);

  res_b_z_(B_z_res, phi1_np1, phi1_n, phi2_np1, phi2_n, A_z_np1, A_z_n, B_z_np1,
           B_z_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, &dR,
           &dt, &dX, &e, &phys_bdy[1], &phys_bdy[0], &phys_bdy[3],
           &phys_bdy[2]);

  l2norm_chi = l2norm_calc(chi_res);
  l2norm_xi = l2norm_calc(xi_res);
  l2norm_phi1 = l2norm_calc(phi1_res);
  l2norm_phi2 = l2norm_calc(phi2_res);
  l2norm_pi1 = l2norm_calc(pi1_res);
  l2norm_pi2 = l2norm_calc(pi2_res);
  l2norm_A_t = l2norm_calc(A_t_res);
  l2norm_Atilde_p = l2norm_calc(Atilde_p_res);
  l2norm_Atilde_phi = l2norm_calc(Atilde_phi_res);
  l2norm_A_z = l2norm_calc(A_z_res);
  l2norm_B_t = l2norm_calc(B_t_res);
  l2norm_Btilde_p = l2norm_calc(Btilde_p_res);
  l2norm_Btilde_phi = l2norm_calc(Btilde_phi_res);
  l2norm_B_z = l2norm_calc(B_z_res);

  if (g_norms[chi_n_gfn - 1] > 0) {
    l2norm += l2norm_chi / g_norms[chi_n_gfn - 1];
  }
  if (g_norms[xi_n_gfn - 1] > 0) {
    l2norm += l2norm_xi / g_norms[xi_n_gfn - 1];
  }
  if (g_norms[phi1_n_gfn - 1] > 0) {
    l2norm += l2norm_phi1 / g_norms[phi1_n_gfn - 1];
  }
  if (g_norms[phi2_n_gfn - 1] > 0) {
    l2norm += l2norm_phi2 / g_norms[phi2_n_gfn - 1];
  }
  if (g_norms[pi1_n_gfn - 1] > 0) {
    l2norm += l2norm_pi1 / g_norms[pi1_n_gfn - 1];
  }
  if (g_norms[pi2_n_gfn - 1] > 0) {
    l2norm += l2norm_pi2 / g_norms[pi2_n_gfn - 1];
  }
  if (g_norms[A_t_n_gfn - 1] > 0) {
    l2norm += l2norm_A_t / g_norms[A_t_n_gfn - 1];
  }
  if (g_norms[Atilde_p_n_gfn - 1] > 0) {
    l2norm += l2norm_Atilde_p / g_norms[Atilde_p_n_gfn - 1];
  }
  if (g_norms[Atilde_phi_n_gfn - 1] > 0) {
    l2norm += l2norm_Atilde_phi / g_norms[Atilde_phi_n_gfn - 1];
  }
  if (g_norms[A_z_n_gfn - 1] > 0) {
    l2norm += l2norm_A_z / g_norms[A_z_n_gfn - 1];
  }
  if (g_norms[B_t_n_gfn - 1] > 0) {
    l2norm += l2norm_B_t / g_norms[B_t_n_gfn - 1];
  }
  if (g_norms[Btilde_p_n_gfn - 1] > 0) {
    l2norm += l2norm_Btilde_p / g_norms[Btilde_p_n_gfn - 1];
  }
  if (g_norms[Btilde_phi_n_gfn - 1] > 0) {
    l2norm += l2norm_Btilde_phi / g_norms[Btilde_phi_n_gfn - 1];
  }
  if (g_norms[B_z_n_gfn - 1] > 0) {
    l2norm += l2norm_B_z / g_norms[B_z_n_gfn - 1];
  }

  // verbose output
  if (0 && my_rank == 0) {
    printf("\n - - - - - - - - - - - - - - - - - - - \n");

    printf("rank %d, level %d\n", my_rank, g_L);

    if (1) // set to 1 or 0 to toggle
    {
      printf("\nGRID FUNCTION L-INFINITY NORMS:\n");
      printf("linfnorm_chi=        %15.13g\n", g_norms[chi_n_gfn - 1]);
      printf("linfnorm_xi=         %15.13g\n", g_norms[xi_n_gfn - 1]);
      printf("linfnorm_phi1=       %15.13g\n", g_norms[phi1_n_gfn - 1]);
      printf("linfnorm_phi2=       %15.13g\n", g_norms[phi2_n_gfn - 1]);
      printf("linfnorm_pi1=        %15.13g\n", g_norms[pi1_n_gfn - 1]);
      printf("linfnorm_pi2=        %15.13g\n", g_norms[pi2_n_gfn - 1]);
      printf("linfnorm_A_t=        %15.13g\n", g_norms[A_t_n_gfn - 1]);
      printf("linfnorm_Atilde_p=   %15.13g\n", g_norms[Atilde_p_n_gfn - 1]);
      printf("linfnorm_Atilde_phi= %15.13g\n", g_norms[Atilde_phi_n_gfn - 1]);
      printf("linfnorm_A_z=        %15.13g\n", g_norms[A_z_n_gfn - 1]);
      printf("linfnorm_B_t=        %15.13g\n", g_norms[B_t_n_gfn - 1]);
      printf("linfnorm_Btilde_p=   %15.13g\n", g_norms[Btilde_p_n_gfn - 1]);
      printf("linfnorm_Btilde_phi= %15.13g\n", g_norms[Btilde_phi_n_gfn - 1]);
      printf("linfnorm_B_z=        %15.13g\n", g_norms[B_z_n_gfn - 1]);
    }

    if (1) // set to 1 or 0 to toggle
    {
      printf("\nGRID FUNCTION L2-NORMS:\n");
      printf("l2norm_chi= %15.13g\n", l2norm_chi);
      printf("l2norm_xi= %15.13g\n", l2norm_xi);
      printf("l2norm_phi1= %15.13g\n", l2norm_phi1);
      printf("l2norm_phi2= %15.13g\n", l2norm_phi2);
      printf("l2norm_pi1=  %15.13g\n", l2norm_pi1);
      printf("l2norm_pi2=  %15.13g\n", l2norm_pi2);
      printf("l2norm_A_t=%15.13g\n", l2norm_A_t);
      printf("l2norm_Atilde_p=%15.13g\n", l2norm_Atilde_p);
      printf("l2norm_Atilde_phi=%15.13g\n", l2norm_Atilde_phi);
      printf("l2norm_A_z=%15.13g\n", l2norm_A_z);
      printf("l2norm_B_t=%15.13g\n", l2norm_B_t);
      printf("l2norm_Btilde_p=%15.13g\n", l2norm_Btilde_p);
      printf("l2norm_Btilde_phi=%15.13g\n", l2norm_Btilde_phi);
      printf("l2norm_B_z=%15.13g\n", l2norm_B_z);
    }

    if (1) // set to 1 or 0 to toggle
    {
      printf("\nTOTAL (NORMALIZED) L2-NORM:\nl2norm= %15.13g\n", l2norm);
    }
  }

  // stop evolution if L2-norm becomes NaN or inf
  if (l2norm != l2norm) {
    AMRD_stop("NaN detected in app_evo_residual()... AMRD_stop", 0);
  }

  return l2norm;
}

/***********************************************************************
 * qball_evolve
 *
 * Performs 1 iteration of the evolution equations.
 *
 **********************************************************************/
void qball_evolve(int iter, int *ifc_mask) {

  // if multigrid data file doesn't exist but mg_toggle turned on, write it and exit
  if ((access(mg_data, F_OK) != 0) && (mg_toggle == 1)) {
    w2f_(V_np1, &NR, &NX); // write to binary file
    AMRD_stop("Multigrid processing complete.", 0);
  }

  ldptr();

  initializer0_(dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, R, X,
                &c, &d);

  update0_(modphi_np1, modphi_n, phi1_np1, phi1_n, phi2_np1, phi2_n, pi1_np1,
           pi1_n, pi2_np1, pi2_n, A_t_np1, A_t_n, Atilde_p_np1, Atilde_p_n,
           Atilde_phi_np1, Atilde_phi_n, A_z_np1, A_z_n, B_t_np1, B_t_n,
           Btilde_p_np1, Btilde_p_n, Btilde_phi_np1, Btilde_phi_n, B_z_np1,
           B_z_n, Eem_np1, Qden_np1, elec_p_np1, elec_phi_np1, elec_z_np1,
           mag_p_np1, mag_phi_np1, mag_z_np1, gaugecond_np1, divE_np1, divB_np1,
           chi_np1, chi_n, xi_np1, xi_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n, pR_n,
           zX_n, &NR, &NX, &dR, &dt, &dX, &beta, &cc, &e, &gg, &h, &mm, &mu,
           &phys_bdy[1], &phys_bdy[0], &Vtype, &phys_bdy[3], &phys_bdy[2]);

  edencalc_(modphi_np1, modphi_n, phi1_np1, phi1_n, phi2_np1, phi2_n, pi1_np1,
            pi1_n, pi2_np1, pi2_n, A_t_np1, A_t_n, Atilde_p_np1, Atilde_p_n,
            Atilde_phi_np1, Atilde_phi_n, A_z_np1, A_z_n, B_t_np1, B_t_n,
            Btilde_p_np1, Btilde_p_n, Btilde_phi_np1, Btilde_phi_n, B_z_np1,
            B_z_n, Eden_np1, Eden_n, chi_np1, chi_n, xi_np1, xi_n, dRdp_n,
            dRdp2_n, dXdz_n, dXdz2_n, pR_n, zX_n, &NR, &NX, R, X, &dR, &dt, &dX,
            &beta, &c, &cc, &d, &e, &gg, &h, &mm, &mu, &Vtype);

  ire_(modphi_np1, modphi_n, phi1_np1, phi1_n, phi2_np1, phi2_n, pi1_np1, pi1_n,
       pi2_np1, pi2_n, A_t_np1, A_t_n, Atilde_p_np1, Atilde_p_n, Atilde_phi_np1,
       Atilde_phi_n, A_z_np1, A_z_n, B_t_np1, B_t_n, Btilde_p_np1, Btilde_p_n,
       Btilde_phi_np1, Btilde_phi_n, B_z_np1, B_z_n, chi_np1, chi_n, xi_np1,
       xi_n, phi1ire_np1, phi1ire_n, phi2ire_np1, phi2ire_n, atire_np1, atire_n,
       apire_np1, apire_n, azire_np1, azire_n, dRdp_n, dRdp2_n, dXdz_n, dXdz2_n,
       pR_n, zX_n, &NR, &NX, R, X, &dR, &dt, &dX, &beta, &c, &cc, &d, &e, &gg,
       &h, &mm, &mu, &Vtype);
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
  int lshape[2];             // number of child grid points per dimension
  int lghost_width[4];       // width of child AMR ghost region per dimension
  int procgs, totgs;         // grid size
  double lbbox[4];           // stores child AMR bounding box data
  real ldR, ldX, ldt;        // child grid spacing in each dimension
  real ldx0[2];              // stores child grid spacings ldR, ldX
  real *lx0[2];              // stores pointers to coordinates x, y, z
  real *gfs0[PAMR_MAX_GFNS]; // stores pointers to child grid function data

  // AMR integrated quantities (charge, energy)
  double ltotalQ = 0.0, proctotalQ = 0.0, grandtotalQ = 0.0;
  double ltotalE = 0.0, proctotalE = 0.0, grandtotalE = 0.0;
  double ltotalEem = 0.0, proctotalEem = 0.0, grandtotalEem = 0.0;

  // AMR fractional integrated quantities (charge)
  double ltotalQhalf_upper = 0.0, proctotalQhalf_upper = 0.0,
         grandtotalQhalf_upper = 0.0;
  double ltotalQhalf_lower = 0.0, proctotalQhalf_lower = 0.0,
         grandtotalQhalf_lower = 0.0;

  // unigrid integrated quantities (base level)
  double total_Qbase = 0.0, proctotal_Qbase = 0.0;
  double total_Ebase = 0.0, proctotal_Ebase = 0.0;
  double total_Eembase = 0.0, proctotal_Eembase = 0.0;
  double total_phi1IRE = 0.0, total_phi2IRE = 0.0;
  double total_atIRE = 0.0, total_apIRE = 0.0, total_azIRE = 0.0;
  double proctotal_phi1IRE = 0.0, proctotal_phi2IRE = 0.0;
  double proctotal_atIRE = 0.0, proctotal_apIRE = 0.0, proctotal_azIRE = 0.0;
  double total_gaugecond = 0.0, proctotal_gaugecond = 0.0;
  double total_divE = 0.0, proctotal_divE = 0.0;

  // pointers for data output
  FILE *fileQtot, *fileEtot, *fileEemtot;
  FILE *fileQtothalf_upper, *fileQtothalf_lower;
  FILE *file_Qtotbase, *file_Etotbase, *file_Eemtotbase;
  FILE *file_phi1IRE, *file_phi2IRE;
  FILE *file_atIRE, *file_apIRE, *file_azIRE;
  FILE *file_gaugecond, *file_divE;

  // pointers for receiving MPI data
  int num_procs;
  real *recvQ = NULL, *recvE = NULL, *recvEem = NULL;
  real *recvQhalf_upper = NULL, *recvQhalf_lower = NULL;
  real *recv_Qbase = NULL, *recv_Ebase = NULL, *recv_Eembase = NULL;
  real *recv_phi1IRE = NULL, *recv_phi2IRE = NULL;
  real *recv_atIRE = NULL, *recv_apIRE = NULL, *recv_azIRE = NULL;
  real *recv_gaugecond = NULL, *recv_divE = NULL;

  // retrieve number of processors and allocate storage for receive buffer
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  if (my_rank == 0) {
    recvQ = (double *)calloc(num_procs, sizeof(double));
    recvE = (double *)calloc(num_procs, sizeof(double));
    recvEem = (double *)calloc(num_procs, sizeof(double));
    recvQhalf_upper = (double *)calloc(num_procs, sizeof(double));
    recvQhalf_lower = (double *)calloc(num_procs, sizeof(double));
    recv_Qbase = (double *)calloc(num_procs, sizeof(double));
    recv_Ebase = (double *)calloc(num_procs, sizeof(double));
    recv_Eembase = (double *)calloc(num_procs, sizeof(double));
    recv_phi1IRE = (double *)calloc(num_procs, sizeof(double));
    recv_phi2IRE = (double *)calloc(num_procs, sizeof(double));
    recv_atIRE = (double *)calloc(num_procs, sizeof(double));
    recv_apIRE = (double *)calloc(num_procs, sizeof(double));
    recv_azIRE = (double *)calloc(num_procs, sizeof(double));
    recv_gaugecond = (double *)calloc(num_procs, sizeof(double));
    recv_divE = (double *)calloc(num_procs, sizeof(double));
  }

  // retrieve max AMR level in the hierarchy
  maxL = PAMR_get_max_lev(PAMR_AMRH);

  // if using AMR, compute base quantities on level 2 (always fully resolved)
  if (maxL == 1) {
    loopL = 1;
    totgs = gridsize[0] * gridsize[1];
  } else {
    loopL = 2;
    totgs = 4 * gridsize[0] * gridsize[1];
  }

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
      ldR = ldx0[0];
      ldX = ldx0[1];

      PAMR_get_g_x(lx0); // retrieve pointers to arrays containing coordinates

      // retrieve pointers to arrays containing grid function data
      PAMR_get_g_gfs(gfs0);

      modphi_n = gfs0[modphi_np1_gfn - 1];
      modphi_np1 = gfs0[modphi_n_gfn - 1];

      phi1_n = gfs0[phi1_np1_gfn - 1];
      phi1_np1 = gfs0[phi1_n_gfn - 1];

      phi2_n = gfs0[phi2_np1_gfn - 1];
      phi2_np1 = gfs0[phi2_n_gfn - 1];

      pi1_n = gfs0[pi1_np1_gfn - 1];
      pi1_np1 = gfs0[pi1_n_gfn - 1];

      pi2_n = gfs0[pi2_np1_gfn - 1];
      pi2_np1 = gfs0[pi2_n_gfn - 1];

      A_t_n = gfs0[A_t_np1_gfn - 1];
      A_t_np1 = gfs0[A_t_n_gfn - 1];

      Atilde_p_n = gfs0[Atilde_p_np1_gfn - 1];
      Atilde_p_np1 = gfs0[Atilde_p_n_gfn - 1];

      Atilde_phi_n = gfs0[Atilde_phi_np1_gfn - 1];
      Atilde_phi_np1 = gfs0[Atilde_phi_n_gfn - 1];

      A_z_n = gfs0[A_z_np1_gfn - 1];
      A_z_np1 = gfs0[A_z_n_gfn - 1];

      B_t_n = gfs0[B_t_np1_gfn - 1];
      B_t_np1 = gfs0[B_t_n_gfn - 1];

      Btilde_p_n = gfs0[Btilde_p_np1_gfn - 1];
      Btilde_p_np1 = gfs0[Btilde_p_n_gfn - 1];

      Btilde_phi_n = gfs0[Btilde_phi_np1_gfn - 1];
      Btilde_phi_np1 = gfs0[Btilde_phi_n_gfn - 1];

      B_z_n = gfs0[B_z_np1_gfn - 1];
      B_z_np1 = gfs0[B_z_n_gfn - 1];

      dRdp_n = gfs0[dRdp_n_gfn - 1];
      dRdp2_n = gfs0[dRdp2_n_gfn - 1];

      dXdz_n = gfs0[dXdz_n_gfn - 1];
      dXdz2_n = gfs0[dXdz2_n_gfn - 1];

      pR_n = gfs0[pR_n_gfn - 1];
      zX_n = gfs0[zX_n_gfn - 1];

      gaugecond_n = gfs0[gaugecond_np1_gfn - 1];
      gaugecond_np1 = gfs0[gaugecond_n_gfn - 1];

      divE_n = gfs0[divE_np1_gfn - 1];
      divE_np1 = gfs0[divE_n_gfn - 1];

      divB_n = gfs0[divB_np1_gfn - 1];
      divB_np1 = gfs0[divB_n_gfn - 1];

      chi_n = gfs0[chi_np1_gfn - 1];
      chi_np1 = gfs0[chi_n_gfn - 1];

      xi_n = gfs0[xi_np1_gfn - 1];
      xi_np1 = gfs0[xi_n_gfn - 1];

      Eden_n = gfs0[Eden_np1_gfn - 1];
      Eden_np1 = gfs0[Eden_n_gfn - 1];

      Eem_n = gfs0[Eem_np1_gfn - 1];
      Eem_np1 = gfs0[Eem_n_gfn - 1];

      Qden_n = gfs0[Qden_np1_gfn - 1];
      Qden_np1 = gfs0[Qden_n_gfn - 1];

      cmask = gfs0[cmask_gfn - 1];

      // set region toggles for fractional integration
      base = 1;
      half = 0;

      // compute desired quantities for each processor on the base level
      proctotal_Qbase =
          qtotcalc_(cmask, Qden_np1, Qden_n, &lshape[0], &lshape[1], lx0[0],
                    lx0[1], &ldR, &ldX, lghost_width, &base, &half, &c, &d);

      proctotal_Ebase =
          etotcalc_(cmask, Eden_np1, Eden_n, &lshape[0], &lshape[1], lx0[0],
                    lx0[1], &ldR, &ldX, lghost_width, &base, &c, &d);

      proctotal_Eembase =
          etotcalc_(cmask, Eem_np1, Eem_n, &lshape[0], &lshape[1], lx0[0],
                    lx0[1], &ldR, &ldX, lghost_width, &base, &c, &d);

      // compute partial L2-norm of desired quantities for each processor on the base level
      // (this ignores ghost points, so will overcount)
      proctotal_phi1IRE =
          gfl2norm_(phi1ire_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_phi2IRE =
          gfl2norm_(phi2ire_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_atIRE =
          gfl2norm_(atire_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_apIRE =
          gfl2norm_(apire_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_azIRE =
          gfl2norm_(azire_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_gaugecond =
          gfl2norm_(gaugecond_np1, &lshape[0], &lshape[1], lghost_width);
      proctotal_divE =
          gfl2norm_(divE_np1, &lshape[0], &lshape[1], lghost_width);

      valid = PAMR_next_g(); // retrieve next grid on the specified level for this processor
    }

    // wait for all processors to reach this point
    MPI_Barrier(MPI_COMM_WORLD);

    // fill recv arrays with values from each processor
    MPI_Gather(&proctotal_phi1IRE, 1, MPI_DOUBLE, recv_phi1IRE, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_phi2IRE, 1, MPI_DOUBLE, recv_phi2IRE, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_atIRE, 1, MPI_DOUBLE, recv_atIRE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_apIRE, 1, MPI_DOUBLE, recv_apIRE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_azIRE, 1, MPI_DOUBLE, recv_azIRE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Qbase, 1, MPI_DOUBLE, recv_Qbase, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Ebase, 1, MPI_DOUBLE, recv_Ebase, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotal_Eembase, 1, MPI_DOUBLE, recv_Eembase, 1, MPI_DOUBLE,
               0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_gaugecond, 1, MPI_DOUBLE, recv_gaugecond, 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&proctotal_divE, 1, MPI_DOUBLE, recv_divE, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);

    // add up totals
    if (my_rank == 0) {
      total_phi1IRE = proctotal_phi1IRE;
      total_phi2IRE = proctotal_phi2IRE;
      total_atIRE = proctotal_atIRE;
      total_apIRE = proctotal_apIRE;
      total_azIRE = proctotal_azIRE;
      total_Qbase = proctotal_Qbase;
      total_Ebase = proctotal_Ebase;
      total_Eembase = proctotal_Eembase;
      total_gaugecond = proctotal_gaugecond;
      total_divE = proctotal_divE;

      for (i = 1; i < num_procs; i++) {
        total_phi1IRE = total_phi1IRE + recv_phi1IRE[num_procs - i];
        total_phi2IRE = total_phi2IRE + recv_phi2IRE[num_procs - i];
        total_atIRE = total_atIRE + recv_atIRE[num_procs - i];
        total_apIRE = total_apIRE + recv_apIRE[num_procs - i];
        total_azIRE = total_azIRE + recv_azIRE[num_procs - i];
        total_Qbase = total_Qbase + recv_Qbase[num_procs - i];
        total_Ebase = total_Ebase + recv_Ebase[num_procs - i];
        total_Eembase = total_Eembase + recv_Eembase[num_procs - i];
        total_gaugecond = total_gaugecond + recv_gaugecond[num_procs - i];
        total_divE = total_divE + recv_divE[num_procs - i];
      }

      // write data to file
      file_phi1IRE = fopen("./IRE_phi1.dat", "a");
      fprintf(file_phi1IRE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_phi1IRE / totgs));
      fclose(file_phi1IRE);

      file_phi2IRE = fopen("./IRE_phi2.dat", "a");
      fprintf(file_phi2IRE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_phi2IRE / totgs));
      fclose(file_phi2IRE);

      file_atIRE = fopen("./IRE_At.dat", "a");
      fprintf(file_atIRE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_atIRE / totgs));
      fclose(file_atIRE);

      file_apIRE = fopen("./IRE_Ap.dat", "a");
      fprintf(file_apIRE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_apIRE / totgs));
      fclose(file_apIRE);

      file_azIRE = fopen("./IRE_Az.dat", "a");
      fprintf(file_azIRE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_azIRE / totgs));
      fclose(file_azIRE);

      file_Qtotbase = fopen("./Qtotbase.dat", "a");
      fprintf(file_Qtotbase, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              total_Qbase);
      fclose(file_Qtotbase);

      file_Etotbase = fopen("./Etotbase.dat", "a");
      fprintf(file_Etotbase, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              total_Ebase);
      fclose(file_Etotbase);

      file_Eemtotbase = fopen("./Eemtotbase.dat", "a");
      fprintf(file_Eemtotbase, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              total_Eembase);
      fclose(file_Eemtotbase);

      file_gaugecond = fopen("./l2norm_gaugecond.dat", "a");
      fprintf(file_gaugecond, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_gaugecond / totgs));
      fclose(file_gaugecond);

      file_divE = fopen("./l2norm_divE.dat", "a");
      fprintf(file_divE, "\n%20.15f %20.15f\n", PAMR_get_time(L),
              sqrt(total_divE / totgs));
      fclose(file_divE);
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
        ldR = ldx0[0];
        ldX = ldx0[1];

        PAMR_get_g_x(lx0); // retrieve pointers to arrays containing coordinates

        // retrieve pointers to arrays containing grid function data
        PAMR_get_g_gfs(gfs0);

        Eden_n = gfs0[Eden_n_gfn - 1];
        Eden_np1 = gfs0[Eden_np1_gfn - 1];

        Eem_n = gfs0[Eem_n_gfn - 1];
        Eem_np1 = gfs0[Eem_np1_gfn - 1];

        Qden_n = gfs0[Qden_n_gfn - 1];
        Qden_np1 = gfs0[Qden_np1_gfn - 1];

        cmask = gfs0[cmask_gfn - 1];

        // set region toggles for fractional integration
        base = 0;
        half = 0;

        // compute desired quantities for each processor
        ltotalQ =
            qtotcalc_(cmask, Qden_np1, Qden_n, &lshape[0], &lshape[1], lx0[0],
                      lx0[1], &ldR, &ldX, lghost_width, &base, &half, &c, &d);

        ltotalE =
            etotcalc_(cmask, Eden_np1, Eden_n, &lshape[0], &lshape[1], lx0[0],
                      lx0[1], &ldR, &ldX, lghost_width, &base, &c, &d);

        ltotalEem =
            etotcalc_(cmask, Eem_np1, Eem_n, &lshape[0], &lshape[1], lx0[0],
                      lx0[1], &ldR, &ldX, lghost_width, &base, &c, &d);

        // integrate in upper half-plane
        half = 1;
        ltotalQhalf_upper =
            qtotcalc_(cmask, Qden_np1, Qden_n, &lshape[0], &lshape[1], lx0[0],
                      lx0[1], &ldR, &ldX, lghost_width, &base, &half, &c, &d);

        // integrate in lower half-plane
        half = -1;
        ltotalQhalf_lower =
            qtotcalc_(cmask, Qden_np1, Qden_n, &lshape[0], &lshape[1], lx0[0],
                      lx0[1], &ldR, &ldX, lghost_width, &base, &half, &c, &d);

        // add contributions from this level to processor total
        proctotalQ += ltotalQ;
        proctotalE += ltotalE;
        proctotalEem += ltotalEem;
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
    MPI_Gather(&proctotalEem, 1, MPI_DOUBLE, recvEem, 1, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Gather(&proctotalQhalf_upper, 1, MPI_DOUBLE, recvQhalf_upper, 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&proctotalQhalf_lower, 1, MPI_DOUBLE, recvQhalf_lower, 1,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // add up totals
    if (my_rank == 0) {
      grandtotalQ += proctotalQ;
      grandtotalE += proctotalE;
      grandtotalEem += proctotalEem;
      grandtotalQhalf_upper += proctotalQhalf_upper;
      grandtotalQhalf_lower += proctotalQhalf_lower;

      for (i = 1; i < num_procs; i++) {
        grandtotalQ = grandtotalQ + recvQ[num_procs - i];
        grandtotalE = grandtotalE + recvE[num_procs - i];
        grandtotalEem = grandtotalEem + recvEem[num_procs - i];
        grandtotalQhalf_upper =
            grandtotalQhalf_upper + recvQhalf_upper[num_procs - i];
        grandtotalQhalf_lower =
            grandtotalQhalf_lower + recvQhalf_lower[num_procs - i];
      }

      // write data to file
      fileQtot = fopen("./Qtot.dat", "a");
      fprintf(fileQtot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalQ);
      fclose(fileQtot);

      fileEtot = fopen("./Etot.dat", "a");
      fprintf(fileEtot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalE);
      fclose(fileEtot);

      fileEemtot = fopen("./Eemtot.dat", "a");
      fprintf(fileEemtot, "%20.15f  %20.15f\n", PAMR_get_time(L), grandtotalEem);
      fclose(fileEemtot);

      fileQtothalf_upper = fopen("./Qtothalf_upper.dat", "a");
      fprintf(fileQtothalf_upper, "%20.15f  %20.15f\n", PAMR_get_time(L),
              grandtotalQhalf_upper);
      fclose(fileQtothalf_upper);

      fileQtothalf_lower = fopen("./Qtothalf_lower.dat", "a");
      fprintf(fileQtothalf_lower, "%20.15f  %20.15f\n", PAMR_get_time(L),
              grandtotalQhalf_lower);
      fclose(fileQtothalf_lower);
    }
  }

  // print elapsed time
  if (my_rank == 0 && L == 1)
    elapsed_time();

  if (my_rank == 0) {
    free(recvQ);
    free(recvE);
    free(recvEem);
    free(recvQhalf_upper);
    free(recvQhalf_lower);
    free(recv_Qbase);
    free(recv_Ebase);
    free(recv_Eembase);
    free(recv_phi1IRE);
    free(recv_phi2IRE);
    free(recv_atIRE);
    free(recv_apIRE);
    free(recv_azIRE);
    free(recv_gaugecond);
    free(recv_divE);
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
  residual_(V_res, V_rhs, V, mask_mg, R, X, &norm, &NR, &NX, &dR, &dX, pR, zX,
            dRdp, dRdp2, dXdz, dXdz2, &e, phi1, phi2, pi1, pi2, modphi, A_t);

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
  relax_(V, V_rhs, mask_mg, phys_bdy, &norm, R, X, &NR, &NX, &dR, &dX, pR, zX,
         dRdp, dRdp2, dXdz, dXdz2, &e, phi1, phi2, pi1, pi2, modphi, A_t);

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
  lop_(V_lop, V, mask_mg, R, X, &NR, &NX, &dR, &dX, pR, zX, dRdp, dRdp2, dXdz,
       dXdz2, &e, phi1, phi2, pi1, pi2, modphi, A_t);
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
       &qball_evo_residual, &qball_MG_residual, &qball_evolve,
       &qball_MG_relax, &qball_L_op, &qball_pre_io_calc, &qball_scale_tre,
       &qball_post_regrid, &qball_post_tstep, &qball_fill_ex_mask,
       &qball_fill_bh_bboxes);

  if (my_rank == 0)
    elapsed_time();

  return 0;
}
