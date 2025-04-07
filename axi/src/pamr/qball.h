#ifndef QBALL_H
#define QBALL_H

void init_qball_(real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *B_t_np1, real *B_t_n, real *Btilde_p_np1, real *Btilde_p_n, real *Btilde_phi_np1, real *Btilde_phi_n, real *B_z_np1, real *B_z_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *R, real *X, real *dR, real *dX, real *c, real *cR, real *cR1, real *cR2, real *cX, real *cX1, real *cX2, real *d, real *e, real *Rmax, real *Rmin, real *w, real *Xmax, real *Xmin);

void initializer0_(real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *R, real *X, real *c, real *d);

void initializer1_(real *phi1ire_n, real *phi2ire_n, real *atire_n, real *apire_n, real *azire_n, int *g1_NR, int *g1_NX);

void initializer3_(real *modphi_n, real *phi1_n, real *phi2_n, real *pi1_n, real *pi2_n, real *A_t_n, real *Eden_n, real *Eem_n, real *Qden_n, real *elec_p_n, real *elec_phi_n, real *elec_z_n, real *mag_p_n, real *mag_phi_n, real *mag_z_n, real *gaugecond_n, real *divE_n, real *divB_n, real *chi_n, real *xi_n, real *dummy_n, real *dRdp_n, real *dXdz_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dX, real *amp, real *delta, real *e, real *r0, real *R_center, real *Rw, real *signum, real *X_center, real *Xw);

void res_chi_(real *chi_res, real *chi_np1, real *chi_n, real *xi_np1, real *xi_n, int *g1_NR, int *g1_NX, real *dt);

void res_xi_(real *xi_res, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *chi_np1, real *chi_n, real *xi_np1, real *xi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *cc, int *Rmax_pb, int *Rmin_pb, int *Xmax_pb, int *Xmin_pb);

void res_phi1_(real *phi1_res, real *phi1_np1, real *phi1_n, real *pi1_np1, real *pi1_n, int *g1_NR, int *g1_NX, real *dt);

void res_phi2_(real *phi2_res, real *phi2_np1, real *phi2_n, real *pi2_np1, real *pi2_n, int *g1_NR, int *g1_NX, real *dt);

void res_pi1_(real *pi1_res, real *modphi_np1, real *modphi_n, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *chi_np1, real *chi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *beta, real *cc, real *e, real *gg, real *h, real *mm, real *mu, int *Rmax_pb, int *Rmin_pb, int *Vtype, int *Xmax_pb, int *Xmin_pb);

void res_pi2_(real *pi2_res, real *modphi_np1, real *modphi_n, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *chi_np1, real *chi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *beta, real *cc, real *e, real *gg, real *h, real *mm, real *mu, int *Rmax_pb, int *Rmin_pb, int *Vtype, int *Xmax_pb, int *Xmin_pb);

void res_a_t_(real *A_t_res, real *A_t_np1, real *A_t_n, real *B_t_np1, real *B_t_n, int *g1_NR, int *g1_NX, real *dt);

void res_atilde_p_(real *Atilde_p_res, real *Atilde_p_np1, real *Atilde_p_n, real *Btilde_p_np1, real *Btilde_p_n, int *g1_NR, int *g1_NX, real *dt);

void res_atilde_phi_(real *Atilde_phi_res, real *Atilde_phi_np1, real *Atilde_phi_n, real *Btilde_phi_np1, real *Btilde_phi_n, int *g1_NR, int *g1_NX, real *dt);

void res_a_z_(real *A_z_res, real *A_z_np1, real *A_z_n, real *B_z_np1, real *B_z_n, int *g1_NR, int *g1_NX, real *dt);

void res_b_t_(real *B_t_res, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *B_t_np1, real *B_t_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *e, int *Rmax_pb, int *Rmin_pb, int *Xmax_pb, int *Xmin_pb);

void res_btilde_p_(real *Btilde_p_res, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *Atilde_p_np1, real *Atilde_p_n, real *Btilde_p_np1, real *Btilde_p_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *e, int *Rmax_pb, int *Rmin_pb, int *Xmax_pb, int *Xmin_pb);

void res_btilde_phi_(real *Btilde_phi_res, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *Btilde_phi_np1, real *Btilde_phi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *e, int *Rmax_pb, int *Rmin_pb, int *Xmax_pb, int *Xmin_pb);

void res_b_z_(real *B_z_res, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *A_z_np1, real *A_z_n, real *B_z_np1, real *B_z_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *e, int *Rmax_pb, int *Rmin_pb, int *Xmax_pb, int *Xmin_pb);

void update0_(real *modphi_np1, real *modphi_n, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *B_t_np1, real *B_t_n, real *Btilde_p_np1, real *Btilde_p_n, real *Btilde_phi_np1, real *Btilde_phi_n, real *B_z_np1, real *B_z_n, real *Eem_np1, real *Qden_np1, real *elec_p_np1, real *elec_phi_np1, real *elec_z_np1, real *mag_p_np1, real *mag_phi_np1, real *mag_z_np1, real *gaugecond_np1, real *divE_np1, real *divB_np1, real *chi_np1, real *chi_n, real *xi_np1, real *xi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *dR, real *dt, real *dX, real *beta, real *cc, real *e, real *gg, real *h, real *mm, real *mu, int *Rmax_pb, int *Rmin_pb, int *Vtype, int *Xmax_pb, int *Xmin_pb);

void ire_(real *modphi_np1, real *modphi_n, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *B_t_np1, real *B_t_n, real *Btilde_p_np1, real *Btilde_p_n, real *Btilde_phi_np1, real *Btilde_phi_n, real *B_z_np1, real *B_z_n, real *chi_np1, real *chi_n, real *xi_np1, real *xi_n, real *phi1ire_np1, real *phi1ire_n, real *phi2ire_np1, real *phi2ire_n, real *atire_np1, real *atire_n, real *apire_np1, real *apire_n, real *azire_np1, real *azire_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *R, real *X, real *dR, real *dt, real *dX, real *beta, real *c, real *cc, real *d, real *e, real *gg, real *h, real *mm, real *mu, int *Vtype);

void edencalc_(real *modphi_np1, real *modphi_n, real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n, real *A_t_np1, real *A_t_n, real *Atilde_p_np1, real *Atilde_p_n, real *Atilde_phi_np1, real *Atilde_phi_n, real *A_z_np1, real *A_z_n, real *B_t_np1, real *B_t_n, real *Btilde_p_np1, real *Btilde_p_n, real *Btilde_phi_np1, real *Btilde_phi_n, real *B_z_np1, real *B_z_n, real *Eden_np1, real *Eden_n, real *chi_np1, real *chi_n, real *xi_np1, real *xi_n, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *pR_n, real *zX_n, int *g1_NR, int *g1_NX, real *R, real *X, real *dR, real *dt, real *dX, real *beta, real *c, real *cc, real *d, real *e, real *gg, real *h, real *mm, real *mu, int *Vtype);

double qtotcalc_(real *cmask, real *Qden_np1, real *Qden_n, int *g1_Np, int *g1_Nz, real *R, real *X, real *dR, real *dX, int *ghost_width, int *base, int *half, real *c, real *d);

double etotcalc_(real *cmask, real *Eden_np1, real *Eden_n, int *g1_NR, int *g1_NX, real *R, real *X, real *dR, real *dX, int *ghost_width, int *base, real *c, real *d);

double gfl2norm_(real *gaugecond_np1, int *g1_NR, int *g1_NX, int *ghost_width);

void initguess_(real *f, real *R, real *X, int *NR, int *NX, real *pR, real *zX, real *A_t);

void lop_(real *LV, real *V, real *cmask, real *R, real *X, int *NR, int *NX, real *dR, real *dX, real *pR, real *zX, real *dRdp, real *dRdp2, real *dXdz, real *dXdz2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *A_t);

void residual_(real *res, real *rhs, real *V, real *cmask, real *R, real *X, real *norm, int *NR, int *NX, real *dR, real *dX, real *pR, real *zX, real *dRdp, real *dRdp2, real *dXdz, real *dXdz2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *A_t);

void relax_(real *V, real *V_rhs, real *cmask, int *phys_bdy, real *norm, real *R, real *X, int *NR, int *NX, real *dR, real *dX, real *pR, real *zX, real *dRdp, real *dRdp2, real *dXdz, real *dXdz2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *A_t);

void init_qball2_(real *V_n, real *R, real *X, int *NR, int *NX, real *pR_n, real *zX_n, real *dR, real *dX, real *dRdp_n, real *dRdp2_n, real *dXdz_n, real *dXdz2_n, real *A_t_n, real *Atilde_p_n, real *A_z_n, real *B_t_n, real *Btilde_p_n, real *B_z_n, real *Rmin, real *Rmax, real *Xmin, real *Xmax);

void w2f_(real *V, int *NR, int *NX);

void rff_(real *V, int *NR, int *NX);

#endif
