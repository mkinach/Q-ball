#ifndef QBALL_H
#define QBALL_H

void set_cmask_child(int L, int hier);  // also in mg.h

void init_qball_(real *n_Xx, real *n_dxdX, real *n_dxdX2, real *n_Yy, real *n_dydY, real *n_dydY2, real *n_Zz, real *n_dzdZ, real *n_dzdZ2, real *np1_phi1, real *n_phi1, real *np1_phi2, real *n_phi2, real *np1_pi1, real *n_pi1, real *np1_pi2, real *n_pi2, real *np1_at, real *n_at, real *np1_ax, real *n_ax, real *np1_ay, real *n_ay, real *np1_az, real *n_az, real *np1_bt, real *n_bt, real *np1_bx, real *n_bx, real *np1_by, real *n_by, real *np1_bz, real *n_bz, int *Nx, int *Ny, int *Nz, real *x, real *y, real *z, real *dx, real *dy, real *dz, real *xmax, real *xmin, real *ymax, real *ymin, real *zmax, real *zmin, real *cx, real *cy, real *cz, real *cx1, real *cx2, real *cy1, real *cy2, real *cz1, real *cz2, real *w, real *v, real *phase, int *anti, int *ini_type, real *e);

void init_xx_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dxdx2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dxdx_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dxdx2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_yy_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dydy2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dydy_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dydy2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_zz_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dzdz2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dzdz_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_dzdz2_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *myzero, double *res);

void init_jac_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *res);

void init_chi_(double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *amp, double *c, double *chiX0, double *chiY0, double *chiZ0, double *d, double *delta, double *r0, double *wwX, double *wwY, double *wwZ, double *res);

void init_xi_(double *n_chi, double *n_dxdX, double *n_dydY, double *n_dzdZ, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *hx, double *hy, double *hz, double *res);

void init_modphi_(double *n_phi1, double *n_phi2, int *Nx, int *Ny, int *Nz, double *res);

void init_qden_(double *n_At, double *n_modphi, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, int *Nx, int *Ny, int *Nz, double *e, double *res);

void init_eden_log_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_Bx, double *n_By, double *n_Bz, double *n_dxdX, double *n_dydY, double *n_dzdZ, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, int *Nx, int *Ny, int *Nz, double *beta, double *e, double *hx, double *hy, double *hz, double *mu, double *res);

void init_eden_poly_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_Bx, double *n_By, double *n_Bz, double *n_dxdX, double *n_dydY, double *n_dzdZ, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, int *Nx, int *Ny, int *Nz, double *e, double *g, double *h, double *hx, double *hy, double *hz, double *m, double *res);

void init_jdenx_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_Bx, double *n_By, double *n_Bz, double *n_dxdX, double *n_dydY, double *n_dzdZ, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *e, double *hx, double *hy, double *hz, double *res);

void init_lorenz_(double *n_Ax, double *n_Ay, double *n_Az, double *n_Bt, double *n_dxdX, double *n_dydY, double *n_dzdZ, int *Nx, int *Ny, int *Nz, double *hx, double *hy, double *hz, double *res);

void init_elecx_(double *n_At, double *n_Bx, double *n_dxdX, int *Nx, int *Ny, int *Nz, double *hx, double *res);

void init_elecy_(double *n_At, double *n_By, double *n_dydY, int *Nx, int *Ny, int *Nz, double *hy, double *res);

void init_elecz_(double *n_At, double *n_Bz, double *n_dzdZ, int *Nx, int *Ny, int *Nz, double *hz, double *res);

void init_magx_(double *n_Ay, double *n_Az, double *n_dydY, double *n_dzdZ, int *Nx, int *Ny, int *Nz, double *hy, double *hz, double *res);

void init_magy_(double *n_Ax, double *n_Az, double *n_dxdX, double *n_dzdZ, int *Nx, int *Ny, int *Nz, double *hx, double *hz, double *res);

void init_magz_(double *n_Ax, double *n_Ay, double *n_dxdX, double *n_dydY, int *Nx, int *Ny, int *Nz, double *hx, double *hy, double *res);

void init_eem_(double *n_elecx, double *n_elecy, double *n_elecz, double *n_magx, double *n_magy, double *n_magz, int *Nx, int *Ny, int *Nz, double *res);

void init_gaussb_(double *n_dxdX, double *n_dydY, double *n_dzdZ, double *n_magx, double *n_magy, double *n_magz, int *Nx, int *Ny, int *Nz, double *hx, double *hy, double *hz, double *res);

void init_gausse_(double *n_At, double *n_Bx, double *n_By, double *n_Bz, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_qden, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, double *res);

void u_phi1_rk_(double *n_pi1, double *np1_phi1, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_phi2_rk_(double *n_pi2, double *np1_phi2, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_pi1_poly_rk_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi2, double *np1_pi1, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *cc, double *e, double *g, double *h, double *hx, double *hy, double *hz, double *m, double *myzero, int *phys_bdy, double *res);

void u_pi2_poly_rk_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *np1_pi2, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *cc, double *e, double *g, double *h, double *hx, double *hy, double *hz, double *m, double *myzero, int *phys_bdy, double *res);

void u_pi1_log_rk_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi2, double *np1_pi1, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *beta, double *cc, double *e, double *hx, double *hy, double *hz, double *mu, double *myzero, int *phys_bdy, double *res);

void u_pi2_log_rk_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *np1_pi2, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *beta, double *cc, double *e, double *hx, double *hy, double *hz, double *mu, double *myzero, int *phys_bdy, double *res);

void u_at_rk_(double *n_Bt, double *np1_At, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_ax_rk_(double *n_Bx, double *np1_Ax, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_ay_rk_(double *n_By, double *np1_Ay, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_az_rk_(double *n_Bz, double *np1_Az, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_bt_rk_(double *n_At, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, double *np1_Bt, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, double *myzero, int *phys_bdy, double *res);

void u_bx_rk_(double *n_Ax, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *np1_Bx, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, double *myzero, int *phys_bdy, double *res);

void u_by_rk_(double *n_Ay, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *np1_By, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, double *myzero, int *phys_bdy, double *res);

void u_bz_rk_(double *n_Az, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *np1_Bz, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, double *myzero, int *phys_bdy, double *res);

void u_chi_rk_(double *n_xi, double *np1_chi, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *myzero, int *phys_bdy, double *res);

void u_xi_rk_(double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *np1_xi, double *x, double *y, double *z, int *Nx, int *Ny, int *Nz, double *cc, double *hx, double *hy, double *hz, double *myzero, int *phys_bdy, double *res);

void u_modphi_(double *np1_phi1, double *np1_phi2, int *Nx, int *Ny, int *Nz, int *phys_bdy, double *res);

void u_qden_(double *np1_At, double *np1_modphi, double *np1_phi1, double *np1_phi2, double *np1_pi1, double *np1_pi2, int *Nx, int *Ny, int *Nz, double *e, int *phys_bdy, double *res);

void u_eden_log_(double *np1_At, double *np1_Ax, double *np1_Ay, double *np1_Az, double *np1_Bx, double *np1_By, double *np1_Bz, double *np1_dxdX, double *np1_dydY, double *np1_dzdZ, double *np1_phi1, double *np1_phi2, double *np1_pi1, double *np1_pi2, int *Nx, int *Ny, int *Nz, double *beta, double *e, double *hx, double *hy, double *hz, double *mu, int *phys_bdy, double *res);

void u_eden_poly_(double *np1_At, double *np1_Ax, double *np1_Ay, double *np1_Az, double *np1_Bx, double *np1_By, double *np1_Bz, double *np1_dxdX, double *np1_dydY, double *np1_dzdZ, double *np1_phi1, double *np1_phi2, double *np1_pi1, double *np1_pi2, int *Nx, int *Ny, int *Nz, double *e, double *g, double *h, double *hx, double *hy, double *hz, double *m, int *phys_bdy, double *res);

void u_jdenx_(double *np1_At, double *np1_Ax, double *np1_Ay, double *np1_Az, double *np1_Bx, double *np1_By, double *np1_Bz, double *np1_dxdX, double *np1_dydY, double *np1_dzdZ, double *np1_phi1, double *np1_phi2, double *np1_pi1, double *np1_pi2, double *y, double *z, int *Nx, int *Ny, int *Nz, double *c, double *d, double *e, double *hx, double *hy, double *hz, int *phys_bdy, double *res);

void u_lorenz_(double *np1_Ax, double *np1_Ay, double *np1_Az, double *np1_Bt, double *np1_dxdX, double *np1_dydY, double *np1_dzdZ, int *Nx, int *Ny, int *Nz, double *hx, double *hy, double *hz, int *phys_bdy, double *res);

void u_elecx_(double *np1_At, double *np1_Bx, double *np1_dxdX, int *Nx, int *Ny, int *Nz, double *hx, int *phys_bdy, double *res);

void u_elecy_(double *np1_At, double *np1_By, double *np1_dydY, int *Nx, int *Ny, int *Nz, double *hy, int *phys_bdy, double *res);

void u_elecz_(double *np1_At, double *np1_Bz, double *np1_dzdZ, int *Nx, int *Ny, int *Nz, double *hz, int *phys_bdy, double *res);

void u_magx_(double *np1_Ay, double *np1_Az, double *np1_dydY, double *np1_dzdZ, int *Nx, int *Ny, int *Nz, double *hy, double *hz, int *phys_bdy, double *res);

void u_magy_(double *np1_Ax, double *np1_Az, double *np1_dxdX, double *np1_dzdZ, int *Nx, int *Ny, int *Nz, double *hx, double *hz, int *phys_bdy, double *res);

void u_magz_(double *np1_Ax, double *np1_Ay, double *np1_dxdX, double *np1_dydY, int *Nx, int *Ny, int *Nz, double *hx, double *hy, int *phys_bdy, double *res);

void u_eem_(double *np1_elecx, double *np1_elecy, double *np1_elecz, double *np1_magx, double *np1_magy, double *np1_magz, int *Nx, int *Ny, int *Nz, int *phys_bdy, double *res);

void u_gaussb_(double *np1_dxdX, double *np1_dydY, double *np1_dzdZ, double *np1_magx, double *np1_magy, double *np1_magz, int *Nx, int *Ny, int *Nz, double *hx, double *hy, double *hz, int *phys_bdy, double *res);

void u_gausse_(double *np1_At, double *np1_Bx, double *np1_By, double *np1_Bz, double *np1_dxdX, double *np1_dxdX2, double *np1_dydY, double *np1_dydY2, double *np1_dzdZ, double *np1_dzdZ2, double *np1_qden, int *Nx, int *Ny, int *Nz, double *e, double *hx, double *hy, double *hz, int *phys_bdy, double *res);

void ire_phi1_log_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi2, double *nm1_phi1, double *np1_phi1, int *Nx, int *Ny, int *Nz, double *beta, double *cc, double *e, double *ht, double *hx, double *hy, double *hz, double *mu, double *res);

void ire_phi1_poly_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi2, double *nm1_phi1, double *np1_phi1, int *Nx, int *Ny, int *Nz, double *cc, double *e, double *g, double *h, double *ht, double *hx, double *hy, double *hz, double *m, double *res);

void ire_phi2_log_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *nm1_phi2, double *np1_phi2, int *Nx, int *Ny, int *Nz, double *beta, double *cc, double *e, double *ht, double *hx, double *hy, double *hz, double *mu, double *res);

void ire_phi2_poly_(double *n_At, double *n_Ax, double *n_Ay, double *n_Az, double *n_chi, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *nm1_phi2, double *np1_phi2, int *Nx, int *Ny, int *Nz, double *cc, double *e, double *g, double *h, double *ht, double *hx, double *hy, double *hz, double *m, double *res);

void ire_at_(double *n_At, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *n_pi1, double *n_pi2, double *nm1_At, double *np1_At, int *Nx, int *Ny, int *Nz, double *e, double *ht, double *hx, double *hy, double *hz, double *res);

void ire_ax_(double *n_Ax, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *nm1_Ax, double *np1_Ax, int *Nx, int *Ny, int *Nz, double *e, double *ht, double *hx, double *hy, double *hz, double *res);

void ire_ay_(double *n_Ay, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *nm1_Ay, double *np1_Ay, int *Nx, int *Ny, int *Nz, double *e, double *ht, double *hx, double *hy, double *hz, double *res);

void ire_az_(double *n_Az, double *n_dxdX, double *n_dxdX2, double *n_dydY, double *n_dydY2, double *n_dzdZ, double *n_dzdZ2, double *n_phi1, double *n_phi2, double *nm1_Az, double *np1_Az, int *Nx, int *Ny, int *Nz, double *e, double *ht, double *hx, double *hy, double *hz, double *res);

double qtotcalc_(real *cmask, real *Qden_np1, real *Qden_n, real *jac_n, int *Nx, int *Ny, int *Nz, real *x, real *y, real *z, real *dx, real *dy, real *dz, int *gw, int *base, int *half, real *c, real *d);

double etotcalc_(real *cmask, real *Eden_np1, real *Eden_n, real *jac_n, int *Nx, int *Ny, int *Nz, real *x, real *y, real *z, real *dx, real *dy, real *dz, int *gw, int *base, real *c, real *d);

void initguess_(real *V, real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, real *Xx, real *Yy, real *Zz, real *At);

void lop_(real *LV, real *V, real *defect, real *cmask, real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, real *dx, real *dy, real *dz, real *Xx, real *Yy, real *Zz, real *dxdX, real *dxdX2, real *dydY, real *dydY2, real *dzdZ, real *dzdZ2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *At);

void residual_(real *res, real *rhs, real *V, real *defect, real *cmask, real *x, real *y, real *z, real *norm, int *Nx, int *Ny, int *Nz, real *dx, real *dy, real *dz, real *Xx, real *Yy, real *Zz, real *dxdX, real *dxdX2, real *dydY, real *dydY2, real *dzdZ, real *dzdZ2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *At);

void relax_(real *V, real *V_rhs, real *defect, real *cmask, int *phys_bdy, real *norm, real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, real *dx, real *dy, real *dz, real *Xx, real *Yy, real *Zz, real *dxdX, real *dxdX2, real *dydY, real *dydY2, real *dzdZ, real *dzdZ2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *At);

void dc_(real *V, real *defect, real *x, real *y, real *z, real *dx, real *dy, real *dz, int *Nx, int *Ny, int *Nz, real *Xx, real *Yy, real *Zz, real *dxdX, real *dxdX2, real *dydY, real *dydY2, real *dzdZ, real *dzdZ2, real *e, real *phi1, real *phi2, real *pi1, real *pi2, real *modphi, real *At);

void init_qball2_(real *V_n, real *x, real *y, real *z, int *Nx, int *Ny, int *Nz, real *dx, real *dy, real *dz, real *Xx_n, real *Yy_n, real *Zz_n, real *dxdX_n, real *dxdX2_n, real *dydY_n, real *dydY2_n, real *dzdZ_n, real *dzdZ2_n, real *At_n, real *Ax_n, real *Ay_n, real *Az_n, real *Bt_n, real *Bx_n, real *By_n, real *Bz_n, real *xmin, real *xmax, real *ymin, real *ymax, real *zmin, real *zmax);

void w2f_(real *V, int *Nx, int *Ny, int *Nz);

void rff_(real *V, int *Nx, int *Ny, int *Nz);

#endif
