#ifndef QBALLPAMR_H
#define QBALLPAMR_H

/*============================================================================= */
/* Global variables and prototypes                                              */
/*============================================================================= */

extern real cx,cy,cx1,cy1,cx2,cy2,B;
extern real *phi1_n,*phi1_np1;
extern real *phi2_n,*phi2_np1;
extern real *pi1_n,*pi1_np1;
extern real *pi2_n,*pi2_np1;
extern real *modphi_n,*modphi_np1;
extern real *qdensity_n,*qdensity_np1;
extern real *x,*y,*z;

extern int phi1_n_gfn,phi1_np1_gfn;
extern int phi2_n_gfn,phi2_np1_gfn;
extern int pi1_n_gfn,pi1_np1_gfn;
extern int pi2_n_gfn,pi2_np1_gfn;
extern int modphi_n_gfn,modphi_np1_gfn;
extern int qdensity_n_gfn,qdensity_np1_gfn;

extern real *phi1_res;
extern real *phi2_res;
extern real *pi1_res;
extern real *pi2_res;

extern int Nx,Ny;
extern int g_L;
extern real dx,dy,dt;

// Declarations for various PAMR-required variables

extern int shape[2],ghost_width[4],phys_bdy[4];
extern int size,g_rank,dim;
extern real base_bbox[4],bbox[4];

void set_gfns(void);
void ldptr(void);
void const_f(real *f, real c);
void zero(real *f);
real l2norm_calc(real *f);

extern int num_procs;
extern real *recv;

extern real *g_norms;

// Prototypes for the Fortran subroutines we will use, which are defined in the
// corresponding Fortran files. Note that most Fortran compilers append a
// trailing underscore to Fortran routine names, so if we wish to call it from
// a C program then we should manually append the underscore
void init_qball_(real *phi1_np1, real *phi1_n, real *phi2_np1, real *phi2_n,
                 real *pi1_np1, real *pi1_n, real *pi2_np1, real *pi2_n,
                 int *g1_Nx, int *g1_Ny, real *x, real *y, real *dx, real *cx,
                 real *cx1, real *cx2, real *cy, real *cy1, real *cy2,
            		 real *xmax, real *xmin, real *ymax, real *ymin);

void initializer1_(real *modphi_n, real *qdensity_n, real *phi1_n, real *phi2_n,
                   real *pi1_n, real *pi2_n, int *g1_Nx, int *g1_Ny);

// Prototypes for updates (as defined in the RNPL-generated file updates.f)
void update0_(real *modphi_np1, real *modphi_n, real *qdensity_np1, real *phi1_np1,
              real *phi1_n, real *phi2_np1, real *phi2_n, real *pi1_np1, real *pi1_n,
	            real *pi2_np1, real *pi2_n,
	            int *g1_Nx, int *g1_Ny, real *dt, real *dx, real *dy, real *B,
              int *xmax_pb, int *xmin_pb, int *ymax_pb, int *ymin_pb);

double qtotcalcpamr_(real *qdensity_np1, real *qdensity_n, int *g1_Nx,
                     int *g1_Ny, real *dx, real *dy, int *ghost_width);

// Prototypes for residuals (as defined in the RNPL-generated file residuals.f)
void res_phi1_(real *phi1_res,
               real *phi1_np1, real *phi1_n,
               real *pi1_np1,  real *pi1_n,
               int *g1_Nx, int *g1_Ny,
               real *dt, real *dx, real *dy,
               int *xmax_pb, int *xmin_pb, int *ymax_pb, int *ymin_pb);

void res_phi2_(real *phi2_res,
               real *phi2_np1, real *phi2_n,
               real *pi2_np1,  real *pi2_n,
               int *g1_Nx, int *g1_Ny,
               real *dt, real *dx, real *dy,
               int *xmax_pb, int *xmin_pb, int *ymax_pb, int *ymin_pb);

void res_pi1_(real *pi1_res,
              real *modphi_np1, real *modphi_n,
              real *phi1_np1,   real *phi1_n,
              real *pi1_np1,    real *pi1_n,
              int *g1_Nx, int *g1_Ny,
              real *dt, real *dx, real *dy, real *B,
              int *xmax_pb, int *xmin_pb, int *ymax_pb, int *ymin_pb);

void res_pi2_(real *pi2_res,
              real *modphi_np1, real *modphi_n,
              real *phi2_np1,   real *phi2_n,
              real *pi2_np1,    real *pi2_n,
              int *g1_Nx, int *g1_Ny,
              real *dt, real *dx, real *dy, real *B,
              int *xmax_pb, int *xmin_pb, int *ymax_pb, int *ymin_pb);

#endif
