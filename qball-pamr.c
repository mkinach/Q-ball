//=============================================================================
// APPLICATION INTERFACE FUNCTIONS 
//
// For more information about built-in PAMR/AMRD functionality, look in the
// PAMR Reference Manual, AMRD Reference Manual, pamr.h, amrd.h, or pamr.c
//
//=============================================================================

#include <stdlib.h>
#include <stdio.h>
#include <pamr.h>
#include <amrd.h>
#include <math.h>
#include <mpi.h>

#include "qball-pamr.h"  // custom header for this project

#include <sys/timeb.h>
void elapsed_time(void);

//=============================================================================
// Model parameters
//=============================================================================

real cx,cy,cx1,cy1,cx2,cy2,B;

//=============================================================================
// Some convenient "local" global variables
//=============================================================================

// pointers to store the memory address of grid function values
real *phi1_n,*phi1_np1;
real *phi2_n,*phi2_np1;
real *pi1_n,*pi1_np1;
real *pi2_n,*pi2_np1;
real *modphi_n,*modphi_np1;
real *qdensity_n,*qdensity_np1;

// pointers to store the memory address of grid functions holding residuals
real *phi1_res;
real *phi2_res;
real *pi1_res;
real *pi2_res;

// pointers to store the memory address of coordinates
real *x,*y,*z;

// variable & array declarations
int shape[2];  // number of grid pts in each dim (shape[0]=Nx, shape[1]=Ny)
int ghost_width[4];  // size of ghost cells (a 2*dim-sized array)
int phys_bdy[4]; // whether the boundary is physical or AMR (size = dim*2)
int size; // size=Nx*Ny
int g_rank; // store the MPI rank of the execution
int dim;
real base_bbox[4],bbox[4]; // stores physical and AMR bounding box data

int Nx,Ny;
int g_L;  // stores grid level in AMR hierarchy
real dx,dy,dt;

real *g_norms; // stores L-infinity norm (see below)

int num_procs;
real *recv;

// variables for storing grid function numbers (gfn)
int phi1_n_gfn,phi1_np1_gfn; 
int phi2_n_gfn,phi2_np1_gfn; 
int pi1_n_gfn,pi1_np1_gfn; 
int pi2_n_gfn,pi2_np1_gfn;
int modphi_n_gfn,modphi_np1_gfn;
int qdensity_n_gfn,qdensity_np1_gfn;
int phi1_res_gfn,phi2_res_gfn,pi1_res_gfn,pi2_res_gfn;

//=============================================================================
// Call this function after variables have been defined. It uses PAMR_get_gfn
// to map a particular grid function (phi1, phi2, pi1, ...) to its grid function
// number (GFN) with error checking. Specifically, PAMR_get_gfn returns the GFN
// of gridfunction "phi" in the hierarchy defined by the PAMR_AMRH variable
// (from /include/pamr.h) at time level 2 or 1. A GFN is basically an index
// into an array of pointers that points to the actual grid function data
//=============================================================================

void set_gfns(void)
{
    if ((phi1_n_gfn   = PAMR_get_gfn("phi1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi1_np1_gfn = PAMR_get_gfn("phi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((phi2_n_gfn   = PAMR_get_gfn("phi2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((phi2_np1_gfn = PAMR_get_gfn("phi2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((pi1_n_gfn   = PAMR_get_gfn("pi1",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((pi1_np1_gfn = PAMR_get_gfn("pi1",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((pi2_n_gfn   = PAMR_get_gfn("pi2",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((pi2_np1_gfn = PAMR_get_gfn("pi2",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((modphi_n_gfn   = PAMR_get_gfn("modphi",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((modphi_np1_gfn = PAMR_get_gfn("modphi",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    if ((qdensity_n_gfn   = PAMR_get_gfn("qdensity",PAMR_AMRH,2))<0) AMRD_stop("set_gnfs error",0);
    if ((qdensity_np1_gfn = PAMR_get_gfn("qdensity",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    
    if ((phi1_res_gfn   = PAMR_get_gfn("phi1_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((phi2_res_gfn   = PAMR_get_gfn("phi2_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((pi1_res_gfn   = PAMR_get_gfn("pi1_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);
    if ((pi2_res_gfn   = PAMR_get_gfn("pi2_res",PAMR_AMRH,1))<0) AMRD_stop("set_gnfs error",0);

    g_norms=AMRD_get_global_norms();

}

//=============================================================================
// Load pointers; call with valid iter to set up globals
//=============================================================================
void ldptr(void)
{
   real dx0[2];  // stores dx, dy spacing 
   real *x0[2],*gfs[PAMR_MAX_GFNS];
   
   static int first=1;
   if (first) 
   {
      first=0; 
      set_gfns(); // see function definition above
      
      // initialize the coordinate bounding box of the base grid
      // (e.g. [x1,x2,y1,y2] in 2D)
      PAMR_get_global_bbox(base_bbox);
   }

   PAMR_get_g_dim(&dim);      
   PAMR_get_g_rank(&g_rank);
   PAMR_get_g_shape(shape);
   PAMR_get_g_bbox(bbox);
   PAMR_get_g_ghost_width(ghost_width);
   PAMR_get_g_level(&g_L);
   PAMR_get_dxdt(g_L,dx0,&dt); // fills dx0 and dt
   dx=dx0[0];
   dy=dx0[1];

   // some more PAMR initializations for detecting physical boundaries
   if ((bbox[0]-base_bbox[0])<dx/2) phys_bdy[0]=1; else phys_bdy[0]=0;
   if ((base_bbox[1]-bbox[1])<dx/2) phys_bdy[1]=1; else phys_bdy[1]=0;
   Nx=shape[0];
   size=Nx;
   Ny=1;
   if ((bbox[2]-base_bbox[2])<dy/2) phys_bdy[2]=1; else phys_bdy[2]=0;
   if ((base_bbox[3]-bbox[3])<dy/2) phys_bdy[3]=1; else phys_bdy[3]=0;
   Ny=shape[1];
   size*=Ny;
   
   // fills the pointer array x0 with pointers to arrays containing coordinates of the grid points
   PAMR_get_g_x(x0);

   x=x0[0];
   y=x0[1]; 

   // similar to PAMR_get_g_x except it's for grid function data (not coords)
   PAMR_get_g_gfs(gfs);

   // make phi1_n, phi1_np1, etc. point to the correct grid function data
   // note that gfs[] is a pointer to a real array
   phi1_n   = gfs[phi1_n_gfn-1];
   phi1_np1 = gfs[phi1_np1_gfn-1];

   phi2_n   = gfs[phi2_n_gfn-1];
   phi2_np1 = gfs[phi2_np1_gfn-1];

   pi1_n    = gfs[pi1_n_gfn-1];
   pi1_np1  = gfs[pi1_np1_gfn-1];

   pi2_n    = gfs[pi2_n_gfn-1];
   pi2_np1  = gfs[pi2_np1_gfn-1];

   modphi_n    = gfs[modphi_n_gfn-1];
   modphi_np1  = gfs[modphi_np1_gfn-1];
   
   qdensity_n    = gfs[qdensity_n_gfn-1];
   qdensity_np1  = gfs[qdensity_np1_gfn-1];

   phi1_res = gfs[phi1_res_gfn-1];
   phi2_res = gfs[phi2_res_gfn-1];
   pi1_res = gfs[pi1_res_gfn-1];
   pi2_res = gfs[pi2_res_gfn-1];
}

//=============================================================================
// Utility routine to define a constant function or zero function
//=============================================================================
void const_f(real *f, real c)
{
   int i;

   for (i=0; i<Nx*Ny; i++) f[i]=c;
}

void zero(real *f)
{
   const_f(f,0);
}

real l2norm_calc(real *f)
{
  int i;
  real norm=0;
  int sum=0;

  for (i=0; i<size; i++)
  {
    sum++;
    norm+=f[i]*f[i];
  }
  
  // if sum=0 somehow, set sum=1 to avoid NaN
  if (!sum) sum=1;

  // return the L2-norm
  return (sqrt(norm/sum));
}

//#############################################################################
// Routines required by AMRD 
//#############################################################################

//=============================================================================
// Calculates the elapsed time for the calculation. If this returns 0, the 
// default mechanism is used for initial hierarchy. If this returns 1, then
// this function is expected to calculate the correct initial hierarchy. Here
// I take advantage of it to compute the elapsed time of the calculation.
//=============================================================================
int qball_id(void)
{
	if( my_rank == 0 ) elapsed_time();
   return 0;
}

//=============================================================================
// Sets custom parameters, variables, etc. Split up into two segments,
// one called before the PAMR built-in parameters are read from the *.param
// file and the base hierarchy has been initialized, and the other after
//=============================================================================

void qball_var_pre_init(char *pfile)
{
   return;
}

void qball_var_post_init(char *pfile)
{
   int i,j;
   char buf[64];

   if (my_rank==0)
   {
      printf("===================================================================\n");
      printf("Reading parameters:\n\n");
   }
  
   // initialize variables before reading from .param file
   cx=cy=cx1=cy1=cx2=cy2=B=0;

   // read the custom parameters from the .param file, which should be written in RNPL
   // format. Note this function expects pointers:
   // void AMRD_real_param(char *pfile, char *name, real *var, int size);
   AMRD_real_param(pfile,"cx",&cx,1);
   AMRD_real_param(pfile,"cy",&cy,1);
   AMRD_real_param(pfile,"cx1",&cx1,1);
   AMRD_real_param(pfile,"cy1",&cy1,1);
   AMRD_real_param(pfile,"cx2",&cx2,1);
   AMRD_real_param(pfile,"cy2",&cy2,1);
   AMRD_real_param(pfile,"B",&B,1);

   if (my_rank==0) printf("===================================================================\n");
   return;
}

//=============================================================================
// Set all the variables in the AMR hierarchy to their "zero" values
//=============================================================================
void qball_AMRH_var_clear(void)
{
   ldptr();

   zero(phi1_n); zero(phi1_np1); 
   zero(phi2_n); zero(phi2_np1); 
   zero(pi1_n);  zero(pi1_np1); 
   zero(pi2_n);  zero(pi2_np1); 
   zero(modphi_n);   zero(modphi_np1); 
   zero(qdensity_n); zero(qdensity_np1); 

   zero(phi1_res); zero(phi2_res); zero(pi1_res); zero(pi2_res); 

   return;
}

//=============================================================================
// Initial data generator
//=============================================================================
void qball_free_data(void)
{
   ldptr();

   // call the relevant Fortran subroutine (defined in initializers.f and 
   // prototyped in qball-pamr.h) by passing *memory addresses* for scalar
   // quantities, since Fortran is call-by-address but C is call-by-reference
   init_qball_(phi1_np1,phi1_n,phi2_np1,phi2_n,pi1_np1,pi1_n,pi2_np1,pi2_n,
               &Nx,&Ny,x,y,&dx,&cx,&cx1,&cx2,&cy,&cy1,&cy2,&base_bbox[1],
               &base_bbox[0],&base_bbox[3],&base_bbox[2]);

   initializer1_(modphi_n,qdensity_n,phi1_n,phi2_n,pi1_n,pi2_n,&Nx,&Ny);

   return;
}  

//=============================================================================
// Initial constraint data: called after each MG iteration if you use multigrid
//=============================================================================
void qball_t0_cnst_data(void)
{
   return;
}

//=============================================================================
// Calculations prior to saving info to disk. This allows you to calculate 
// diagnostic grid functions outside of the evolution loop.
// NOTE: at this point, the time sequence is: n,nm1,np1
//=============================================================================
void qball_pre_io_calc(void)
{
   return;
}

//=============================================================================
// Returns some norm of the residual for the evolution variables,
// called after an evolution iteration.
//=============================================================================
real qball_evo_residual(void)
{
  real l2norm=0;
  real l2norm_phi1;
  real l2norm_phi2;
  real l2norm_pi1;
  real l2norm_pi2;

  ldptr();

  // calulate residuals using subroutines defined in residuals.f
  res_phi1_(phi1_res,
            phi1_np1,phi1_n,
            pi1_np1,pi1_n,
            &Nx,&Ny,&dt,&dx,&dy,
            &phys_bdy[1],&phys_bdy[0],&phys_bdy[3],&phys_bdy[2]);

  res_phi2_(phi2_res,
            phi2_np1,phi2_n,
            pi2_np1,pi2_n,
            &Nx,&Ny,&dt,&dx,&dy,
            &phys_bdy[1],&phys_bdy[0],&phys_bdy[3],&phys_bdy[2]);

  res_pi1_(pi1_res,
           modphi_np1,modphi_n,
           phi1_np1,phi1_n,
           pi1_np1,pi1_n,
           &Nx,&Ny,&dt,&dx,&dy,&B,
           &phys_bdy[1],&phys_bdy[0],&phys_bdy[3],&phys_bdy[2]);

  res_pi2_(pi2_res,
           modphi_np1,modphi_n,
           phi2_np1,phi2_n,
           pi2_np1,pi2_n,
           &Nx,&Ny,&dt,&dx,&dy,&B,
           &phys_bdy[1],&phys_bdy[0],&phys_bdy[3],&phys_bdy[2]);

  l2norm_phi1=l2norm_calc(phi1_res);
  l2norm_phi2=l2norm_calc(phi2_res);
  l2norm_pi1=l2norm_calc(pi1_res);
  l2norm_pi2=l2norm_calc(pi2_res);

  // Normalize the computed L2-norms
  if ( g_norms[phi1_n_gfn-1] > 0 ) { l2norm+=l2norm_phi1/g_norms[phi1_n_gfn-1]; }
  if ( g_norms[phi2_n_gfn-1] > 0 ) { l2norm+=l2norm_phi2/g_norms[phi2_n_gfn-1]; }
  if ( g_norms[pi1_n_gfn-1] > 0 )  { l2norm+=l2norm_pi1/g_norms[pi1_n_gfn-1];   }   
  if ( g_norms[pi1_n_gfn-1] > 0 )  { l2norm+=l2norm_pi2/g_norms[pi2_n_gfn-1];   }

  // Verbose output for debugging
  if (0 && my_rank==0) 
  {
    printf("\n - - - - - - - - - - - - - - - - - - - \n");

    if (1)
    {
      printf("\nGRID FUNCTION L-INFINITY NORMS:\n");
      printf("linfnorm_phi1=       %15.13g\n", g_norms[phi1_n_gfn-1]);
      printf("linfnorm_phi2=       %15.13g\n", g_norms[phi2_n_gfn-1]);
      printf("linfnorm_pi1=        %15.13g\n", g_norms[pi1_n_gfn-1]);
      printf("linfnorm_pi2=        %15.13g\n", g_norms[pi2_n_gfn-1]);
    }

    if (1)
    {
      printf("\nGRID FUNCTION L2 NORMS:\n");
      printf("l2norm_phi1= %15.13g\n",      l2norm_phi1);
      printf("l2norm_phi2= %15.13g\n",      l2norm_phi2);
      printf("l2norm_pi1=  %15.13g\n",      l2norm_pi1);
      printf("l2norm_pi2=  %15.13g\n",      l2norm_pi2);
    }
    
    if (1)
    {
      printf("\nTOTAL (NORMALIZED) L2 NORM:\nl2norm= %15.13g\n",l2norm);
    }

    printf("xmax =%15.13g\n",bbox[0]);
    printf("xmin =%15.13g\n",bbox[1]);
    printf("ymax =%15.13g\n",bbox[2]);
    printf("ymin =%15.13g\n",bbox[3]);
  }
  return l2norm;
}

//=============================================================================
// Returns some norm of the residual for the MG variables, *AND* 
// stores the point-wise residual in "f_res" for each MG variable "f" (for
// computing new RHS's). Required by the FAS multigrid algorithm to solve 
// elliptic equations
//=============================================================================
real qball_MG_residual(void)
{
   return 0;
}

//=============================================================================
// Performs 1 iteration of the evolution equations
//=============================================================================
void qball_evolve(int iter, int *ifc_mask)
{
   ldptr();

   // defined in updates.f
   update0_(modphi_np1,modphi_n,qdensity_np1,phi1_np1,phi1_n,phi2_np1,phi2_n,
            pi1_np1,pi1_n,pi2_np1,pi2_n,&Nx,&Ny,&dt,&dx,&dy,&B,
            &phys_bdy[1],&phys_bdy[0],&phys_bdy[3],&phys_bdy[2]);

   return;
}

//=============================================================================
// sets excision mask (NO ITERATOR, SO DON'T LOAD POINTERS!!!) This is meant to
// set points within an excised region defined by a grid function called 'mask'
// to some value 'excised', which is mostly useful for evolutions that contain
// black holes
//=============================================================================
void qball_fill_ex_mask(real *mask, real *mask_c, int dim, int *shape, int *shape_c, real *bbox, real excised)
{
}

//=============================================================================
// Called after each regridding
//=============================================================================
void qball_fill_bh_bboxes(real *bbox, int *num, int max_num)
{
}

//=============================================================================
// Called after each evolution step has been completed
//=============================================================================
void qball_post_tstep(int L)
{
   FILE *fp;
   double total;
   double mysum;
   int i;

   mysum=qtotcalcpamr_(qdensity_np1,qdensity_n,&Nx,&Ny,&dx,&dy,ghost_width);

   // find out total number of processes
   MPI_Comm_size(MPI_COMM_WORLD,&num_procs);

   // allocate an array of doubles of size num_procs
   recv = (double*) calloc(num_procs, sizeof(double));

   MPI_Gather(&mysum, 1, MPI_DOUBLE, recv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   if ( my_rank == 0 ) 
   {
     total = mysum; // add contribution from rank 0 processor
     for (i = 1; i < num_procs; i++) total = total + recv[num_procs-i];
     if ( L == 1 ) 
     {
       fp = fopen("./qtot.dat","a"); // print value to file
       //fprintf(fp,"%20.15f  %20.15f\n", PAMR_get_time(L), total);
       fclose(fp);
     }
   }

   return;
}

//=============================================================================
// Performs 1 relaxation sweep of the MG equations, and returns an estimate
// of the norm of the residual. Only relevant for multigrid
//=============================================================================
real qball_MG_relax(void)
{
   return 0;
}

//=============================================================================
// Computes the differential operator L acting on the current grid,
// in a region specified by the grid function "mask". Stores the result
// in "f_lop" for each MG variable "f". Only relevant for multigrid
//=============================================================================
void qball_L_op(void)
{
   return;
}

//=============================================================================
// Called after calculating the TRE for all variables. Allows you to modify
// the standard truncation error estimates for TRE_vars. After this is called,
// the pointwise l2-norm of the (potentially modified) TRE_vars is computed and
// used by AMRD to determine where to place child grids. See AMRD ref manual.
//=============================================================================
void qball_scale_tre(void)
{
   return;
}

//=============================================================================
// Called after regridding to reinitialize any time-independent (constant)
// functions that are not updated via evolution equations
//=============================================================================
void qball_post_regrid(void)
{
}

//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
// Main function
//#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#
int main(int argc, char **argv)
{
   amrd(argc,argv,&qball_id,&qball_var_pre_init,
        &qball_var_post_init, &qball_AMRH_var_clear,
        &qball_free_data, &qball_t0_cnst_data,
        &qball_evo_residual, &qball_MG_residual,
        &qball_evolve, &qball_MG_relax, &qball_L_op, 
        &qball_pre_io_calc, &qball_scale_tre, 
        &qball_post_regrid, &qball_post_tstep,
        &qball_fill_ex_mask, &qball_fill_bh_bboxes);
   if (my_rank==0) elapsed_time();
}

//=============================================================================
// Maintains/reports elapsed wall-clock time.
//=============================================================================
void elapsed_time(void) {
   static int    first = 1;
   struct        timeb t;
   static double msinit;
   double        mscurr, mselapsed;

   ftime(&t);
   mscurr = 1000.0 * t.time + t.millitm;
   if( first ) {
      msinit = mscurr;
      first = 0;
   }
	mselapsed = mscurr - msinit;
   printf("elapsed_time: Seconds since initial call: %12.3f\n",
         mselapsed / 1000.0);
}
