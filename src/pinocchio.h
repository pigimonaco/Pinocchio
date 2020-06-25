/*****************************************************************
 *                        PINOCCHIO  V4.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco
 Copyright (C) 2016
 
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <fftw3-mpi.h>

#define NYQUIST 1.
#define PI      3.14159265358979323846 
#define BLENGTH 200
#define SBLENGTH 100
#define GBYTE 1073741824.0
#define MBYTE 1048576.0
#define alloc_verbose 0
#define MAXOUTPUTS 100
#define SPEEDOFLIGHT ((double)299792.458)
#define GRAVITY ((double)4.30200e-9)      /*  (M_sun^-1 (km/s)^2 Mpc)  */

#if defined(THREE_LPT) && !defined(TWO_LPT)
#define TWO_LPT
#endif

#if defined(KEEP_DENSITY) && defined(THREE_LPT)
#error KEEP_DENSITY and THREE_LPT should not be used together
#endif

extern int ThisTask,NTasks;

typedef struct
{
  double Fmax;
  float Vel[3];
  int Rmax;
#ifdef TWO_LPT
  float Vel_2LPT[3];
#ifdef THREE_LPT
  float Vel_3LPT_1[3],Vel_3LPT_2[3];
#endif
#endif
#ifdef KEEP_DENSITY
  float Vel_after[3];
#ifdef TWO_LPT
  float Vel_2LPT_after[3];
#ifdef THREE_LPT
  float Vel_3LPT_1_after[3],Vel_3LPT_2_after[3];
#endif
#endif
#endif
} product_data;

extern void *main_memory, *wheretoplace_mycat;

extern product_data *products, *frag;

extern unsigned int **seedtable;

extern double **kdensity;
extern double **density;
extern double ***first_derivatives;
extern double ***second_derivatives;

extern double **VEL_for_displ;

#ifdef TWO_LPT
extern double *kvector_2LPT;
extern double *source_2LPT;
extern double **VEL2_for_displ;
#ifdef THREE_LPT
extern double *kvector_3LPT_1,*kvector_3LPT_2;
extern double *source_3LPT_1,*source_3LPT_2;
#endif
#endif

extern double Rsmooth;
typedef struct 
{
  int Nsmooth;
  double *Radius, *Variance, *TrueVariance;
} smoothing_data;
extern smoothing_data Smoothing;

extern int Ngrids;
typedef struct
{
  unsigned int GSglobal_x, GSglobal_y, GSglobal_z;
  unsigned int GSlocal_x, GSlocal_y, GSlocal_z;
  unsigned int GSstart_x, GSstart_y, GSstart_z;
  unsigned int GSlocal_k_x, GSlocal_k_y, GSlocal_k_z;
  unsigned int GSstart_k_x, GSstart_k_y, GSstart_k_z;
  unsigned int total_local_size,total_local_size_fft;
  unsigned int off;
  double lower_k_cutoff, upper_k_cutoff, norm, BoxSize, CellSize;
  fftw_plan forward_plan, reverse_plan;
} grid_data;
extern grid_data *MyGrids;

extern fftw_complex **cvector_fft;
extern double **rvector_fft;

#ifdef SCALE_DEPENDENT_GROWTH
typedef struct
{
  double ReferenceRedshift;
  int ReferenceOutput,Nkbins, NCAMB, ReferenceScale;
  char MatterFile[BLENGTH], TransferFile[BLENGTH], RunName[SBLENGTH], RedshiftsFile[BLENGTH];
  double *Logk, *LogPkref, D2ref, *Scalef, *RefGM;
} camb_data;
#endif

#ifdef SCALE_DEPENDENT
typedef struct
{
  int order;
  double redshift;
} ScaleDep_data;
extern ScaleDep_data ScaleDep;
#endif

typedef struct
{
  double Omega0, OmegaLambda, Hubble100, Sigma8, OmegaBaryon, DEw0, DEwa, 
    PrimordialIndex, InterPartDist, BoxSize, BoxSize_htrue, BoxSize_h100, ParticleMass, 
    StartingzForPLC, LastzForPLC, InputSpectrum_UnitLength_in_cm, WDM_PartMass_in_kev, 
    BoundaryLayerFactor, Largest, MaxMemPerParticle, k_for_GM, PLCAperture,
    PLCCenter[3], PLCAxis[3];
  char RunFlag[SBLENGTH],DataDir[SBLENGTH],TabulatedEoSfile[BLENGTH],ParameterFile[BLENGTH],
    OutputList[BLENGTH],FileWithInputSpectrum[BLENGTH];
  int GridSize[3],WriteRmax, WriteFmax, WriteVmax,
    CatalogInAscii, DoNotWriteCatalogs, DoNotWriteHistories, WriteSnapshot, 
    OutputInH100, RandomSeed, MaxMem, NumFiles, MinMassForCat, 
    BoxInH100, simpleLambda, AnalyticMassFunction, MinHaloMass, PLCProvideConeData;
#ifdef SCALE_DEPENDENT_GROWTH
  camb_data camb;
#endif
} param_data;
extern param_data params;

typedef struct
{
  int n;
  double F[MAXOUTPUTS],z[MAXOUTPUTS],zlast,Flast;
} output_data;
extern output_data outputs;


typedef struct
{
  unsigned int safe, Npart;
  int nbox_x,  nbox_y,  nbox_z_thisslice, nbox_z_allslices;
  int mybox_x, mybox_y, mybox_z;
  int Lgrid_x, Lgrid_y, Lgrid_z; 
  int Lgwbl_x, Lgwbl_y, Lgwbl_z; 
  int start_x, start_y, start_z;
  int stabl_x, stabl_y, stabl_z;
  int safe_x,  safe_y,  safe_z;
  int pbc_x,   pbc_y,   pbc_z;
  double SafetyBorder,overhead;
} subbox_data;
extern subbox_data subbox;

typedef struct
{
  double init,total,dens,fft,coll,vel,lpt,fmax,distr,sort,group,frag,io
#ifdef PLC
    ,plc
#endif
    ;
} cputime_data;
extern cputime_data cputime;

extern int WindowFunctionType;

typedef struct
{
  int Mass;
  double Pos[3],Vel[3];
#ifdef TWO_LPT
  double Vel_2LPT[3];
#ifdef THREE_LPT
  double Vel_3LPT_1[3], Vel_3LPT_2[3];
#endif
#endif
#ifdef KEEP_DENSITY
  double Vel_after[3];
#ifdef TWO_LPT
  double Vel_2LPT_after[3];
#ifdef THREE_LPT
  double Vel_3LPT_1_after[3], Vel_3LPT_2_after[3];
#endif
#endif
#endif
  int ll, halo_app, mass_at_merger, merged_with, point, bottom, good;
  double t_appear, t_peak, t_merge;
  unsigned long long int name;
  int trackT,trackC;
#ifdef PLC
  double Flast;
#endif
} group_data;
extern group_data *groups;

#ifdef PLC
typedef struct
{
  int i,j,k;
  double F1,F2;
} replication_data;

typedef struct
{
  int Nreplications, Nmax, Nstored, Nhalotot;
  double Fstart,Fstop,center[3];
  double xvers[3],yvers[3],zvers[3];
  replication_data *repls;
} plc_data;
extern plc_data plc;

typedef struct
{
  int Mass;
  unsigned long long int name;
  double z,x[3],v[3],rhor,theta,phi;
} plcgroup_data;
extern plcgroup_data *plcgroups;
#endif

extern char date_string[25];

extern int *indices,*group_ID,*linking_list;
extern int NSlices,ThisSlice;

/* fragmentation parameters */
extern double f_m, f_rm, espo, f_a, f_ra, f_200, sigmaD0;

#ifdef SCALE_DEPENDENT_GROWTH
typedef struct
{
  int flag, ismooth;
  double radius;
} SDGM_data;
extern SDGM_data SDGM;
#endif

#define NWINT 1000
extern gsl_integration_workspace *workspace;
extern gsl_rng *random_generator;

typedef struct
{
  unsigned long long int name;
  int nick, ll, mw, mass, mam;
  double zme, zpe, zap;
}  histories_data;

#define DELTAM 0.08

typedef struct
{
  int NBIN;
  double mmin,mmax,vol,hfactor,hfactor4;
  int *ninbin,*ninbin_local;
  double *massinbin,*massinbin_local;
} mf_data;
extern mf_data mf;


/* splines for interpolations */
extern gsl_spline **SPLINE;
extern gsl_interp_accel **ACCEL;

#ifdef MOD_GRAV_FR
extern double H_over_c;
#endif

/* prototypes for functions defined in collapse_times.c */
int compute_collapse_times(int);
int compute_velocities(int);
#ifdef TABULATED_CT
int initialize_collapse_times(int);
int reset_collapse_times(int);
#endif

/* prototypes for functions defined in fmax-fftw.c */
int set_one_grid(int);
int initialize_fft();
double forward_transform(int);
double reverse_transform(int);
int finalize_fftw();
int compute_derivative(int, int, int);
void write_in_cvector(int, double *);
void write_from_cvector(int, double *);
void write_in_rvector(int, double *);
void write_from_rvector(int, double *);

/* prototypes for functions defined in allocations.c */
int allocate_main_memory(void);
int deallocate_fft_vectors(int);
int reallocate_memory_for_fragmentation_1();
int reallocate_memory_for_fragmentation_2(int);


/* prototypes for functions defined in GenIC.c */
int GenIC(int);

/* prototypes for functions defined in initialization.c */
int initialization();
int find_start(int, int, int);
int find_length(int, int, int);
int set_parameters(void);
int set_grids(void);

/* prototypes in write_fields.c */
int write_fields(void);
int write_density(int);
int write_product(int, char*);

/* prototypes in write_snapshot.c */
int write_snapshot(int);
int write_LPT_snapshot(double);

/* prototypes for functions defined in cosmo.c */
int initialize_cosmology();
int initialize_MassVariance();
double OmegaMatter(double);
double OmegaLambda(double);
double Hubble(double);
double Hubble_Gyr(double);
double fomega(double,double);
double fomega_2LPT(double,double);
double fomega_3LPT_2(double,double);
double fomega_3LPT_1(double,double);
double CosmicTime(double);
double InverseCosmicTime(double);
double GrowingMode(double,double);
double GrowingMode_2LPT(double,double);
double GrowingMode_3LPT_1(double,double);
double GrowingMode_3LPT_2(double,double);
double InverseGrowingMode(double);
double ComovingDistance(double);
double InverseComovingDistance(double);
double dComovingDistance_dz(double);
double PowerSpectrum(double);
double MassVariance(double);
double dMassVariance_dr(double);
double DisplVariance(double);
double Radius(double);
double SizeForMass(double);
double MassForSize(double);
//double TypicalCollapsingMass(double);
double dOmega_dVariance(double, double);
double AnalyticMassFunction(double, double);
double WindowFunction(double);
double my_spline_eval(gsl_spline *, double, gsl_interp_accel *);
int jac(double, const double [], double *, double [], void *);

/* prototypes for functions defined in ReadParamFile.c */
int read_parameter_file();

/* prototypes for functions defined in fmax.c */
int compute_fmax(void);
int compute_displacements(void);
char *fdate(void);

#ifdef TWO_LPT
/* prototypes for functions defined in LPT.c */
int compute_LPT_displacements(int);
#endif

/* prototypes for functions defined in distribute.c */
int distribute(void);

/* prototypes for functions defined in fragment.c */
int fragment(void);

/* prototypes for functions defined in build_groups.c */
int build_groups(int,double);

#ifdef SCALE_DEPENDENT_GROWTH
int read_power_table_from_CAMB(void);
int initialize_ScaleDependentGrowth(void);
double MatterGrowingMode(double);
double VelGrowingMode(double);
double InverseMatterGrowingMode(double);
double VelfOmega(double);
double VelfOmega2(double);
#endif

#ifdef WHITENOISE
int read_white_noise(void);
#endif
