/* ######HEADER###### */

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
#include <pfft.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef USE_GPERFTOOLS
#include <gperftools/profiler.h>
#endif

/* this library is used to vectorize the computation of collapse times */
/* #if !(defined(__aarch64__) || defined(__arm__)) */
/* #include <immintrin.h> */
/* #endif */

/* Defines */
#define NYQUIST 1.
#define PI      3.14159265358979323846 
#define LBLENGTH 400
#define SBLENGTH 100
#define GBYTE 1073741824.0
#define MBYTE 1048576.0
#define alloc_verbose 0
#define MAXOUTPUTS 100
#define SPEEDOFLIGHT ((double)299792.458) /* km/s */
#define GRAVITY ((double)4.30200e-9)      /*  (M_sun^-1 (km/s)^2 Mpc)  */
#define NBINS 210     /* number of time bins in cosmological quantities */
#define FRAGFIELDS 6

#define NSIGMA ((double)6.0)
#define STEP_VAR ((double)0.3)  //0.2)   /* this sets the spacing for smoothing radii */

#define NV 6
#define FILAMENT 1
#define SHIFT 0.5

#define ORDER_FOR_GROUPS 2
#define ORDER_FOR_CATALOG 3

#define ALIGN 32     /* for memory alignment */
#define UINTLEN 32   /* 8*sizeof(unsigned int) */

/* these templates define how to pass from coordinates to indices */
#define INDEX_TO_COORD(I,X,Y,Z,L) ({Z=(I)%L[_z_]; int _KK_=(I)/L[_z_]; Y=_KK_%L[_y_]; X=_KK_/L[_y_];})
#define COORD_TO_INDEX(X,Y,Z,L) ((Z) + L[_z_]*((Y) + L[_y_]*(X)))

/* coordinates */
#define _x_ 0
#define _y_ 1
#define _z_ 2

#define DECOMPOSITION_LIMIT_FACTOR_2D 1    /* smallest allowed side lenght of rectangular pencils in */
                                           /* 2D decomposition of FFT */

/* debug levels */
#define dprintf(LEVEL, TASK, ...) do{if( ((LEVEL) <= internal.verbose_level) && (ThisTask == (TASK))) fprintf(stdout, __VA_ARGS__);} while(1 == 0)
#define VDBG  4   // verbose level for debug
#define VDIAG 2   // verbose level for diagnostics
#define VMSG  1   // verbose level for flow messages
#define VXX   0   // essential messages
#define VERR  VXX // non letal errors
#define VXERR -1  // letal errors

#define SWAP_INT( A, B ) (A) ^= (B), (B) ^= (A), (A) ^= (B);

/* checks of compiler flags */
#if defined(THREE_LPT) && !defined(TWO_LPT)
#define TWO_LPT
#endif

#if !defined(ELL_SNG) && !defined(ELL_CLASSIC)
#define ELL_CLASSIC
#endif

#if (defined(READ_PK_TABLE) || defined(MOD_GRAV_FR)) && !defined(SCALE_DEPENDENT)
#define SCALE_DEPENDENT
#endif

#if defined(READ_PK_TABLE) && defined(MOD_GRAV_FR)
#error REAK_PK_TABLE and MOD_GRAV_FR cannot be chosen together
#endif

#if defined(MOD_GRAV_FR) && !defined(FR0)
#error Please set a value to FR0 when you choose MOD_GRAV_FR
#endif


/* vectorialization */
#define DVEC_SIZE 4

typedef double dvec __attribute__ ((vector_size (DVEC_SIZE*sizeof(double))));
typedef long int ivec __attribute__ ((vector_size (DVEC_SIZE*sizeof(long int))));

typedef union
{
  dvec   V;
  double v[DVEC_SIZE];
} dvec_u;

typedef union
{
  ivec V;
  int  v[DVEC_SIZE];
} ivec_u;


/* variables and type definitions */
extern int ThisTask,NTasks;
/* pfft-related variables */
/* extern int pfft_flags_c2r, pfft_flags_r2c; */
extern MPI_Comm FFT_Comm;


typedef struct
{
  int tasks_subdivision_dim;            /* 1, 2 or 3 to divide in slabs, pencils and volumes */
  int tasks_subdivision_3D[4];          /* ??? */
  int constrain_task_decomposition[3];  /* constraints on the number of subdivisions for each dimension */
  int verbose_level;                    /* for dprintf */
  int mimic_original_seedtable;         /* logical, set to 1 to reproduce exactly GenIC */
  int dump_vectors;                     /* logical, dump vectors to files */
  int dump_seedplane;                   /* logical, dump seedplane to files */
  int dump_kdensity;                    /* logical, dump Fourier-space density to files */
  int large_plane;                      /* select the new generation of ICs */
  int nthreads_omp;                     /* number of OMP threads */
  int nthreads_fft;                     /* number of FFT threads */
} internal_data;
extern internal_data internal;

typedef unsigned int uint;
//typedef unsigned long long int UL;  // mi pare non ci sia


#ifdef DOUBLE_PRECISION_PRODUCTS
#define MPI_PRODFLOAT MPI_DOUBLE
typedef double PRODFLOAT;
#else
#define MPI_PRODFLOAT MPI_FLOAT
typedef float PRODFLOAT;
#endif

typedef struct  // RIALLINEARE?
{
  int Rmax;
  PRODFLOAT Fmax,Vel[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1[3],Vel_3LPT_2[3];
#endif
#endif

#ifdef SNAPSHOT
  PRODFLOAT zacc;
  int group_ID;
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT Vel_after[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT_after[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1_after[3],Vel_3LPT_2_after[3];
#endif
#endif
#endif

} product_data __attribute__((aligned (ALIGN)));  // VERIFICARE

extern char *main_memory, *wheretoplace_mycat;

extern product_data *products, *frag;

extern unsigned int *cubes_ordering;

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
#ifdef SCALE_DEPENDENT
  double *Rad_GM, *k_GM_dens, *k_GM_displ, *k_GM_vel;
#endif
} smoothing_data;
extern smoothing_data Smoothing;

extern int Ngrids;
typedef struct
{
  unsigned int       total_local_size, total_local_size_fft;
  unsigned int       off, ParticlesPerTask;
  ptrdiff_t          GSglobal[3];
  ptrdiff_t          GSlocal[3];
  ptrdiff_t          GSstart[3];
  ptrdiff_t          GSlocal_k[3];
  ptrdiff_t          GSstart_k[3];
  double             lower_k_cutoff, upper_k_cutoff, norm, BoxSize, CellSize;
  pfft_plan          forward_plan, reverse_plan;
  unsigned long long Ntotal;
} grid_data;
extern grid_data *MyGrids;

extern pfft_complex **cvector_fft;
extern double **rvector_fft;

#ifdef READ_PK_TABLE
typedef struct
{
  int Nkbins, NCAMB;
  char MatterFile[SBLENGTH], TransferFile[SBLENGTH], RunName[SBLENGTH], RedshiftsFile[LBLENGTH];
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
    BoundaryLayerFactor, Largest, MaxMemPerParticle, k_for_GM, PredPeakFactor, PLCAperture,
    PLCCenter[3], PLCAxis[3];
  char RunFlag[SBLENGTH],DumpDir[SBLENGTH],TabulatedEoSfile[LBLENGTH],ParameterFile[LBLENGTH],
    OutputList[LBLENGTH],FileWithInputSpectrum[LBLENGTH],CTtableFile[LBLENGTH];
  int GridSize[3],WriteProducts,WriteDensity,DumpProducts,ReadProductsFromDumps,
    CatalogInAscii, DoNotWriteCatalogs, DoNotWriteHistories, WriteTimelessSnapshot,
    OutputInH100, RandomSeed, MaxMem, NumFiles, 
    BoxInH100, simpleLambda, AnalyticMassFunction, MinHaloMass, PLCProvideConeData, ExitIfExtraParticles,
    use_transposed_fft, use_inplace_fft;
#ifdef READ_PK_TABLE
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
  unsigned int Npart, Ngood, Nstored, PredNpeaks, maplength;
  unsigned int Nalloc, Nneeded;
  int nbox[3];
  int mybox[3];
  int Lgrid[3]; 
  int Lgwbl[3]; 
  int start[3];
  int stabl[3];
  int safe[3];
  int pbc[3];
  double SafetyBorder,overhead;
} subbox_data;
extern subbox_data subbox;

typedef struct
{
  double init,total, dens, fft, coll, invcoll, ell, vel, lpt, fmax, distr, sort, group, frag, io,
    deriv, mem_transf, partial, set_subboxes, set_plc, memory_allocation, fft_initialization
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
  PRODFLOAT Pos[3],Vel[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1[3], Vel_3LPT_2[3];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT Vel_after[3];
#ifdef TWO_LPT
  PRODFLOAT Vel_2LPT_after[3];
#ifdef THREE_LPT
  PRODFLOAT Vel_3LPT_1_after[3], Vel_3LPT_2_after[3];
#endif
#endif
#endif
  int ll, halo_app, mass_at_merger, merged_with, point, bottom, good;
  PRODFLOAT t_appear, t_peak, t_merge;
  unsigned long long int name;
  int trackT,trackC;
#ifdef PLC
  PRODFLOAT Flast;
#endif
} group_data;
extern group_data *groups;

#ifdef PLC
typedef struct
{
  int i,j,k;
  PRODFLOAT F1,F2;
} replication_data;

typedef struct
{
  int Nreplications, Nmax, Nstored, Nstored_last, Nhalotot;
  double Nexpected;
  double Fstart,Fstop,center[3];
  double xvers[3],yvers[3],zvers[3];
  replication_data *repls;
  int nzbins;
  double delta_z,*nz;
} plc_data;
extern plc_data plc;

typedef struct
{
  int Mass;
  unsigned long long int name;
  PRODFLOAT z,x[3],v[3]; //,rhor,theta,phi;
} plcgroup_data;
extern plcgroup_data *plcgroups;
#endif

extern char date_string[25];

extern int *frag_pos,*indices,*indicesY,*sorted_pos,*group_ID,*linking_list;

extern unsigned int *frag_map, *frag_map_update;

/* fragmentation parameters */
extern double f_m, f_rm, espo, f_a, f_ra, f_200, sigmaD0;

#define NWINT 1000
extern gsl_integration_workspace *workspace;
extern gsl_rng *random_generator;
#define TOLERANCE ((double)1.e-4)

typedef struct
{
  unsigned long long int name;
  int nick, ll, mw, mass, mam;
  PRODFLOAT zme, zpe, zap;
}  histories_data;

#define DELTAM 0.05

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
#if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC)
extern gsl_spline **SPLINE_INVGROW;
extern gsl_interp_accel **ACCEL_INVGROW;
#endif

#ifdef MOD_GRAV_FR
extern double H_over_c;
#endif

typedef struct
{
  size_t prods, fields, fields_to_keep, fft, first_allocated, fmax_total, 
    frag_prods, frag_arrays, groups, frag_allocated, frag_total, all_allocated, all;
} memory_data;
extern memory_data memory;

extern int ngroups;

typedef struct
{
  int M,i;
  PRODFLOAT R,q[3],v[3],D,Dv;
  double z;
#ifdef TWO_LPT
  PRODFLOAT D2,D2v,v2[3];
#ifdef THREE_LPT
  PRODFLOAT D31,D31v,v31[3],D32,D32v,v32[3];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
  PRODFLOAT w;
  PRODFLOAT v_aft[3];
#ifdef TWO_LPT
  PRODFLOAT v2_aft[3];
#ifdef THREE_LPT
  PRODFLOAT v31_aft[3],v32_aft[3];
#endif
#endif
#endif
} pos_data;
extern pos_data obj, obj1, obj2;

typedef struct
{
  unsigned long long int name;
  PRODFLOAT M,x[3],v[3];
  PRODFLOAT q[3];
#ifndef LIGHT_OUTPUT
  int n;
  int pad;
#endif
}  catalog_data;

typedef struct
{
  unsigned long long int name;
#ifndef LIGHT_OUTPUT
  PRODFLOAT red,x,y,z,vx,vy,vz,Mass,theta,phi,v_los,obsz;
#else
  PRODFLOAT red,Mass,theta,phi,obsz;
#endif
} plc_write_data;


typedef struct
{
  int n, mine;
  double z[MAXOUTPUTS];
} Segment_data;
extern Segment_data Segment;

/* prototypes for functions defined in collapse_times.c */
int compute_collapse_times(int);
#ifdef TABULATED_CT
int initialize_collapse_times(int, int);
int reset_collapse_times(int);
#endif

/* prototypes for functions defined in fmax-fftw.c */
int set_one_grid(int);
int compute_fft_plans();
double forward_transform(int);
double reverse_transform(int);
int finalize_fft();
int compute_derivative(int, int, int);
void write_in_cvector(int, double *);
void write_from_cvector(int, double *);
void write_in_rvector(int, double *);
void write_from_rvector(int, double *);
int store_velocities();

// PROBABILMENTE DA TOGLIERE DOPO IL DEBUG
void dump_cvector(double*, int, int, ptrdiff_t *, ptrdiff_t *,  char *, int);
void dump_rvector(double*, int, ptrdiff_t *, ptrdiff_t *,  char *, int);

/* prototypes for functions defined in allocations.c */
int organize_main_memory(void);
int allocate_main_memory(void);
int deallocate_fft_vectors(int);
int reallocate_memory_for_fragmentation(void);
//int rearrange_memory(int);

/* prototypes for functions defined in GenIC.c */
int GenIC(int);
int GenIC_large(int);  // NE BASTA UNA?
//double VarianceOnGrid(int, double); //, double);

/* prototypes for functions defined in initialization.c */
int initialization();
int find_start(int, int, int);
int find_length(int, int, int);
int set_parameters(void);
int set_grids(void);
int set_subboxes(void);
int set_smoothing(void);
void greetings(void);
int check_parameters_and_directives(void);

/* prototypes in write_snapshot.c */
#ifdef SNAPSHOT
int write_density(int);
int write_LPT_snapshot(double);
int write_timeless_snapshot(void);
#endif

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
double InverseGrowingMode(double,int);
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
int dump_products(void);
int read_dumps(void);

#ifdef TWO_LPT
/* prototypes for functions defined in LPT.c */
int compute_LPT_displacements();
#endif

/* prototypes for functions defined in distribute.c */
int distribute(void);
int distribute_back(void);

/* prototypes for functions defined in fragment.c */
int fragment_driver(void);
int get_mapup_bit(unsigned int);
int get_map_bit(int, int, int);
void set_mapup_bit(int, int, int);
int estimate_file_size(void);
double compute_Nhalos_in_PLC(double, double);

/* prototypes for functions defined in build_groups.c */
int build_groups(int,double,int);
int quick_build_groups(int);
int update_map(unsigned int *);

// RIMETTERE LA LETTURA DEL WHITE NOISE
//#ifdef WHITENOISE
//int read_white_noise(void);
//#endif


/* fragmentation prototypes */
void condition_for_accretion(int, int, int, int, int, PRODFLOAT, int, double *, double *); // LEVARE primo argomento
void condition_for_merging(PRODFLOAT, int, int, int *);
void set_obj(int, PRODFLOAT, pos_data *);
void set_obj_vel(int, PRODFLOAT, pos_data *);
void set_point(int, int, int, int, PRODFLOAT, pos_data *);
void set_group(int, pos_data *);
PRODFLOAT q2x(int, pos_data *, int, double, int);
PRODFLOAT vel(int, pos_data *);
PRODFLOAT distance(int, pos_data *, pos_data *);
void clean_list(int *);
PRODFLOAT virial(int, PRODFLOAT, int);
void merge_groups(int, int, PRODFLOAT);
void update_history(int, int, PRODFLOAT);
void accretion(int, int, int, int, int, PRODFLOAT);
void update(pos_data *, pos_data *);
int write_catalog(int);
#ifdef ONLYFORFIRSTBHS
int write_histories(int);
#else
int write_histories(void);
#endif
int compute_mf(int);
int find_location(int, int, int);
#ifdef PLC
int write_PLC();
void coord_transformation_cartesian_polar(PRODFLOAT *, double *, double *, double *);
#endif

/* SCOREP */
#if defined (_SCOREP)
   /* remove inlining */
   #define FORCE_INLINE
#else
   #define FORCE_INLINE inline
#endif /* _SCOREP */
