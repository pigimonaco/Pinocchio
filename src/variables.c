/* ######HEADER###### */

#include "pinocchio.h"

int ThisTask,NTasks;

//int            pfft_flags_c2r, pfft_flags_r2c;
MPI_Comm        FFT_Comm;
internal_data   internal;

char *main_memory, *wheretoplace_mycat;
product_data *products, *frag;


#if defined(GPU_OMP)
int hostID;
int devID;
gpu_product_data gpu_products, host_products;
#elif defined(GPU_OMP_FULL)
int hostID;
int devID;
gpu_product_data gpu_products;
#endif // end GPU flags. We are working directly with the data structure allocated and initilialized withing Pinocchio 

unsigned int **seedtable;  // QUESTO RIMANE?
unsigned int   *cubes_ordering;
double **kdensity;
double **density;
double ***first_derivatives;
double ***second_derivatives;
#if defined(GPU_OMP) 
gpu_second_derivatives_data gpu_second_derivatives;
#endif // GPU_OMP
double **VEL_for_displ;

#ifdef TWO_LPT
double *kvector_2LPT;
double *source_2LPT;
double **VEL2_for_displ;
#ifdef THREE_LPT
double *kvector_3LPT_1,*kvector_3LPT_2;
double *source_3LPT_1,*source_3LPT_2;
#endif
#endif

double Rsmooth;
smoothing_data Smoothing;

grid_data *MyGrids;
int Ngrids;

pfft_complex **cvector_fft;
double **rvector_fft;

param_data params={0};
output_data outputs;
subbox_data subbox;
#ifdef PLC
plc_data plc;
plcgroup_data *plcgroups;
#endif

cputime_data cputime={0.0};

#if defined(GPU_OMP) || defined(GPU_OMP_FULL)
gputime_data gputime;
#endif // GPU_OMP || GPU_OMP_FULL

int WindowFunctionType;

group_data *groups;

char date_string[25];

int *frag_pos,*indices,*indicesY,*sorted_pos,*group_ID,*linking_list;
unsigned int *frag_map, *frag_map_update;
double f_m, f_rm, espo, f_a, f_ra, f_200, sigmaD0;

#ifdef SCALE_DEPENDENT
ScaleDep_data ScaleDep;
#endif

gsl_integration_workspace * workspace;
gsl_rng *random_generator;

mf_data mf;

gsl_spline **SPLINE;
gsl_interp_accel **ACCEL;

#if defined(CUSTOM_INTERPOLATION) || defined(GPU_OMP)
CubicSpline *host_spline;
#endif

#if defined(GPU_OMP)
CubicSpline gpu_spline;
#elif defined(GPU_OMP_FULL)
CubicSpline *gpu_spline;
#endif

#if defined(TABULATED_CT) && !defined(CUSTOM_INTERPOLATION) && !defined(GPU_OMP_FULL)
gsl_spline ***CT_Spline;
#elif defined(TABULATED_CT) && (defined(CUSTOM_INTERPOLATION) || defined(GPU_OMP_FULL))
CubicSpline **CT_Spline;
#endif

#if defined(TABULATED_CT)
int    Ncomputations, start, length;
double *CT_table = 0x0;
double bin_x;
double *delta_vector;
gsl_interp_accel *accel = 0x0;
FILE *CTtableFilePointer = NULL;
#endif


#if defined(SCALE_DEPENDENT) && defined(ELL_CLASSIC)
gsl_spline **SPLINE_INVGROW;
gsl_interp_accel **ACCEL_INVGROW;
#endif

#ifdef MOD_GRAV_FR
double H_over_c;
#endif

memory_data memory;

int ngroups;
pos_data obj, obj1, obj2;
Segment_data Segment;
