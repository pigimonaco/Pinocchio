/* ######HEADER###### */

#include "pinocchio.h"

int ThisTask,NTasks;

//int            pfft_flags_c2r, pfft_flags_r2c;
MPI_Comm        FFT_Comm;
internal_data   internal;

char *main_memory, *wheretoplace_mycat;
product_data *products, *frag;
unsigned int **seedtable;  // QUESTO RIMANE?
unsigned int   *cubes_ordering;
double **kdensity;
double **density;
double ***first_derivatives;
double ***second_derivatives;
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
