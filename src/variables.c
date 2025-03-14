/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan, 
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025
 
 github: https://github.com/pigimonaco/Pinocchio
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
int map_to_be_used;
double f_m, f_rm, espo, f_a, f_ra, f_200, sigmaD0;

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
extern pos_data obj, obj1, obj2;
ScaleDep_data ScaleDep;
