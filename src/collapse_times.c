/*****************************************************************
 *                        PINOCCHI0  V4.0                        *
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

// LUCA: qui bisogna rimettere a posto la vettorializzazione, che pero` va spenta
//       se si usano le tabelle
// LUCA: inoltre, va ompizzato (e vettorializzato?) il calcolo delle tabelle di collapse times


#include "pinocchio.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

// LUCA: per il momento ho tolto l'inline, non mi compilava
#ifdef VETTORIALIZZA
#ifndef TABULATED_CT
int     v_inverse_collapse_time ( int, dvec dtensor[6], dvec_u , dvec , dvec_u *, dvec_u *, dvec_u *, dvec_u *); //__attribute__((always_inline));
#endif
#endif
double  inverse_collapse_time   ( int, double *, double *, double *, double *, int *); //__attribute__((always_inline));
double  ell                     ( int, double, double, double); //__attribute__((always_inline));
double  ell_classic             ( int, double, double, double); //__attribute__((always_inline));
double  ell_sng                 ( int, double, double, double); //__attribute__((always_inline));

void ord(double *,double *,double *);

int ismooth;

#if defined( _OPENMP )

double cputime_ell;
#pragma omp threadprivate(cputime_ell)
int fails;
#pragma omp threadprivate(fails)
#else

#define cputime_ell       cputime.ell

#endif

// PER LUCA: controlliamo che sto aggiungendo bene le funzioni
#ifdef TABULATED_CT
double  interpolate_collapse_time(double, double, double) //__attribute__((always_inline));

// quanto segue alla fine potra` essere spostato giu` dove inizia la parte TABULATED_CT

//#define DEBUG
//#define PRINTJUNK

//#define ALL_SPLINE
//#define TRILINEAR
//#define HISTO
#define BILINEAR_SPLINE

#define CT_NBINS_XY (50)  // 50
#define CT_NBINS_D (100)  // 100
#define CT_SQUEEZE (1.2)  // 1.2
#define CT_EXPO   (1.75)  // 1.75
#define CT_RANGE_D (7.0)  // 7.0
#define CT_RANGE_X (3.5)
#define CT_DELTA0 (-1.0)
double *CT_table=0x0;
double bin_x;
static gsl_interp_accel *accel=0x0;
gsl_spline ***CT_Spline;
double *delta_vector;
int Ncomputations, start, length;

#ifdef TRILINEAR
#ifdef HISTO
int *CT_histo=0x0;
#endif
#endif

#ifdef DEBUG
double ave_diff=0.0,var_diff=0.0;
int counter=0;
#ifdef PRINTJUNK
FILE *JUNK;
#endif
#endif

#endif /* TABULATED_CT */


int compute_collapse_times(int ismooth)
{

  double local_average,local_variance;
  dvec_u  vlocal_average, vlocal_variance;
  dvec    vzero = {0, 0, 0, 0};

  /* initialize the variance */
  local_variance = 0.0;
  local_average  = 0.0;
  vlocal_average.V  = vzero;
  vlocal_variance.V = vzero; 

  /* initialize Fmax if it is the first smoothing */
  if (!ismooth)
#if defined( _OPENMP )
#pragma omp parallel for
#endif
    for (int i = 0; i < MyGrids[0].total_local_size; i++)
      {
	products[i].Fmax   = -10.0;
	products[i].Rmax   = -1;
	products[i].Vel[0] = 0.0;
	products[i].Vel[1] = 0.0;
	products[i].Vel[2] = 0.0;
#ifdef TWO_LPT
	products[i].Vel_2LPT[0] = 0.0;
	products[i].Vel_2LPT[1] = 0.0;
	products[i].Vel_2LPT[2] = 0.0;
#ifdef THREE_LPT
	products[i].Vel_3LPT_1[0] = 0.0;
	products[i].Vel_3LPT_1[1] = 0.0;
	products[i].Vel_3LPT_1[2] = 0.0;
	products[i].Vel_3LPT_2[0] = 0.0;
	products[i].Vel_3LPT_2[1] = 0.0;
	products[i].Vel_3LPT_2[2] = 0.0;
#endif
#endif
      }

  /* loop on all particles */

#if defined( _OPENMP )
  double invcoll_update = 0, ell_update = 0;
#pragma omp parallel
  {
    double  mylocal_average, mylocal_variance;
    double  cputime_invcoll;
    dvec_u  vmylocal_average, vmylocal_variance;
    dvec    vzero = {0, 0, 0, 0};
    int     fails = 0;
    int     pid = omp_get_thread_num();
    
    /* initialize the variance */
    mylocal_average     = 0;
    mylocal_variance    = 0;
    cputime_invcoll     = 0;
    cputime_ell         = 0;
    vmylocal_average.V  = vzero;
    vmylocal_variance.V = vzero;

#else
    
#define mylocal_average   local_average
#define mylocal_variance  local_variance
#define vmylocal_average  vlocal_average
#define vmylocal_variance vlocal_variance
#define cputime_invcoll   cputime.invcoll

#endif

    /* NOTE: the first index to change is the x-index in the fragmentation part;
       but this is not relevant in this part of the code */
    /* for (int local_x = 0; local_x < MyGrids[0].GSlocal[_x_]; local_x++) */
    /*   { */
    /* 	int idx_x = local_x * MyGrids[0].GSlocal[_y_]; */
    /* 	for (int local_y = 0; local_y < MyGrids[0].GSlocal[_y_]; local_y++) */
    /* 	  { */
    /* 	    int idx_y = (idx_x + local_y) * MyGrids[0].GSlocal[_z_]; */
    /* 	    int mysize = DVEC_SIZE * MyGrids[0].GSlocal[_z_] / DVEC_SIZE; */
    /* 	    int local_z = 0; */


// LOOP NON VETTORIALIZZATO
    double tmp = MPI_Wtime();
#ifdef _OPENMP
#pragma omp for nowait
#endif
    for (int index=0; index<MyGrids[0].total_local_size; index++)
      {

	double diff_ten[6];
	for (int i = 0; i < 6; i++)
	  diff_ten[i] = second_derivatives[0][i][index];

	/* Computation of the variance of the linear density field */
	double delta = diff_ten[0]+diff_ten[1]+diff_ten[2];
	mylocal_average  += delta;
	mylocal_variance += delta*delta;

	/* Computation of ellipsoidal collapse */
	double lambda1,lambda2,lambda3;
	int    fail;
	double Fnew = inverse_collapse_time(ismooth, diff_ten, &lambda1, &lambda2, &lambda3, &fail);

	if (fail)
	  {
	    printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
	    fflush(stdout);
#if !defined( _OPENMP )		    		    
	    return 1;
#else
	    fails = 1;
#endif
	  }

	if (products[index].Fmax < Fnew)
	  {
	    products[index].Fmax = Fnew;
	    products[index].Rmax = ismooth;
	  }	      

      }
    
    cputime_invcoll += MPI_Wtime() - tmp;


//#ifdef VETTORIALIZZA
//#ifndef TABULATED_CT
//	    for ( ; local_z < mysize; local_z += DVEC_SIZE)
//	      {
//		int index = idx_y + local_z;
//		/* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */
//
//		/* sum all contributions */
//		dvec elements[6] __attribute__((aligned(32)));
//
//// LUCA: quella funzione non mi compila	      
////#if !(defined(__aarch64__) || defined(__arm__))
////		for (int i = 0; i < 6; i++)
////		  elements[i] = _mm256_load_pd(&second_derivatives[0][i][index]);
////#else
//		dvec_u ld_elements[6] __attribute__((aligned(32)));
//		for (int i = 0; i < 6; i++)
//		  {
//		    ld_elements[i].v[0] = second_derivatives[0][i][index ];
//		    ld_elements[i].v[1] = second_derivatives[0][i][index + 1];
//		    ld_elements[i].v[2] = second_derivatives[0][i][index + 2];
//		    ld_elements[i].v[3] = second_derivatives[0][i][index + 3];
//		    elements[i] = ld_elements[i].V;
//		  }
////#endif
//
//		/* Computation of the variance of the linear density field */
//		dvec_u delta __attribute__((aligned(32)));
//		delta.V = elements[0] + elements[1] + elements[2];
//		dvec delta_2 __attribute__((aligned(32)));
//		delta_2 = delta.V*delta.V;
//	      
//		vmylocal_average.V  += delta.V;
//		vmylocal_variance.V += delta_2;
//
//		/* Computation of ellipsoidal collapse */
//		dvec_u lambda1,lambda2,lambda3 __attribute__((aligned(32)));
//		dvec_u Fnew __attribute__((aligned(32)));
//		double tmp = MPI_Wtime();
//		if(v_inverse_collapse_time(ismooth, elements, delta, delta_2, &lambda1, &lambda2, &lambda3, &Fnew))
//		  {
//		    printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
//		    fflush(stdout);
//#if !defined( _OPENMP )		    
//		    return 1;
//#else
//		    fails = 1;
//#endif
//		  }
//		cputime_invcoll += MPI_Wtime() - tmp;
//
//		for(int ii = 0; ii < DVEC_SIZE; ii++)
//		  {
//
//		    if (products[index + ii].Fmax < Fnew.v[ii])
//		      {
//			products[index].Fmax = Fnew.v[ii];
//			products[index].Rmax = ismooth;
//		      }
//		  }
//	      }
//#endif
//#endif

	    /* /\* this is for all particles if vectorialization is switched off  */
	    /*    or for remaining particles otherwise *\/ */
	    /* for ( ; local_z < MyGrids[0].GSlocal[_z_]; local_z++) */
	    /*   { */
		/* int index = idx_y + local_z; */
   
    
#if defined( _OPENMP )

#pragma omp atomic
    local_variance += mylocal_variance;
#pragma omp atomic
    local_average += mylocal_average;    
    for ( int i = 0; i < DVEC_SIZE; i++)
      {
#pragma omp atomic
	vlocal_average.v[i]  += vmylocal_average.v[i];
#pragma omp atomic    
	vlocal_variance.v[i] += vmylocal_variance.v[i];
      }
    
#pragma omp atomic
    invcoll_update += cputime_invcoll;
#pragma omp atomic	
    ell_update     += cputime_ell;
#endif
        
#if defined( _OPENMP )
  }
#endif

  int all_fails = 0;
  
  /* Computes the obtained variance */

#if defined (_OPENMP)
#pragma omp atomic
  cputime.invcoll += invcoll_update/internal.nthreads_omp;
#pragma omp atomic	
  cputime.ell     += ell_update/internal.nthreads_omp;
#pragma omp atomic
  all_fails       += fails;
#endif

  if (all_fails)
    {
      printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
      fflush(stdout);
      return 1;
    }

    
  for(int ii = 0; ii < DVEC_SIZE; ii++)
    {
      local_variance += vlocal_variance.v[ii];
      local_average  += vlocal_average.v[ii];
    }
      
  double global_variance = 0.0;
  double global_average  = 0.0;

  global_variance = 0.0;
  global_average  = 0.0;

  MPI_Reduce(&local_variance, &global_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_average , &global_average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      global_variance /= (double)MyGrids[0].Ntotal;
      global_average  /= (double)MyGrids[0].Ntotal;
    }

  MPI_Bcast(&global_variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  Smoothing.TrueVariance[ismooth] = global_variance;

  return 0;
}


#ifdef VETTORIALIZZA
#ifndef TABULATED_CT
/* vectorialization is not used for the interpolation */
int v_inverse_collapse_time(int ismooth, dvec dtensor[6], dvec_u delta, dvec delta_2, dvec_u * restrict x1, dvec_u * restrict x2, dvec_u * restrict x3, dvec_u * restrict res) // inline
{

  dvec   mu2;
  dvec_u q, mu3;
  dvec vzero = {0.0, 0.0, 0.0, 0.0};  
  /* mu1 = delta, mu2 and mu3 are the principal invariants of the 3x3 tensor 
     of second derivatives. */

  /* diagonalization of the tensor */
  
  mu2 = 0.5 * delta_2;

  dvec dtensor_2[3];
  dtensor_2[0] = dtensor[0] * dtensor[0];
  dtensor_2[1] = dtensor[1] * dtensor[1];
  dtensor_2[2] = dtensor[2] * dtensor[2];
  mu2 -= 0.5 * (dtensor_2[0] + dtensor_2[1] + dtensor_2[2]);

  dtensor_2[0] = dtensor[3] * dtensor[3];
  dtensor_2[1] = dtensor[4] * dtensor[4];
  dtensor_2[2] = dtensor[5] * dtensor[5];
  mu2 -= dtensor_2[0] + dtensor_2[1] + dtensor_2[2];
    
  mu3.V = dtensor[0]*dtensor[1]*dtensor[2] +
    2.*dtensor[3]*dtensor[4]*dtensor[5] -
    dtensor[0] * dtensor_2[2] -
    dtensor[1] * dtensor_2[1] -
    dtensor[2] * dtensor_2[0];
  
  q.V = ( delta_2 - 3.0*mu2 ) /9.0;

  ivec_u compare;
  compare.V = (q.V < vzero);
  if ( compare.v[0] || compare.v[1] || compare.v[2] || compare.v[3] )
    {
      printf("\t***task %d: %g %g %g %g\n", ThisTask, q.v[0], q.v[1], q.v[2], q.v[3]);
      return 1;
    }

  dvec_u r, q3;

  r.V  = -(2.*delta_2*delta.V - 9.0*delta.V*mu2 + 27.0*mu3.V)/54.;
  q3.V = (q.V * q.V) * q.V;

  compare.V = (q3.V < r.V*r.V);
  if( compare.v[0] || compare.v[1] || compare.v[2] || compare.v[3] )
    {
      printf("\t+++task %d: %g %g %g %g\n", ThisTask, q3.v[0], q3.v[1], q3.v[2], q3.v[3]);
      return 1;
    }

  double inv_3 = 1.0/3.0;
  dvec   vinv_3 = {1.0/3.0, 1.0/3.0, 1.0/3.0, 1.0/3.0};
  dvec_u sq;
//#if !(defined(__aarch64__) || defined(__arm__) ) LUCA: questa chiamata non me la compila
//  sq.V = 2.0 * _mm256_sqrt_pd(q.V);
//#else
  sq.v[0] = 2.0 * sqrt(q.v[0]);
  sq.v[1] = 2.0 * sqrt(q.v[1]);
  sq.v[2] = 2.0 * sqrt(q.v[2]);
  sq.v[3] = 2.0 * sqrt(q.v[3]);
//#endif
  dvec_u delta_inv_3;
  delta_inv_3.V = delta.V * vinv_3;
  
  // the following assignations will be valid in case
  // the tensor is already diagonal (q == 0)
  (*x1).V = dtensor[0];
  (*x2).V = dtensor[1];
  (*x3).V = dtensor[2];
  
  compare.V = (q.V != vzero);
  for(int ii = 0; ii < DVEC_SIZE; ii++)
    if(compare.v[ii])
      {
	double t = acos(2*r.v[ii] / q.v[ii] / sq.v[ii]);
	(*x1).v[ii] = -sq.v[ii]*cos(t * inv_3) + delta_inv_3.v[ii];
	(*x2).v[ii] = -sq.v[ii]*cos((t + 2.*PI) * inv_3 )+ delta_inv_3.v[ii];
	(*x3).v[ii] = -sq.v[ii]*cos((t + 4.*PI) * inv_3 )+ delta_inv_3.v[ii];
      }

  /* ordering and inverse collapse time */
  //  double tmp = MPI_Wtime();
  for(int ii = 0; ii < DVEC_SIZE; ii++)
    {
      double _x1 = (*x1).v[ii];
      double _x2 = (*x2).v[ii];
      double _x3 = (*x3).v[ii];
      ord(&_x1, &_x2, &_x3);
#ifdef TABULATED_CT
      double t=interpolate_collapse_time(*x1,*x2,*x3);
#else
      double t = ell(ismooth, _x1, _x2, _x3);
#endif
      if(t > 0)
	(*res).v[ii] = 1.0 / t;
      else
	(*res).v[ii] = -10.0;
    }
  // cputime_ell += MPI_Wtime() - tmp;

  return 0;
}

#endif
#endif

double inverse_collapse_time(int ismooth, double * restrict deformation_tensor, double * restrict x1, double * restrict x2, double * restrict x3, int * restrict fail) // inline
{

  double mu1, mu2, mu3;
  double dtensor[6] = {deformation_tensor[0], deformation_tensor[1],
		       deformation_tensor[2], deformation_tensor[3],
		       deformation_tensor[4], deformation_tensor[5] };
  
  /* mu1, mu2 and mu3 are the principal invariants of the 3x3 tensor 
     of second derivatives. */

  *fail=0;

  /* diagonalization of the tensor */
  mu1 = dtensor[0] + dtensor[1] + dtensor[2];
  double mu1_2 = mu1 * mu1;
  
  mu2 = 0.5 * mu1_2;

  {
    double add[3];
    add[0] = dtensor[0]*dtensor[0];
    add[1] = dtensor[1]*dtensor[1];
    add[2] = dtensor[2]*dtensor[2];
    mu2 -= 0.5 * (add[0] + add[1] + add[2]);
  }

  double add[3];
  add[0] = dtensor[3]*dtensor[3];
  add[1] = dtensor[4]*dtensor[4];
  add[2] = dtensor[5]*dtensor[5];
  mu2 -= add[0] + add[1] + add[2];
    
  mu3 =  dtensor[0]*dtensor[1]*dtensor[2] +
    2.*dtensor[3]*dtensor[4]*dtensor[5] -
    dtensor[0]*add[2] -
    dtensor[1]*add[1] -
    dtensor[2]*add[0];

  double q;
  q = ( mu1_2 - 3.0*mu2 ) /9.0;

  if (q == 0.)
    {
      /* in this case the tensor is already diagonal */
      *x1 = dtensor[0];
      *x2 = dtensor[1];
      *x3 = dtensor[2];
    }
  else
    {      
      double r = -(2.*mu1_2*mu1 - 9.0*mu1*mu2 + 27.0*mu3)/54.;

      if (q*q*q < r*r || q<0.0)
	{
	  //*fail=1;
	  return -10.0;
	}

      double sq = 2 * sqrt(q);
      double t = acos(2*r/q/sq);
      double inv_3 = 1.0/3.0;
      *x1 = -sq*cos(t * inv_3) + mu1 * inv_3;
      *x2 = -sq*cos((t + 2.*PI) * inv_3 )+ mu1 * inv_3;
      *x3 = -sq*cos((t + 4.*PI) * inv_3 )+ mu1 * inv_3;
    }

  /* ordering and inverse collapse time */
  ord(x1, x2, x3);
  //double tmp = MPI_Wtime();
#ifdef TABULATED_CT
  double t = interpolate_collapse_time(*x1,*x2,*x3);

#ifdef DEBUG
  double t2=ell(ismooth,*x1,*x2,*x3);
  if (t2>0.9)
    {
      ave_diff+=t-t2;
      var_diff+=pow(t-t2,2.0);
      counter+=1;
#ifdef PRINTJUNK
      double ampl = sqrt(Smoothing.Variance[ismooth]);
      fprintf(JUNK,"%d  %f %f %f   %f %f %f   %20g %20g %20g\n",ismooth, (*x1+*x2+*x3)/ampl, (*x1-*x2)/ampl, (*x2-*x3)/ampl, *x1, *x2, *x3, t, t2, t-t2);
#endif
    }
#endif

#else
  double t = ell(ismooth,*x1,*x2,*x3);
#endif

  //cputime_ell += MPI_Wtime() - tmp;
  
  return t;
}



#ifdef TABULATED_CT

FILE *CTtableFilePointer;

int check_CTtable_header();
void write_CTtable_header();


int initialize_collapse_times(int ismooth, int onlycompute)
{
  double l1, l2, l3, del, ampl, x, y, interval, ref_interval, deltaf;
  int id,ix,iy,i,j, dummy, fail;
  char fname[LBLENGTH];

  /* the bin size for delta varies like (delta - CT_DELTA0)^CT_EXPO, 
     so it is smallest around CT_DELTA0 if CT_EXPO>1. The bin size 
     is never smaller than CT_MIN_INTERVAL.
     The delta vector is computed at the first smoothing radius. */

  if (!ismooth)
    {
      delta_vector=(double*)malloc(CT_NBINS_D * sizeof(double));

      /* sets the sampling in density delta */
      if (CT_EXPO==1)
	{
	  /* in this case the spacing is even */
	  interval=2.*CT_RANGE_D/(double)CT_NBINS_D;
	  for (id=0; id<CT_NBINS_D; id++)
	    delta_vector[id]=id*interval -CT_RANGE_D;
	}
      else
	{
	  /* sampling is finer around CT_DELTA0 */
	  deltaf = pow(CT_SQUEEZE/CT_EXPO,1./(CT_EXPO-1.));
	  if (CT_EXPO==2)
	    ref_interval=( ( log((CT_RANGE_D-CT_DELTA0)/deltaf) + log((CT_RANGE_D+CT_DELTA0)/deltaf) )/ 
			   CT_EXPO + 2.*deltaf/CT_SQUEEZE ) / (CT_NBINS_D-2.0);
	  else
	    ref_interval=( ( pow(CT_RANGE_D-CT_DELTA0,2.-CT_EXPO) + pow(CT_RANGE_D+CT_DELTA0,2.-CT_EXPO) 
			     - 2.*pow(deltaf,2.-CT_EXPO) )/ 
			   CT_EXPO/(2.-CT_EXPO) + 2.*deltaf/CT_SQUEEZE ) / (CT_NBINS_D-2.0);
	  id=0;
	  del=-CT_RANGE_D;
	  do 
	    {
	      delta_vector[id]=del;
	      interval = CT_EXPO * ref_interval * pow(fabs(del-CT_DELTA0),CT_EXPO-1.0);
	      interval = (interval/ref_interval<CT_SQUEEZE? ref_interval*CT_SQUEEZE : interval);
	      del+=interval;
	      id++;
	    }  
	  while (id<CT_NBINS_D);
	}

      if (!ThisTask)
	printf("[%s] Grid for interpolating collapse times: CT_NBINS_D=%d, CT_NBINS_XY=%d\n",fdate(),CT_NBINS_D,CT_NBINS_XY);

      /* computations are divided among tasks */
      Ncomputations = CT_NBINS_D * CT_NBINS_XY * CT_NBINS_XY;
      bin_x = CT_RANGE_X /(double)(CT_NBINS_XY);

      /* allocates spline for interpolations */
      accel = gsl_interp_accel_alloc();
      CT_Spline=(gsl_spline***)calloc(CT_NBINS_XY,sizeof(gsl_spline**));
      for (i=0; i<CT_NBINS_XY; i++)
	{
	  CT_Spline[i]=(gsl_spline**)calloc(CT_NBINS_XY,sizeof(gsl_spline*));
	  for (j=0; j<CT_NBINS_XY; j++)
	    CT_Spline[i][j]=gsl_spline_alloc(gsl_interp_cspline,CT_NBINS_D);
	}
      CT_table=(double *)calloc(Ncomputations, sizeof(double));
#ifdef TRILINEAR
#ifdef HISTO
      CT_histo=(int *)calloc(Ncomputations, sizeof(int));
#endif
#endif

#ifdef DEBUG
#ifdef PRINTJUNK
      char filename[50];
      sprintf(filename,"debug.task%d.txt",ThisTask);
      JUNK=fopen(filename,"w");
#endif
#endif

    }

  /* In this case collapse times are read from a file */
  if (strcmp(params.CTtableFile,"none") && !onlycompute)
    {
      /* opening of file and consistency tests */
      if (!ismooth)
	{
	  /* Task 0 reads the file */
	  if (!ThisTask)
	    {

	      CTtableFilePointer=fopen(params.CTtableFile,"r");
	      fail=check_CTtable_header(CTtableFilePointer);

	    }
	  MPI_Bcast(&fail, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
	  if (fail)
	    return 1;
	}

      if (!ThisTask)
	{
	  fread(&dummy,sizeof(int),1,CTtableFilePointer);
	  fread(CT_table,sizeof(double),Ncomputations,CTtableFilePointer);
	}

      MPI_Bcast(CT_table, Ncomputations, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (!ThisTask && ismooth==Smoothing.Nsmooth-1)
	fclose(CTtableFilePointer);

    }
  else

  /* The collapse time table is computed for each smoothing radius */
    {

      ampl = sqrt(Smoothing.Variance[ismooth]);

      /* map of computations of all tasks */
      int *counts=(int*)malloc(NTasks*sizeof(int));
      int *displs=(int*)malloc(NTasks*sizeof(int));
      int bunch = Ncomputations / NTasks;
      int remainder = Ncomputations % NTasks;

      for (i=0, ix=0; i<NTasks; i++)
	{
	  counts[i]=bunch + (i < remainder);
	  displs[i]=bunch*i + ((remainder > i)? i : remainder );
	}

      /* computations for this task */
      for (i=displs[ThisTask]; i<displs[ThisTask]+counts[ThisTask]; i++)
	{
	  id=i%CT_NBINS_D;
	  ix=(i/CT_NBINS_D)%CT_NBINS_XY;
	  iy=i/CT_NBINS_D/CT_NBINS_XY;

	  x   = ix * bin_x;
	  y   = iy * bin_x;
	  l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
	  l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
	  l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

	  CT_table[i]= ell(l1,l2,l3);
	}

      MPI_Allgatherv( MPI_IN_PLACE, counts[ThisTask], MPI_DOUBLE, CT_table, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

      free(displs);
      free(counts);

      /* Writes the table on a file */
      if (!ThisTask)
	{
	  if (onlycompute)
	    strcpy(fname,params.CTtableFile);
	  else
	    sprintf(fname,"pinocchio.%s.CTtable.out",params.RunFlag);
	  if (!ismooth)
	    {
	      CTtableFilePointer=fopen(fname,"w");
	      write_CTtable_header(CTtableFilePointer);
	    }
	  else
	    CTtableFilePointer=fopen(fname,"a");
#ifdef ASCII
	  fprintf(CTtableFilePointer,"################\n# smoothing %2d #\n################\n",ismooth);

	  for (i=0; i<Ncomputations; i++)
	    {
	      id=i%CT_NBINS_D;
	      ix=(i/CT_NBINS_D)%CT_NBINS_XY;
	      iy=i/CT_NBINS_D/CT_NBINS_XY;

	      x   = ix * bin_x;
	      y   = iy * bin_x;
	      l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
	      l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
	      l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

	      fprintf(CTtableFilePointer," %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20g\n",
		      ismooth,id,ix,iy,delta_vector[id],x,y,l1,l2,l3,CT_table[i]);

	    }
#else
	  fwrite(&ismooth,sizeof(int),1,CTtableFilePointer);
	  fwrite(CT_table,sizeof(double),Ncomputations,CTtableFilePointer);
#endif

	  fclose(CTtableFilePointer);
	}

    }

  /* now initialize the splines */
  for (i=0; i<CT_NBINS_XY; i++)
    for (j=0; j<CT_NBINS_XY; j++)
	gsl_spline_init(CT_Spline[i][j],delta_vector,&CT_table[i*CT_NBINS_D+j*CT_NBINS_D*CT_NBINS_XY],CT_NBINS_D);
  
  return 0;
}


int reset_collapse_times(int ismooth)
{
  

#ifdef DEBUG
  if (counter)
    {
#ifdef TRILINEAR
      int flag=1;
#endif
#ifdef BILINEAR_SPLINE
      int flag=2;
#endif
#ifdef ALL_SPLINE
      int flag=3;
#endif

      double aa=ave_diff/(double)counter;
      printf("%d   %2d  %10d   %3d %3d %5.2f %5.2f   %20g  %20g DEBUG\n",
	     flag, ismooth, counter, CT_NBINS_D, CT_NBINS_XY, CT_EXPO, CT_SQUEEZE,
	     aa,sqrt(var_diff/(double)counter)-aa*aa);
    }
  /* else */
  /*   printf("N=0, no average or std dev\n"); */
  ave_diff=var_diff=0.0;
  counter=0;

#endif

#ifdef TRILINEAR
#ifdef HISTO
  // scrive l'istogramma della tabella

  // QUESTA PARTE ANDREBBE SCRITTA E COMMENTATA MEGLIO

  char fname[50];
  FILE *fd;
  int i,id,ix,iy;
  double x,y,l1,l2,l3;
  sprintf(fname,"popolazione.txt");
  if (!ismooth)
    fd=fopen(fname,"w");
  else
    fd=fopen(fname,"a");
  fprintf(fd,"################\n# smoothing %2d #\n################\n",ismooth);
  double ampl = sqrt(Smoothing.Variance[ismooth]);
  for (i=0; i<Ncomputations; i++)
    {
      id=i%CT_NBINS_D;
      ix=(i/CT_NBINS_D)%CT_NBINS_XY;
      iy=i/CT_NBINS_D/CT_NBINS_XY;

      x   = ix * bin_x;
      y   = iy * bin_x;
      l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
      l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
      l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

      fprintf(fd," %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20d\n",
	      ismooth,id,ix,iy,delta_vector[id],x,y,l1,l2,l3,CT_histo[i]);
    }
  fclose(fd);
#endif
#endif

  return 0;
}


int compare_search(const void *A, const void *B)
{
  if ( *(double*)A >= *(double*)B && *(double*)A < *((double*)B+1))
    return 0;
  if ( *(double*)A < *(double*)B )
    return -1;
  return 1;
}



double interpolate_collapse_time(double l1, double l2, double l3) //inline
{
  int ix, iy;
  double ampl, d, x, y;

  ampl = sqrt(Smoothing.Variance[ismooth]);
  d = (l1+l2+l3)/ampl;
  x = (l1-l2)/ampl;
  y = (l2-l3)/ampl;
  ix = (int)(x/bin_x);
  iy = (int)(y/bin_x);

  /* in the unlikely event the point is beyond the limits */
  if (ix>=CT_NBINS_XY-1)
    ix=CT_NBINS_XY-2;
  if (ix<0)
    ix=0;
  if (iy>=CT_NBINS_XY-1)
    iy=CT_NBINS_XY-2;
  if (iy<0)
    iy=0;

#ifdef ALL_SPLINE

  int ixx, iyy, ixstart, iystart, N=4;
  double xls[4],yls[4],zls[16];
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  if (ix==0)
    ixstart=0;
  else if (ix>=CT_NBINS_XY-2)
    ixstart=CT_NBINS_XY-4;
  else
    ixstart=ix-1;

  if (iy==0)
    iystart=0;
  else if (iy>=CT_NBINS_XY-2)
    iystart=CT_NBINS_XY-4;
  else
    iystart=iy-1;

  /* little 4x4 array for interpolation */
  for (ixx=0; ixx<4; ixx++)
    {
      xls[ixx]=(ixx+ixstart)*bin_x;
      yls[ixx]=(ixx+iystart)*bin_x;
    }
  for (ixx=0; ixx<4; ixx++)
    for (iyy=0; iyy<4; iyy++)
      zls[ixx+iyy*4]=my_spline_eval(CT_Spline[ixx+ixstart][iyy+iystart], d, accel);

  /* initialize interpolation */
  gsl_spline2d *LittleSpline = gsl_spline2d_alloc(T, N, N);
  gsl_spline2d_init(LittleSpline, xls, yls, zls, N, N);

  double a=gsl_spline2d_eval(LittleSpline, x, y, xacc, yacc);
  gsl_spline2d_free(LittleSpline);

  return a;

#endif

#ifdef TRILINEAR

  int id;
  if (d<=delta_vector[0])
    id=0;
  else if (d>=delta_vector[CT_NBINS_D-1])
    id = CT_NBINS_D-2;
  else
    id = (double*)bsearch((void*)&d, (void*)delta_vector, (size_t)(CT_NBINS_D-1), sizeof(double), compare_search)-delta_vector;

#ifdef HISTO
  CT_histo[id + ix*CT_NBINS_D + iy*CT_NBINS_D*CT_NBINS_XY]++;
#endif

  double dd = (d-delta_vector[id])/(delta_vector[id+1]-delta_vector[id]);
  double dx = x/bin_x-ix;
  double dy = y/bin_x-iy;

  return ((1.-dd)*(1.-dx)*(1.-dy)*CT_table[id   + (ix  )*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(1.-dx)*(1.-dy)*CT_table[id+1 + (ix  )*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(   dx)*(1.-dy)*CT_table[id   + (ix+1)*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(   dx)*(1.-dy)*CT_table[id+1 + (ix+1)*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(1.-dx)*(   dy)*CT_table[id   + (ix  )*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(1.-dx)*(   dy)*CT_table[id+1 + (ix  )*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(   dx)*(   dy)*CT_table[id   + (ix+1)*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(   dx)*(   dy)*CT_table[id+1 + (ix+1)*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]);

#endif

#ifdef BILINEAR_SPLINE

  double dx = x/bin_x-ix;
  double dy = y/bin_x-iy;

  return ((1.-dx)*(1.-dy)*my_spline_eval(CT_Spline[ix  ][iy  ], d, accel)+
	  (   dx)*(1.-dy)*my_spline_eval(CT_Spline[ix+1][iy  ], d, accel)+
	  (1.-dx)*(   dy)*my_spline_eval(CT_Spline[ix  ][iy+1], d, accel)+
	  (   dx)*(   dy)*my_spline_eval(CT_Spline[ix+1][iy+1], d, accel));

#endif

}

int check_CTtable_header()
{

  /* checks that the information contained in the header is consistent with the run 
     NB: this is done by Task 0 */

  int fail=0, dummy;
  double fdummy;

  fread(&dummy,sizeof(int),1,CTtableFilePointer);
#ifdef ELL_CLASSIC 
  if (dummy!=1)   /* classic ellipsoidal collapse, tabulated */
    {
      fail=1;
      printf("ERROR: CT table not constructed for ELL_CLASSIC, %d\n",
	     dummy);
    }
#endif
#ifdef ELL_SNG
#ifndef MOD_GRAV_FR
  if (dummy!=3)   /* standard gravity, numerical ellipsoidal collapse */
    {
      fail=1;
      printf("ERROR: CT table not constructed for ELL_SNG and standard gravity, %d\n",
	     dummy);
    }
#else
  if (dummy!=4)   /* f(R) gravity, numerical ellipsoidal collapse */
    {
      fail=1;
      printf("ERROR: CT table not constructed for ELL_SNG and MOD_GRAV_FR, %d\n",
	     dummy);
    }
#endif
#endif
  fread(&fdummy,sizeof(double),1,CTtableFilePointer);
  if (fabs(fdummy-params.Omega0)>1.e-10)
    {
      fail=1;
      printf("ERROR: CT table constructed for the wrong Omega0, %f in place of %f\n",
	     fdummy,params.Omega0);
    }
  fread(&fdummy,sizeof(double),1,CTtableFilePointer);
  if (fabs(fdummy-params.OmegaLambda)>1.e-10)
    {
      fail=1;
      printf("ERROR: CT table constructed for the wrong OmegaLambda, %f in place of %f\n",
	     fdummy,params.OmegaLambda);
    }
  fread(&fdummy,sizeof(double),1,CTtableFilePointer);
  if (fabs(fdummy-params.Hubble100)>1.e-10)
    {
      fail=1;
      printf("ERROR: CT table constructed for the wrong Hubble100, %f in place of %f\n",
	     fdummy,params.Hubble100);
    }

  fread(&dummy,sizeof(int),1,CTtableFilePointer);
  if (dummy != Ncomputations)
    {
      fail=1;
      printf("ERROR: CT table has the wrong size, %d in place of %d\n",
	     dummy,Ncomputations);
    }
  fread(&dummy,sizeof(int),1,CTtableFilePointer);
  if (dummy != CT_NBINS_D)
    {
      fail=1;
      printf("ERROR: CT table has the wrong density sampling, %d in place of %d\n",
	     dummy,CT_NBINS_D);
    }
  fread(&dummy,sizeof(int),1,CTtableFilePointer);
  if (dummy != CT_NBINS_XY)
    {
      fail=1;
      printf("ERROR: CT table has the wrong x and y sampling, %d in place of %d\n",
	     dummy,CT_NBINS_XY);
    }

  return fail;
}


void write_CTtable_header()
{

  int dummy;

#ifdef ELL_CLASSIC 
  dummy=1;   /* classic ellipsoidal collapse, tabulated */
  fwrite(&dummy,sizeof(int),1,CTtableFilePointer);
#endif
#ifdef ELL_SNG
#ifndef MOD_GRAV_FR
  dummy=3;   /* standard gravity, numerical ellipsoidal collapse */
  fwrite(&dummy,sizeof(int),1,CTtableFilePointer);
#else
  dummy=4;   /* f(R) gravity, numerical ellipsoidal collapse */
  fwrite(&dummy,sizeof(int),1,CTtableFilePointer);
#endif
#endif
  fwrite(&params.Omega0,sizeof(double),1,CTtableFilePointer);
  fwrite(&params.OmegaLambda,sizeof(double),1,CTtableFilePointer);
  fwrite(&params.Hubble100,sizeof(double),1,CTtableFilePointer);

  fwrite(&Ncomputations,sizeof(int),1,CTtableFilePointer);
  dummy=CT_NBINS_D;
  fwrite(&dummy,sizeof(int),1,CTtableFilePointer);
  dummy=CT_NBINS_XY;
  fwrite(&dummy,sizeof(int),1,CTtableFilePointer);

}

#endif /* TABULATED_CT */


#define SMALL ((double)1.e-20)

double ell(int ismooth, double l1, double l2, double l3) //inline
{
#ifdef ELL_CLASSIC
  double bc=ell_classic(ismooth, l1, l2, l3);
  if (bc>0.0)
    return 1.+InverseGrowingMode(bc,ismooth);
  else
    return 0.0;
#endif

#ifdef ELL_SNG
  double bc=ell_sng(ismooth, l1, l2, l3);
  if (bc>0.0)
    return 1./bc;
  else
    return 0.0;
#endif

}



/* implementation of ellipsoidal collapse following Nadkarni-Ghosh & Singhal (2016) */

#ifdef MOD_GRAV_FR

double ForceModification(double size, double a, double delta)
{

  double ff = 4.*params.OmegaLambda/params.Omega0;
  double thickness = FR0 / params.Omega0 / pow(H_over_c * size,2.0)
  * pow(a,7.)  * pow((1.+delta),-1./3.)
  * ( pow((1.0 + ff) / (1.0 + ff*pow(a,3.)) ,2.0) - pow((1.0 + ff) / (1.0 + delta + ff*pow(a,3.)) ,2.0));
  double F3=(thickness * (3. + thickness * (-3. + thickness)));
  if (F3<0.)
    F3=0.;

  return (F3<1. ? F3/3. : 1./3);
}

#endif

/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */

int sng_system(double t, const double y[], double f[], void *sng_par)
{

  int i,j;
  double sum;
  double omegam=OmegaMatter(1./t-1.);
  double omegal=OmegaLambda(1./t-1.);
  double delta=y[6]+y[7]+y[8];

  for (i=0; i<3; i++) 
    {
      sum=0.;
      for (j=0; j<3; j++) 
	{
	  if (i==j || y[i] == y[j])
	    continue;
	  else
	    sum += (y[j+6] - y[i+6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i+3]) 
					- (1. - y[j]) * (1. - y[j]) * (1. + y[j+3])) 
	      / ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
	}
      f[i]   = (y[i+3]*(y[i]-1.0)) / t;
      f[i+3] = (0.5*(y[i+3]*(omegam - 2.0*omegal - 2.0) 
#ifdef MOD_GRAV_FR
		     - 3.0*omegam*y[i+6] * (1. + ForceModification(*(double*)sng_par, t, delta))
#else
		     - 3.0*omegam*y[i+6] 
#endif
		     - 2.0*y[i+3]*y[i+3])) / t;
      f[i+6] = ((5./6. + y[i+6]) * 
		((3. + y[3] + y[4] + y[5]) 
		 - (1. + delta) / (2.5 + delta) * (y[3] + y[4] + y[5]))
		- (2.5 + delta) * (1. + y[i+3]) + sum) / t;

    }
  return GSL_SUCCESS;
}


double  ell_sng(int ismooth, double l1, double l2, double l3) //inline
{
  /* call of ODE solution for SNG ellipsoidal collapse */
  double ode_param, hh=1.e-6;
  double amin=1.e-5, amax=5.0;
  double mya = amin;

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9);
  gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);
  gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);
  gsl_odeiv2_system   ode_sys   = {sng_system, jac, 9, (void*)&ode_param};

  /* ICs */
#ifdef SCALE_DEPENDENT
  double D_in = GrowingMode(1./amin-1.,params.k_for_GM/Smoothing.Radius[ismooth]*params.InterPartDist);
#else
  double D_in = GrowingMode(1./amin-1.,0.);
#endif
  double y[9] = {l1*D_in, l2*D_in, l3*D_in,
		 l1*D_in/(l1*D_in - 1.), l2*D_in/(l2*D_in - 1.), l3*D_in/(l3*D_in - 1.),
		 l1*D_in, l2*D_in, l3*D_in};

#ifdef MOD_GRAV_FR
  if (ismooth<Smoothing.Nsmooth-1)
    ode_param=Smoothing.Radius[ismooth];
  else
    ode_param=Smoothing.Radius[ismooth-1];
#endif

  /* Loop to call ODE integrator */
  double olda=mya;
  double oldlam=y[0];

  while (mya<amax)
    {
      int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &mya, amax, &hh, y);
      if (status!=GSL_SUCCESS)
	{
	  printf("ERROR on task %d: integration of cosmological quantities failed\n",ThisTask);
	  fflush(stdout);
	  return -1;
	}
      if (y[0]>=0.99999)
	return olda + (1.-oldlam) * (mya - olda) / (y[0] - oldlam);
    }

  /* in this case the ellipsoid does not collapse */
  return 0;
}


double  ell_classic(int ismooth, double l1, double l2, double l3) //inline
{

  /*
    This rouine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
  */


  double ell;

  double del=l1+l2+l3;
  double det=l1*l2*l3;

  if (fabs(l1) < SMALL)
    ell = -0.1;

  else
    {
      double den = det / 126. + 5.*l1*del*(del-l1) / 84.;      

      if (fabs(den) < SMALL)
	{
	  if (fabs(del-l1) < SMALL)
	    {
	      if ( l1 > 0.0)
		ell = 1./l1;
	      else
		ell = -.1;
	    }
	  else
	    {
	      double dis = 7.*l1* (l1+6.*del);
	      if (dis < 0.0)
		ell = -.1;
	      else
		{
		  ell = (7.*l1 - sqrt(dis)) / (3.*l1*(l1-del));
		  if (ell < 0.) 
		    ell = -.1;
		}
	    }
	}
      else
	{
	  double rden = 1.0 / den;
	  
	  double a1   = 3.*l1*(del-l1) /14. * rden;
	  double a1_2 = a1*a1;
	  
	  double a2 = l1 * rden;
	  double a3 = -1.0 * rden;
	  double q  = (a1_2 - 3.*a2) / 9.;
	  double r  = (2.* a1_2*a1 - 9.*a1*a2 + 27.*a3) / 54.;

	  double r_2_q_3 = r*r-q*q*q;
	  
	  if (r_2_q_3 > 0)
	    {
	      double fabs_r = fabs(r);
	      double sq  = pow(sqrt(r_2_q_3)+fabs_r,0.333333333333333);
	      ell = -fabs_r / r * (sq+q/sq) - a1/3.;
	      if (ell < 0.) 
		ell = -.1;
	    }
	  else
	    {
	      double sq    = 2*sqrt(q);
	      double inv_3 = 1.0/3;
	      double t     = acos(2*r/q/sq);
	      double s1 = -sq*cos (t * inv_3)-a1 * inv_3;
	      double s2 = -sq*cos ((t+2.*PI) * inv_3) - a1 * inv_3;
	      double s3 = -sq*cos ((t+4.*PI) * inv_3) - a1 * inv_3;
	      
	      if (s1 < 0.) 
		s1 =1.e10;
	       
	      if (s2 < 0.) 
		s2 = 1.e10;
	      
	      if (s3 < 0.) 
		s3 = 1.e10;
	      
	      ell = (s1 < s2? s1 : s2);
	      ell = (s3 <ell? s3 : ell);
	      if (ell == 1.e10) 
		ell =- .1;
	    }
	}
    }
  
  if (del > 0. && ell > 0.)
    {
      double inv_del = 1.0/del;
      ell+=-.364 * inv_del * exp(-6.5*(l1-l2)*inv_del-2.8*(l2-l3)*inv_del);
    }

  return ell;
}


#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

void ord(double * restrict a, double * restrict b, double * restrict c)
{
  /* orders a,b,c in decreasing order a>b>c */
  double lo,hi;

  hi = max( max(*a,*b), *c );
  lo = min( min(*a,*b), *c );

  *b = *a + *b+ *c - lo - hi;
  *a = hi;
  *c = lo;
}
