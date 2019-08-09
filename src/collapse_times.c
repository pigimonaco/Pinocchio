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

#include "pinocchio.h"

int    inline v_inverse_collapse_time ( dvec dtensor[6], dvec_u , dvec , dvec_u *, dvec_u *, dvec_u *, dvec_u *) __attribute__((always_inline));
double inline inverse_collapse_time   ( double *, double *, double *, double *, int *) __attribute__((always_inline));
double inline ell                     ( double, double, double, double, double) __attribute__((always_inline));

void ord(double *,double *,double *);

#if defined( _OPENMP )

double cputime_ell;
#pragma omp threadprivate(cputime_ell)
int fails;
#pragma omp threadprivate(fails)
#else

#define cputime_ell       cputime.ell

#endif

int compute_collapse_times(int ismooth)
{

  double local_average,local_variance;
  dvec_u  vlocal_average, vlocal_variance;
  dvec    vzero = {0, 0, 0, 0};
  
  /* initialize the variance of the slice */
  local_variance = 0.0;
  local_average  = 0.0;
  vlocal_average.V  = vzero;
  vlocal_variance.V = vzero; 

  /* initialize Fmax if it is the first smoothing */
  if (!ismooth)
#pragma omp parallel for
    for (int i = 0; i < MyGrids[0].total_local_size; i++)
      {
	products[i].Fmax   = -10.0;
	products[i].Rmax   = -1;
	products[i].Vmax[0] = 0.0;
	products[i].Vmax[1] = 0.0;
	products[i].Vmax[2] = 0.0;
#ifdef TWO_LPT
	products[i].Vmax_2LPT[0] = 0.0;
	products[i].Vmax_2LPT[1] = 0.0;
	products[i].Vmax_2LPT[2] = 0.0;
#ifdef THREE_LPT
	products[i].Vmax_3LPT_1[0] = 0.0;
	products[i].Vmax_3LPT_1[1] = 0.0;
	products[i].Vmax_3LPT_1[2] = 0.0;
	products[i].Vmax_3LPT_2[0] = 0.0;
	products[i].Vmax_3LPT_2[1] = 0.0;
	products[i].Vmax_3LPT_2[2] = 0.0;
#endif
#endif
      }

  /* loop on all particles in the slice */

#if defined( _OPENMP )
  double invcoll_update = 0, ell_update = 0;
#endif
#pragma omp parallel
  {
#if defined( _OPENMP )
    double  mylocal_average, mylocal_variance;
    double  cputime_invcoll;
    dvec_u  vmylocal_average, vmylocal_variance;
    dvec    vzero = {0, 0, 0, 0};
    int     fails = 0;
    int     pid = omp_get_thread_num();
    
    /* initialize the variance of the slice */
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
    
#pragma omp for nowait
    for (int local_x = 0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
      {
	int idx_x = local_x * MyGrids[0].GSlocal[_y_];
	for (int local_y = 0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	  {
	    int idx_y = (idx_x + local_y) * MyGrids[0].GSlocal[_z_];
	    int mysize = DVEC_SIZE * MyGrids[0].GSlocal[_z_] / DVEC_SIZE;
	    int local_z = 0;
	    for ( ; local_z < mysize; local_z += DVEC_SIZE)
	      {
		int index = idx_y + local_z;
		/* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */

		/* sum all contributions */
		dvec elements[6] __attribute__((aligned(32)));
	      
#if !(defined(__aarch64__) || defined(__arm__))
		for (int i = 0; i < 6; i++)
		  elements[i] = _mm256_load_pd(&second_derivatives[0][i][index]);
#else
		dvec_u ld_elements[6] __attribute__((aligned(32)));
		for (int i = 0; i < 6; i++)
		  {
		    ld_elements[i].v[0] = second_derivatives[0][i][index ];
		    ld_elements[i].v[1] = second_derivatives[0][i][index + 1];
		    ld_elements[i].v[2] = second_derivatives[0][i][index + 2];
		    ld_elements[i].v[3] = second_derivatives[0][i][index + 3];
		    elements[i] = ld_elements[i].V;
		  }
#endif

		/* Computation of the variance of the linear density field */
		dvec_u delta __attribute__((aligned(32)));
		delta.V = elements[0] + elements[1] + elements[2];
		dvec delta_2 __attribute__((aligned(32)));
		delta_2 = delta.V*delta.V;
	      
		vmylocal_average.V  += delta.V;
		vmylocal_variance.V += delta_2;

		/* Computation of ellipsoidal collapse */
		dvec_u lambda1,lambda2,lambda3 __attribute__((aligned(32)));
		dvec_u Fnew __attribute__((aligned(32)));
		double tmp = MPI_Wtime();
		if(v_inverse_collapse_time(elements, delta, delta_2, &lambda1, &lambda2, &lambda3, &Fnew))
		  {
		    printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
		    fflush(stdout);
#if !defined( _OPENMP )		    
		    return 1;
#else
		    fails = 1;
#endif
		  }
		cputime_invcoll += MPI_Wtime() - tmp;

#ifdef SCALE_DEPENDENT_GROWTH
		SDGM.flag=0;
		SDGM.ismooth=ismooth;
#endif	      
		for(int ii = 0; ii < DVEC_SIZE; ii++)
		  {
		    /* comment these two lines out if F is the inverse growing mode */
		    if (Fnew.v[ii] > 0)
		      Fnew.v[ii] = 1. + InverseGrowingMode(1./Fnew.v[ii]);
		  
		    if (products[index + ii].Fmax < Fnew.v[ii])
		      {
			products[index].Fmax = Fnew.v[ii];
			products[index].Rmax = ismooth;
		      }
		  }
	      }
	    for ( ; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	      {
		int index = idx_y + local_z;
		/* sum all contributions */
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
		double tmp = MPI_Wtime();
		double Fnew = inverse_collapse_time(diff_ten, &lambda1, &lambda2, &lambda3, &fail);
		cputime_invcoll += MPI_Wtime() - tmp;
		
#ifdef SCALE_DEPENDENT_GROWTH
		SDGM.flag=0;
		SDGM.ismooth=ismooth;
#endif
		/* comment these two lines out if F is the inverse growing mode */
		if (Fnew>0)
		  Fnew = 1.+InverseGrowingMode(1./Fnew);

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
	  }
      }

    
    
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
        
  }

  int all_fails = 0;
  
  /* Computes the obtained variance */

#if defined (_OPENMP)
#pragma omp atomic
    cputime.invcoll += invcoll_update/num_omp_th;
#pragma omp atomic	
    cputime.ell     += ell_update/num_omp_th;
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
      double Np = (double)MyGrids[0].GSglobal[_x_] * (double)MyGrids[0].GSglobal[_y_] *(double)MyGrids[0].GSglobal[_z_];
      global_variance /= Np;
      global_average  /= Np;
    }

  MPI_Bcast(&global_variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  Smoothing.TrueVariance[ismooth] = global_variance;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=0;
  SDGM.ismooth=ismooth;
  Smoothing.TrueVariance[ismooth] *= GrowingMode(params.camb.ReferenceRedshift);
#endif


  return 0;
}



int inline v_inverse_collapse_time(dvec dtensor[6], dvec_u delta, dvec delta_2, dvec_u * restrict x1, dvec_u * restrict x2, dvec_u * restrict x3, dvec_u * restrict res)
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
#if !(defined(__aarch64__) || defined(__arm__) )  
  sq.V = 2.0 * _mm256_sqrt_pd(q.V);
#else
  sq.v[0] = 2.0 * sqrt(q.v[0]);
  sq.v[1] = 2.0 * sqrt(q.v[1]);
  sq.v[2] = 2.0 * sqrt(q.v[2]);
  sq.v[3] = 2.0 * sqrt(q.v[3]);
#endif
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
  /* ord(&x1,&x2,&x3); */
  /* t=ell(x1,x2,x3,mu1,mu3); */
  double tmp = MPI_Wtime();
  for(int ii = 0; ii < DVEC_SIZE; ii++)
    {
      double _x1 = (*x1).v[ii];
      double _x2 = (*x2).v[ii];
      double _x3 = (*x3).v[ii];
      ord(&_x1, &_x2, &_x3);
      double t = ell(_x1, _x2, _x3, delta.v[ii], mu3.v[ii]);
      if(t > 0)
	(*res).v[ii] = 1.0 / t;
      else
	(*res).v[ii] = -10.0;
    }
  cputime_ell += MPI_Wtime() - tmp;

  return 0;
}



double inline inverse_collapse_time(double * restrict deformation_tensor, double * restrict x1, double * restrict x2, double * restrict x3, int * restrict fail)
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
	  *fail=1;
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
  /* ord(&x1,&x2,&x3); */
  /* t=ell(x1,x2,x3,mu1,mu3); */
  ord(x1, x2, x3);
  double tmp = MPI_Wtime();
  double t = ell(*x1, *x2, *x3, mu1, mu3);
  cputime_ell += MPI_Wtime() - tmp;
  
  /* return value */
  if (t > 0.)
    return 1. / t;
  else
    return -10.;
}





double inline ell(double l1,double l2,double l3,double del,double det)
{

  /*
    This rouine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
  */


  double ell;


  if (l1 == 0.)
    ell = -0.1;
  else
    {
      double den = det / 126. + 5.*l1*del*(del-l1) / 84.;      

      if (den == 0.)
	{
	  if (del == l1)
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


int compute_velocities(int ismooth)
{

#ifdef SMOOTH_VELOCITIES
  int final=(ismooth==Smoothing.Nsmooth-1);
#endif
  /* loop on all particles in the slice */
  for (int local_x = 0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
    {
      int idx_x = local_x * MyGrids[0].GSlocal[_y_];
      for (int local_y=0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	{
	  int idx_y = (idx_x + local_y)*MyGrids[0].GSlocal[_z_];
	  for (int local_z=0; local_z < MyGrids[0].GSlocal[_x_]; local_z++)	    
	    {
	      int index = local_z + idx_y;

		/* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */
		
#ifdef SMOOTH_VELOCITIES
		if (products[index].Rmax==ismooth || (products[index].Fmax<0.0 && final) )
#endif
		  for (int i = 0; i < 3; i++)
		    products[index].Vmax[i]=first_derivatives[0][i][index];
	    }
	}
    }
  return 0;
}
