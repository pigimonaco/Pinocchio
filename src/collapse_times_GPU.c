/*****************************************************************
 *                        PINOCCHI0  V5.0                        *
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

#if defined(GPU_OMP)

#include "pinocchio.h"
#include <assert.h>

/*------------------------------------------------------- Macros declaration --------------------------------------------------------*/
#define SMALL ((double)1.e-20)
#define INV_3 (1.0 / 3.0)
#define _MIN_(a,b) ((a) < (b) ? (a) : (b))
#define _MAX_(a,b) ((a) > (b) ? (a) : (b))
/*-----------------------------------------------------------------------------------------------------------------------------------*/

/* Orders a,b,c in decreasing order a>b>c */
void ord_gpu(double *const restrict a,
	     double *const restrict b,
	     double *const restrict c)
{
  const double hi = _MAX_(_MAX_(*a, *b), *c);
  const double lo = _MIN_(_MIN_(*a, *b), *c);
  *b = *a + *b + *c - lo - hi;
  *a = hi;
  *c = lo;

  return;
}

/*------------------------------------------------------- Functions implementation --------------------------------------------------------*/

/* Classical ellipsoidal collapse solution */
double ell_classic_gpu(const int    ismooth,
		       const double l1,
		       const double l2,
		       const double l3)
{
  /* The actual implementation solves the branch thread-divergence */
  
  /* Local variables declaration */
  const double del = (l1 + l2 + l3);
  const double det = (l1 * l2 * l3);

  double ell = 0.0;

  /* Vanishing lambda1 eigenvalue case */
  const unsigned int mask_l1 = ((l1 > -SMALL) && (l1 < SMALL));
  ell += (mask_l1 * -0.1);

  const double den = det / 126. + 5. * l1 * del * (del - l1) / 84.;

  const unsigned int mask_den = ((den > -SMALL) && (den < SMALL));
  /* Check 1st perturbative order conditions */

  const unsigned int mask_del_l1 = (((del - l1) > -SMALL) && ((del - l1) < SMALL));
  ell += (mask_del_l1 * mask_den * !mask_l1) * ((l1 > 0.0) ? (1.0 / l1) : -0.1); /* Zel'dovich approximation */
  /* Check 2nd perturbative order conditions */
  const double dis = (7.0 * l1 * (l1 + 6.0 * del));
  ell += (!mask_del_l1 * mask_den * !mask_l1) * ((dis < 0.0) ? -0.1 : (7. * l1 - sqrt(dis)) / (3. * l1 * (l1 - del)));
  /* 3rd order perturbative solution. For more details about the equations implemented, see Monaco 1996a */
            
  /* Intermediate values */
  const double rden = (mask_den ? 0.0 : (1.0 / den));
  const double a1 = 3. * l1 * (del - l1) / 14. * rden;
  const double a1_2 = a1 * a1;
  const double a2 = l1 * rden;
  const double a3 = -1.0 * rden;

  /* The collapse time b_c will be a combination of R, Q, and D == den */
  const double q = (a1_2 - 3. * a2) / 9.;
  const double r = (2. * a1_2 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
  const double r_2_q_3 = r * r - q * q * q;

  /* Check 3rd perturbative order conditions */

  /* ---------------- Case 1 --------------- */
  /* If R^2 - Q^2 > 0, which is valid for spherical and quasi-spherical perturbations */

  const unsigned int mask_r_2_q_3 = (r_2_q_3 > 0.0);
  /* 3rd order solution */
  const double fabs_r = ((r > 0.0) ? r : -r);
  const double inv_r = (((r > -SMALL) && (r < SMALL)) ? 0.0 : (1.0 / r));
  const double sq = pow(sqrt(mask_r_2_q_3 * r_2_q_3) + fabs_r, 0.333333333333333);
  const double inv_sq = (((sq > -SMALL) && (sq < SMALL)) ? 0.0 : (1.0 / sq));
  ell += (mask_r_2_q_3 * !mask_den * !mask_l1) * ((-fabs_r * inv_r) * (sq + (q * inv_sq)) - (a1 * INV_3));

  /* ---------------- Case 2 --------------- */
  /* The solution has to be chosen as the smallest non-negative one between s1, s2, and s3 */            
  const double sq_ = !mask_r_2_q_3 * (2.0 * sqrt((q > 0.0) * q));
  const double inv_q = (((q > -SMALL) && (q < SMALL)) ? 0.0 : (1.0 / q));
  const double inv_sq_ = (((sq_ > -SMALL) && (sq_ < SMALL)) ? 0.0 : (1.0 / sq_));
  const double t = !mask_r_2_q_3 * acos(2.0 * r * inv_q * inv_sq_);
  const double a1_inv_3 = !mask_r_2_q_3 * (a1 * INV_3);

  double s1 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos(t * INV_3) - a1_inv_3);
  double s2 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos((t + 2. * PI) * INV_3) - a1_inv_3);
  double s3 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos((t + 4. * PI) * INV_3) - a1_inv_3);

  ord_gpu(&s1, &s2, &s3);

  ell += (s3 > 0.0) * s3;
  ell += ((s3 < 0.0) && (s2 > 0.0)) * s2;
  ell += ((s3 < 0.0) && (s2 < 0.0) && (s1 > 0.0)) * s1;      

  const unsigned int mask_del_ell = ((del > 0.0) && (ell > 0.0));
  const double inv_del = (mask_del_ell ? (1.0 / del) : 0.0);
  ell += mask_del_ell * (-0.364 * inv_del * exp(-6.5 * (l1 - l2) * inv_del - 2.8 * (l2 - l3) * inv_del));  
  
  return ell;
}

/*---------------------------------------- Calculation of b_c == growing mode at collapse time ---------------------------------------*/

double ell_gpu(const int ismooth,
	       const double l1,
	       const double l2,
	       const double l3)
{

#ifdef ELL_CLASSIC

    const double bc = ell_classic_gpu(ismooth, l1, l2, l3);

    return ((bc > 0.0) * (1.0 + InverseGrowingMode(bc, ismooth)));

#else

#error "GPU works with ELL_CLASSIC only!"    

#endif // ELL_CLASSIC
}

/* ------------------------------------  Computation of collapse time i.e. F = 1 + z_collapse and variance ------------------------------------ */

double inverse_collapse_time_gpu(const int                     ismooth,
				 const double * const restrict dtensor,
				       int    * const restrict fail)
{  
  /* Local variables declaration */
  /* mu1, mu2 and mu3 are the principal invariants of the 3x3 tensor of second derivatives */
  /*  Performs diagonalization of the tensor by calculating the values of mu1, mu2, and mu3 */
  const double add[6] = {dtensor[0] * dtensor[0],
			 dtensor[1] * dtensor[1],
			 dtensor[2] * dtensor[2],
			 dtensor[3] * dtensor[3],
			 dtensor[4] * dtensor[4],
			 dtensor[5] * dtensor[5]};
  
  const double mu1   = (dtensor[0] + dtensor[1] + dtensor[2]);
  const double mu1_2 = (mu1 * mu1);
  const double mu2   = ((0.5 * mu1_2) - (0.5 * (add[0] + add[1] + add[2])) - (add[3] + add[4] + add[5]));
  
  /* mu3 calculation */
  const double mu3 = ((dtensor[0] * dtensor[1] * dtensor[2])       +
		      (2.0 * dtensor[3] * dtensor[4] * dtensor[5]) -
		      (dtensor[0] * add[5])                        -
		      (dtensor[1] * add[4])                        -
		      (dtensor[2] * add[3]));

  /* Check if the tensor is already diagonal */
  const double q = (mu1_2 - 3.0 * mu2) / 9.0;

  // q == 0.0
  const unsigned int mask_q0 = ((q <= EPSILON) && (q >= -EPSILON));
  const double x_q0[3] = {dtensor[0], dtensor[1], dtensor[2]};
  const double inv_q = (mask_q0 ? 0.0 : (1.0 / q));
  
  // q > 0.0
  const double r  = -(((2.0 * mu1_2 * mu1) - (9.0 * mu1 * mu2) + (27.0 * mu3)) / 54.0);
  const unsigned int mask = (((q * q * q) < (r * r)) || (q < 0.0));
  *fail = (mask ? 1 : 0);
  if (mask) // kernel abort
    return -10.0;
    
  const double sq = (2.0 * sqrt(q));
  const double inv_sq = (mask_q0 ? 0.0 : (1.0 / sq));
  const double t = acos(2.0 * r * inv_q * inv_sq);
  const double x_q_gt_0[3] = {((-sq * cos(t * INV_3)) + (mu1 * INV_3)),
			      ((-sq * cos((t + 2.0 * PI) * INV_3)) + (mu1 * INV_3)),
			      ((-sq * cos((t + 4.0 * PI) * INV_3)) + (mu1 * INV_3))};

  /* Ordering and inverse collapse time */
  double x1 = (mask_q0 * x_q0[0]) + (!mask_q0 * x_q_gt_0[0]);
  double x2 = (mask_q0 * x_q0[1]) + (!mask_q0 * x_q_gt_0[1]);
  double x3 = (mask_q0 * x_q0[2]) + (!mask_q0 * x_q_gt_0[2]);

  /* Ordering and inverse collapse time */
  ord_gpu(&x1, &x2, &x3);

  const double ret = ell_gpu(ismooth, x1, x2, x3);

  return ret;
}

/* Function: common_initialization */
/* Performed once by the host */
void common_initialization_gpu(const unsigned int size)
{
  /* init on the GPU gpu_products */
  /* synch kernel */

#if defined(GPU_OMP_DEBUG)
 #pragma omp target device(internal.device.devID)
#else
 #pragma omp target teams distribute parallel for device(internal.device.devID)
#endif // GPU_OMP_DEBUG
  for (unsigned int i=0 ; i<size ; i++)
    {
      gpu_products.Rmax[i] = -1;
      gpu_products.Fmax[i] = (PRODFLOAT)-10.0;
    }
  
#pragma omp parallel for
  for (unsigned int i=0 ; i<size ; i++)
    {
      /*----------- Common initialization ----------- */
      products[i].Rmax = -1;
      products[i].Fmax = (PRODFLOAT)-10.0;
            
      /*----------- Zel'dovich case ----------- */
      products[i].Vel[_x_] = (PRODFLOAT)0.0;
      products[i].Vel[_y_] = (PRODFLOAT)0.0;
      products[i].Vel[_z_] = (PRODFLOAT)0.0;

#ifdef TWO_LPT

      products[i].Vel_2LPT[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_2LPT[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_2LPT[_z_] = (PRODFLOAT)0.0;

#ifdef THREE_LPT

      products[i].Vel_3LPT_1[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_1[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_1[_z_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_x_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_y_] = (PRODFLOAT)0.0;
      products[i].Vel_3LPT_2[_z_] = (PRODFLOAT)0.0;

#endif // THREE_LPT
#endif // TWO_LPT
    } // loop over MyGrids[0].total_local_size
  
  return;
}

int compute_collapse_times_gpu(int ismooth)
{
  /*----------------------------------------------------------------*/
  // Size of the products array
  const unsigned int total_size = MyGrids[0].total_local_size;
  
  if (!ismooth)
    {
      /*----------- Common initialization ----------- */
      common_initialization_gpu(total_size);
    }

  /* timing the main loop of 'compute_collapse_times' function */
  double cputmp = MPI_Wtime();  
  /*--------------------- GPU memory movements ----------------------------------------*/

  double tmp = MPI_Wtime();
  
  /* copy second_derivatives to the GPU using multiple threads */
 #pragma omp parallel
  {
    const int tid = omp_get_thread_num();
    const int thr = omp_get_num_threads();
    
    for (int i=tid ; i<6 ; i+=thr)
      {
	omp_target_memcpy((void *)internal.device.gpu_main_memory, 
			  (void *)second_derivatives[0][i],
			  (total_size * sizeof(double)),
			  (internal.device.memory_second_derivatives.offset + internal.device.memory_second_derivatives.tensor[i]),
			  0,
			  internal.device.devID,
			  internal.device.hostID);
      }
  } // omp parallel
  
  cputime.gpu_mem_transf += (MPI_Wtime() - tmp);
  
  /*-----------------------------------------------------------------------------------*/
  
  /*--------------------- GPU CALCULATION OF COLLAPSE TIME ----------------------------*/
  /*----------------- Calculation of variance, average, and collapse time -------------*/
  
  /* Local average and variance declaration */
  double local_average = 0.0, local_variance = 0.0;
  int all_fails = 0; 
  
  tmp = MPI_Wtime();

#if defined(GPU_OMP_DEBUG)
   #pragma omp target map(tofrom: local_average, local_variance, all_fails) device(internal.device.devID)
#else
   #pragma omp target teams distribute parallel for reduction(+: local_average, local_variance, all_fails) device(internal.device.devID)
#endif // GPU_OMP_DEBUG  
  for (unsigned int index=0 ; index<total_size ; index++)
    {
      /* Computation of second derivatives of the potential i.e. the gravity Hessian */
      double diff_ten[6]; 
      for (int i=0 ; i<6 ; i++)
	{
	  diff_ten[i] = gpu_second_derivatives.tensor[i][index];
	}
        
      /* Computation of the variance of the linear density field */
      const double delta = (diff_ten[0] + diff_ten[1] + diff_ten[2]);
      local_average  += delta;
      local_variance += (delta * delta);

      /* Computation of the collapse time */
      int    fail;	
      /* inverse_collapse_time(funzione ell) ---- qui si usa le GPU spline */
      const double Fnew = inverse_collapse_time_gpu(ismooth, diff_ten, &fail);
      all_fails += fail;
      
      /* Updating collapse time */
      const int update = (gpu_products.Fmax[index] < Fnew);
      gpu_products.Rmax[index] = (update ? ismooth : gpu_products.Rmax[index]);
      gpu_products.Fmax[index] = (update ? Fnew    : gpu_products.Fmax[index]);
      
    } // target region

  /*-------------------- END OF GPU CALCULATION -----------------------------------------------*/
  cputime.gpu_computation += (MPI_Wtime() - tmp);
    
  /* Fail check during computation of the inverse collapse time */
  /* If there were failures, an error message is printed and the function returns 1 */
  if (all_fails)
    {
      printf("ERROR on task %d: failure in inverse_collapse_time\n", ThisTask);
      fflush(stdout);
      return 1;
    }
  
  /*--------------------- GPU memory movements ----------------------------------------*/

  tmp = MPI_Wtime();

  #pragma omp parallel
  {
    #pragma omp single nowait
    {
      /* copy products.Rmax from the GPU */
      omp_target_memcpy((void *)host_products.Rmax,
			(void *)internal.device.gpu_main_memory,
			(total_size * sizeof(int)),
			0,
			(internal.device.memory_products.offset + internal.device.memory_products.Rmax),
			internal.device.hostID,
			internal.device.devID);
    } // omp single nowait

    #pragma omp single nowait
    {
      /* copy products.Fmax from the GPU */
      omp_target_memcpy((void *)host_products.Fmax,
			(void *)internal.device.gpu_main_memory,
			(total_size * sizeof(PRODFLOAT)),
			0,
			(internal.device.memory_products.offset + internal.device.memory_products.Fmax),
			internal.device.hostID,
			internal.device.devID);
    } // omp single nowait

    #pragma omp barrier

    #pragma omp single nowait
    {
      cputime.gpu_mem_transf += (MPI_Wtime() - tmp);
    }
  
    /* update products */
    #pragma omp for nowait
    for (unsigned int index=0 ; index<total_size ; index++)
      {
	products[index].Rmax = host_products.Rmax[index];
	products[index].Fmax = host_products.Fmax[index];
      }
  } // omp parallel
  
  /* ------------- Updates cpu collapse time for a single threads -----------------------------*/
  
  /* Calculating obtained variance and avarage */		
  double global_variance = 0.0;
  double global_average  = 0.0;
    
  /* Calculates the global variance and average by reducing the values of local_variance and local_average across all MPI tasks */
  MPI_Allreduce(&local_variance, &global_variance, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&local_average , &global_average , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  global_variance /= (double)MyGrids[0].Ntotal;
  global_average  /= (double)MyGrids[0].Ntotal;
	
  /* Stores the true variance*/
  Smoothing.TrueVariance[ismooth] = global_variance;

  /* CPU collapse time */	
  cputime.coll += (MPI_Wtime() - cputmp);
  
  return 0;

} // compute_collapse_times

#endif // GPU_OMP
