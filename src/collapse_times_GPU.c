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

#if defined(GPU_OMP) || defined(GPU_OMP_FULL)

#include "pinocchio.h"
#include <assert.h>

/*------------------------------------------------------- Macros declaration --------------------------------------------------------*/
#define SMALL ((double)1.e-20)
#define INV_3 (1.0 / 3.0)
#define _MIN_(a,b) ((a) < (b) ? (a) : (b))
#define _MAX_(a,b) ((a) > (b) ? (a) : (b))

#ifdef TABULATED_CT  // Tabulated collapse time calculation

#define CT_NBINS_XY (50) 
#define CT_NBINS_D (100) 
#define CT_SQUEEZE (1.2) 
#define CT_EXPO   (1.75)  
#define CT_RANGE_D (7.0)  
#define CT_RANGE_X (3.5)
#define CT_DELTA0 (-1.0)

#endif  // end TABULATED_CT
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
  ell                        += (mask_l1 * -0.1);

  const double den = det / 126. + 5. * l1 * del * (del - l1) / 84.;

  const unsigned int mask_den = ((den > -SMALL) && (den < SMALL));
  /* Check 1st perturbative order conditions */

  const unsigned int mask_del_l1 = (((del - l1) > -SMALL) && ((del - l1) < SMALL));
  ell                            += (mask_del_l1 * mask_den * !mask_l1) * ((l1 > 0.0) ? (1.0 / l1) : -0.1); /* Zel'dovich approximation */
  /* Check 2nd perturbative order conditions */
  const double dis = (7.0 * l1 * (l1 + 6.0 * del));
  ell              += (!mask_del_l1 * mask_den * !mask_l1) * ((dis < 0.0) ? -0.1 : (7. * l1 - sqrt(dis)) / (3. * l1 * (l1 - del)));
  /* 3rd order perturbative solution. For more details about the equations implemented, see Monaco 1996a */
            
  /* Intermediate values */
  const double rden = (mask_den ? 0.0 : (1.0 / den));
  const double a1   = 3. * l1 * (del - l1) / 14. * rden;
  const double a1_2 = a1 * a1;
  const double a2   = l1 * rden;
  const double a3   = -1.0 * rden;

  /* The collapse time b_c will be a combination of R, Q, and D == den */
  const double q       = (a1_2 - 3. * a2) / 9.;
  const double r       = (2. * a1_2 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
  const double r_2_q_3 = r * r - q * q * q;

  /* Check 3rd perturbative order conditions */

  /* ---------------- Case 1 --------------- */
  /* If R^2 - Q^2 > 0, which is valid for spherical and quasi-spherical perturbations */

  const unsigned int mask_r_2_q_3 = (r_2_q_3 > 0.0);
  /* 3rd order solution */
  const double fabs_r = ((r > 0.0) ? r : -r);
  const double inv_r  = (((r > -SMALL) && (r < SMALL)) ? 0.0 : (1.0 / r));
  const double sq     = pow(sqrt(mask_r_2_q_3 * r_2_q_3) + fabs_r, 0.333333333333333);
  const double inv_sq = (((sq > -SMALL) && (sq < SMALL)) ? 0.0 : (1.0 / sq));
  ell                 += (mask_r_2_q_3 * !mask_den * !mask_l1) * ((-fabs_r * inv_r) * (sq + (q * inv_sq)) - (a1 * INV_3));

  /* ---------------- Case 2 --------------- */
  /* The solution has to be chosen as the smallest non-negative one between s1, s2, and s3 */            
  const double sq_      = !mask_r_2_q_3 * (2.0 * sqrt((q > 0.0) * q));
  const double inv_q    = (((q > -SMALL) && (q < SMALL)) ? 0.0 : (1.0 / q));
  const double inv_sq_  = (((sq_ > -SMALL) && (sq_ < SMALL)) ? 0.0 : (1.0 / sq_));
  const double t        = !mask_r_2_q_3 * acos(2.0 * r * inv_q * inv_sq_);
  const double a1_inv_3 = !mask_r_2_q_3 * (a1 * INV_3);

  double s1 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos(t * INV_3) - a1_inv_3);
  double s2 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos((t + 2. * PI) * INV_3) - a1_inv_3);
  double s3 = (!mask_r_2_q_3 * !mask_den * !mask_l1) * (-sq_ * cos((t + 4. * PI) * INV_3) - a1_inv_3);

  ord_gpu(&s1, &s2, &s3);

  ell += (s3 > 0.0) * s3;
  ell += ((s3 < 0.0) && (s2 > 0.0)) * s2;
  ell += ((s3 < 0.0) && (s2 < 0.0) && (s1 > 0.0)) * s1;      

  const unsigned int mask_del_ell = ((del > 0.0) && (ell > 0.0));
  const double inv_del            = (mask_del_ell ? (1.0 / del) : 0.0);
  ell                             += mask_del_ell * (-0.364 * inv_del * exp(-6.5 * (l1 - l2) * inv_del - 2.8 * (l2 - l3) * inv_del));  
  
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
  const double x_q0[3]       = {dtensor[0], dtensor[1], dtensor[2]};
  const double inv_q         = (mask_q0 ? 0.0 : (1.0 / q));
  
  // q > 0.0
  const double r          = -(((2.0 * mu1_2 * mu1) - (9.0 * mu1 * mu2) + (27.0 * mu3)) / 54.0);
  const unsigned int mask = (((q * q * q) < (r * r)) || (q < 0.0));
  *fail                   = (mask ? 1 : 0);
  if (mask) // kernel abort
    return -10.0;
    
  const double sq          = (2.0 * sqrt(q));
  const double inv_sq      = (mask_q0 ? 0.0 : (1.0 / sq));
  const double t           = acos(2.0 * r * inv_q * inv_sq);
  const double x_q_gt_0[3] = {((-sq * cos(t * INV_3)) + (mu1 * INV_3)),
			                        ((-sq * cos((t + 2.0 * PI) * INV_3)) + (mu1 * INV_3)),
			                        ((-sq * cos((t + 4.0 * PI) * INV_3)) + (mu1 * INV_3))};

  /* Ordering and inverse collapse time */
  double x1 = (mask_q0 * x_q0[0]) + (!mask_q0 * x_q_gt_0[0]);
  double x2 = (mask_q0 * x_q0[1]) + (!mask_q0 * x_q_gt_0[1]);
  double x3 = (mask_q0 * x_q0[2]) + (!mask_q0 * x_q_gt_0[2]);

  /* Ordering and inverse collapse time */
  ord_gpu(&x1, &x2, &x3);

  #ifdef TABULATED_CT
  const double ret = interpolate_collapse_time(ismooth,x1,x2,x3);
  #else 
  const double ret = ell_gpu(ismooth, x1, x2, x3);
  #endif

  return ret;
}

/* Function: common_initialization */
/* Performed once by the host */
void common_initialization_gpu(const unsigned int size)
{
  /* init on the GPU gpu_products */
  /* synch kernel */
  
  #pragma omp target teams distribute parallel for device(devID)
  for (unsigned int i = 0; i < size; i++) {
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

  /* PMT measure */
  PMT_CPU_START("collapse_time_CPU", ThisTask);
  PMT_GPU_START("collapse_time_GPU", devID, ThisTask);
  /* timing the main loop of 'compute_collapse_times' function */
  double cputmp, tmp;  
  cputmp = tmp = MPI_Wtime();  
  
  /*--------------------- GPU memory movements ----------------------------------------*/

  #if defined(GPU_OMP) 
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
			  devID,
			  hostID);
      }
  } // omp parallel

  #elif defined(GPU_OMP_FULL) 

  // Update each individual each individual second_derivatives
  for (int igrid = 0; igrid < Ngrids; igrid++) {
    for (int i = 0; i < 6; i++) {
        #pragma omp target update to(second_derivatives[igrid][i][0:total_size]) device(devID)
    }
  }

  #endif

  gputime.memory_transfer.collapse_times += (MPI_Wtime() - tmp);
  
  /*-----------------------------------------------------------------------------------*/
  
  /*--------------------- GPU CALCULATION OF COLLAPSE TIME ----------------------------*/
  /*----------------- Calculation of variance, average, and collapse time -------------*/

  /* Local average and variance declaration */
  double local_average = 0.0, local_variance = 0.0;
  int all_fails        = 0; 
  
  /* const size_t Nblocks = ((total_size + GPU_OMP_BLOCK - 1) / GPU_OMP_BLOCK); */
  
  tmp = MPI_Wtime();

#if defined(GPU_OMP_DEBUG)
   #pragma omp target map(tofrom: local_average, local_variance, all_fails) device(devID)
#else
   #pragma omp target teams distribute parallel for reduction(+: local_average, local_variance, all_fails) device(devID)
#endif // GPU_OMP_DEBUG  
  for (unsigned int index=0 ; index<total_size ; index++)
    {
      /* Computation of second derivatives of the potential i.e. the gravity Hessian */
      double diff_ten[6]; 
      for (int i=0 ; i<6 ; i++)
	    { 
        #if !defined(GPU_OMP_FULL)
	      diff_ten[i] = gpu_second_derivatives.tensor[i][index];
        #else
        diff_ten[i] = second_derivatives[0][i][index];
        #endif
	    }
        
      /* Computation of the variance of the linear density field */
      const double delta = (diff_ten[0] + diff_ten[1] + diff_ten[2]);
      local_average      += delta;
      local_variance     += (delta * delta);

      /* Computation of the collapse time */
      int fail;	
      /* inverse_collapse_time(funzione ell) ---- qui si usa le GPU spline */
      const double Fnew = inverse_collapse_time_gpu(ismooth, diff_ten, &fail);
      all_fails         += fail;
      
      /* Updating collapse time */      
      const int update         = (gpu_products.Fmax[index] < Fnew);
      gpu_products.Rmax[index] = (update ? ismooth : gpu_products.Rmax[index]);
      gpu_products.Fmax[index] = (update ? Fnew    : gpu_products.Fmax[index]);
    } // target region

  /*-------------------- END OF GPU CALCULATION -----------------------------------------------*/
  gputime.computation.collapse_times += (MPI_Wtime() - tmp);
    
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
  
  #if defined(GPU_OMP)
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
			hostID,
			devID);
    } // omp single nowait

    #pragma omp single nowait
    {
      /* copy products.Fmax from the GPU */
      omp_target_memcpy((void *)host_products.Fmax,
			(void *)internal.device.gpu_main_memory,
			(total_size * sizeof(PRODFLOAT)),
			0,
			(internal.device.memory_products.offset + internal.device.memory_products.Fmax),
			hostID,
			devID);
    } // omp single nowait

    #pragma omp barrier

    #pragma omp master
    {
      gputime.memory_transfer.collapse_times += (MPI_Wtime() - tmp);
    }
  
    /* update products */
    #pragma omp for nowait
    for (unsigned int index=0 ; index<total_size ; index++)
      {
	      products[index].Rmax = host_products.Rmax[index];
	      products[index].Fmax = host_products.Fmax[index];
      }
  } // omp parallel

  #elif defined(GPU_OMP_FULL) 
  
  #pragma omp target update from(gpu_products.Rmax[0:total_size], \
                                 gpu_products.Fmax[0:total_size]) device(devID)
  
  gputime.memory_transfer.collapse_times += (MPI_Wtime() - tmp);

  #pragma omp parallel for
  for (unsigned int index=0 ; index<total_size ; index++)
    {
      products[index].Rmax = gpu_products.Rmax[index];
      products[index].Fmax = gpu_products.Fmax[index];
    }
  
  #endif // end GPU_OMP_FULL
  
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

  /* PMT measures */
  PMT_CPU_STOP("collapse_time_CPU", ThisTask);
  PMT_GPU_STOP("collapse_time_GPU", devID, ThisTask);
  
  return 0;

} // compute_collapse_times

/* ------------- Tabulated collapse time calculation -----------------------------*/
#if defined(TABULATED_CT)

int initialize_collapse_times(int ismooth, int onlycompute){	

	/* Variable declarations */
	double l1, l2, l3, del, ampl, x, y, interval, ref_interval, deltaf;

	/* Variables used for loop indices, counters, and flags */
	int	id, ix, iy, i, j, dummy, fail; 

	/* File name */
	char fname[LBLENGTH];

	/* The bin size for delta varies like (delta - CT_DELTA0)^CT_EXPO, 
	so it is smallest around CT_DELTA0 if CT_EXPO>1. The bin size 
	is never smaller than CT_MIN_INTERVAL. */

	/*-------------------------- FIRST SMOOTHING RADIUS CASE --------------------------*/ 
	
	if (!ismooth){
		
		/* Dynamically allocation of memory for the delta_vector */
		delta_vector = (double*)malloc(CT_NBINS_D * sizeof(double));

    #if defined(GPU_OMP_FULL)
    #pragma omp target enter data map(alloc: delta_vector[0:CT_NBINS_D]) device(devID)
    #endif

		/* Sets the sampling in density delta */
		/* CT_EXPO value determines the sampling method*/

		if (CT_EXPO == 1)
		  {
		    /* In this case the spacing is even */
		    /* The interval between the elements is 2 times the range of the delta vector divided by the number of bins in the delta vector */
		    interval = 2. * CT_RANGE_D / (double)CT_NBINS_D;
		    for (int id = 0; id < CT_NBINS_D; id++)
		      {
			delta_vector[id] = id * interval - CT_RANGE_D;
		      }
		  } 
		else 
		  {
		    /* In this case the sampling is finer around CT_DELTA0 */		
		    deltaf = pow(CT_SQUEEZE / CT_EXPO, 1. / (CT_EXPO - 1.));
		    if (CT_EXPO == 2)
		      {
			ref_interval = ((log((CT_RANGE_D - CT_DELTA0) / deltaf) + log((CT_RANGE_D + CT_DELTA0) / deltaf))/ CT_EXPO + 
					2. * deltaf / CT_SQUEEZE ) / (CT_NBINS_D - 2.0);
		      } 
		    else 
		      {
			ref_interval = ((pow(CT_RANGE_D - CT_DELTA0, 2. - CT_EXPO) + pow(CT_RANGE_D + CT_DELTA0, 2. - CT_EXPO) 
					 - 2. *pow(deltaf, 2. - CT_EXPO) ) / CT_EXPO / (2. - CT_EXPO) + 2. *deltaf / CT_SQUEEZE ) / (CT_NBINS_D - 2.0);
		      }
	 	
		    del = -CT_RANGE_D;
		    int id = 0;

		    do 
		      {
			/* Assign current value to delta_vector */
			delta_vector[id] = del;

			/* Adjust the interval based on the difference from CT_DELTA0 */
			interval = CT_EXPO * ref_interval * pow(fabs(del - CT_DELTA0), CT_EXPO -1.0);
				
			/* Ensure interval is within the desired range */
			interval = (interval / ref_interval < CT_SQUEEZE ? ref_interval * CT_SQUEEZE : interval);
				
			/* Increment the value */
			del += interval;
			id++;
		      }  
		    while (id<CT_NBINS_D);
		  }
		/* If the current task is the first task */
		if (!ThisTask)
		  {
		    printf("[%s] Grid for interpolating collapse times: CT_NBINS_D=%d, CT_NBINS_XY=%d\n",fdate(),CT_NBINS_D,CT_NBINS_XY);
		  }

		/* Computations are divided among tasks */
		Ncomputations = CT_NBINS_D * CT_NBINS_XY * CT_NBINS_XY;
		bin_x         = CT_RANGE_X /(double)(CT_NBINS_XY);

    #if defined(GPU_OMP_FULL)
    /* Allocate memory for cubic spline array on the host */
    CT_Spline = (CubicSpline **)malloc(CT_NBINS_XY * CT_NBINS_XY * sizeof(CubicSpline *));

    /* Initialize each spline and allocate memory for its internal arrays on the host */
    for (int i = 0; i < CT_NBINS_XY; ++i) {
        for (int j = 0; j < CT_NBINS_XY; ++j) {
            int index = i * CT_NBINS_XY + j; // Flattened index calculation
            CT_Spline[index] = custom_cubic_spline_alloc(CT_NBINS_D); // Allocate
        }
    }

    /* Allocate memory for cubic spline array on the GPU */
    #pragma omp target enter data map(alloc: CT_Spline[0:CT_NBINS_XY * CT_NBINS_XY]) device(devID)

    /* Allocate CT spline and struct members on the GPU */
    for (int i = 0; i < CT_NBINS_XY; ++i) 
    {
      for (int j = 0; j < CT_NBINS_XY; ++j) 
      {
        int index = i * CT_NBINS_XY + j;

        /* Allocate struct members on device */
        #pragma omp target enter data map(alloc: CT_Spline[index][0:1],		              \
					         CT_Spline[index]->d2y_data[0:CT_NBINS_D],    \
                                                 CT_Spline[index]->coeff_a[0:CT_NBINS_D - 1], \
                                                 CT_Spline[index]->coeff_b[0:CT_NBINS_D - 1], \
                                                 CT_Spline[index]->coeff_c[0:CT_NBINS_D - 1], \
                                                 CT_Spline[index]->coeff_d[0:CT_NBINS_D - 1]) device(devID)
        }
    }

  #endif // GPU_OMP_FULL

    /* Allocate memory for CT_table array */
    CT_table = (double*)calloc(Ncomputations, sizeof(double));

    #if defined(GPU_OMP_FULL)
    #pragma omp target enter data map(alloc: CT_table[0:Ncomputations]) device(devID)
    #endif

    #ifdef DEBUG
    #ifdef PRINTJUNK
		/* Create a debug file for writing debug output */
		char filename[50];
		sprintf(filename, "debug.task%d.txt", ThisTask);
		JUNK = fopen(filename, "w");
    #endif
    #endif

	} // fine (!ismooth)

	/*--------------------------  NOT FIRST SMOOTHING RADIUS CASE --------------------------*/
	
	/* Check if collapse times should be read from a file and not in the "onlycompute" mode */
	if (strcmp(params.CTtableFile,"none") && !onlycompute)
	  {
	    /* In the case where ismooth is false (i.e., the current smoothing iteration is not being performed), Task 0 reads the file */
	    if (!ismooth)
	      {
		if (!ThisTask)
		  {
		    /* Performs consistency checks using check_CTtable_header */
		    CTtableFilePointer = fopen(params.CTtableFile , "r");
		    fail               = check_CTtable_header(CTtableFilePointer);
		  }

		/* Broadcast the fail status to all tasks to ensure consistency across the parallel execution */
		MPI_Bcast(&fail, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		if (fail)
		  return 1;
	      }
	    if (!ThisTask)
	      {
		/* Read collapse times from the file */
		fread(&dummy, sizeof(int), 1, CTtableFilePointer);
		fread(CT_table, sizeof(double), Ncomputations, CTtableFilePointer);
	      }
			
	    /* Broadcast CT_table array to all tasks */
	    MPI_Bcast(CT_table, Ncomputations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	    /* If the current task is task 0 and ismooth is true for the last smoothing iteration */
	    if (!ThisTask && ismooth == Smoothing.Nsmooth - 1)
	      fclose(CTtableFilePointer);
	  } 
	else 
	  {
	    /* The collapse time table is computed for each smoothing radius */
	    ampl = sqrt(Smoothing.Variance[ismooth]);

	    /* Create a map of computations for all tasks by dividing the total number of computations (Ncomputations) among the tasks */
	    int *counts = (int*)malloc(NTasks*sizeof(int));
	    int *displs = (int*)malloc(NTasks*sizeof(int));
    
    #if defined(GPU_OMP_FULL)
    #pragma omp target enter data map(alloc: counts[0:NTasks], \
                                             displs[0:NTasks]) device(devID)
    #endif

	    int bunch     = Ncomputations / NTasks;
	    int remainder = Ncomputations % NTasks;

	    for (int i = 0, ix = 0; i < NTasks; i++){
	      counts[i] = bunch + (i < remainder);
	      displs[i] = bunch*i + ((remainder > i) ? i : remainder );
		}

    #if defined(GPU_OMP_FULL)
    
    /* Update variables on the GPU */
    #pragma omp target update to(delta_vector[0:CT_NBINS_D], displs[0:NTasks], counts[0:NTasks]) device(devID)
    #endif

    /* Computations for this task */
    /* Offload the loop to the GPU */
    #pragma omp target teams distribute parallel for private(ix, iy, x, y, l1, l2, l3) device(devID)
	    for (int i = displs[ThisTask]; i < displs[ThisTask] + counts[ThisTask]; ++i)
	      {
		int id            = i % CT_NBINS_D;
		ix                = (i / CT_NBINS_D) % CT_NBINS_XY;
		iy                = i / CT_NBINS_D / CT_NBINS_XY;
		double bin_x_test = CT_RANGE_X /(double)(CT_NBINS_XY);
		x                 = ix * bin_x_test;
		y                 = iy * bin_x_test;
		l1                = (delta_vector[id] + 2.*x +    y) / 3.0*ampl;
		l2                = (delta_vector[id] -    x +    y) / 3.0*ampl;
		l3                = (delta_vector[id] -    x - 2.*y) / 3.0*ampl;
		CT_table[i]       = ell_gpu(ismooth,l1,l2,l3);
	      } // end omp kernel

           #if defined(GPU_OMP_FULL)
             #pragma omp target update from(CT_table[0:Ncomputations]) device(devID)
	   #endif
	   
	    /* Synchronize computed CT_table across all tasks */
	    MPI_Allgatherv(MPI_IN_PLACE, counts[ThisTask], MPI_DOUBLE, CT_table, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

           #if defined (GPU_OMP_FULL)
              #pragma omp target update to(CT_table[0:Ncomputations]) device(devID)
	   #endif 
	    
	    /* Free memory */
	    free(displs);
	    free(counts);

            /* Clean up GPU data mapping */
            #if defined(GPU_OMP_FULL)
                #pragma omp target exit data map(delete: displs[0:NTasks], counts[0:NTasks]) device(devID)
            #endif 
	    
	    /* Write the table on a file */			
	    if (!ThisTask)
	      {
		if (onlycompute)
		  {
		    strcpy(fname, params.CTtableFile);
		  } 
		else 
		  {
		    sprintf(fname,"pinocchio.%s.CTtable.out",params.RunFlag);
		  }
		if (!ismooth)
		  {
		    /* Create a new file and write header if ismooth is false */
		    CTtableFilePointer = fopen(fname,"w");
		    write_CTtable_header(CTtableFilePointer);
		  } 
		else 
		  {
		    /* Append to the existing file if ismooth is true */
		    CTtableFilePointer = fopen(fname,"a");
		  }
      
#ifdef ASCII_TABLE
      /* Write table contents in ASCII format */
			fprintf(CTtableFilePointer,"################\n# smoothing %2d #\n################\n",ismooth);

			for (i = 0; i < Ncomputations; i++){
				id = i % CT_NBINS_D;
				ix = (i / CT_NBINS_D) % CT_NBINS_XY;
				iy = i / CT_NBINS_D / CT_NBINS_XY;
				x  = ix * bin_x;
				y  = iy * bin_x;
				l1 = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
				l2 = (delta_vector[id] -    x +    y)/3.0*ampl;
				l3 = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

				fprintf(CTtableFilePointer, " %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20g\n",
					    ismooth, id, ix, iy, delta_vector[id], x, y, l1, l2, l3, CT_table[i]);
			}
#else
			/* Write table contents in binary format */
			fwrite(&ismooth, sizeof(int), 1, CTtableFilePointer);
			fwrite(CT_table, sizeof(double), Ncomputations, CTtableFilePointer);
#endif
			fclose(CTtableFilePointer);
		}

	}

#if defined(GPU_OMP_FULL)

	// Fill each spline  and offload their internal values on the device
  for (int i = 0; i < CT_NBINS_XY; ++i) 
  {
    for (int j = 0; j < CT_NBINS_XY; ++j) 
    {
      int index = i * CT_NBINS_XY + j; // Flattened index calculation
      // Initialize the cubic spline with the allocated structure
      custom_cubic_spline_init(CT_Spline[index], delta_vector, &CT_table[i * CT_NBINS_D + j * CT_NBINS_D * CT_NBINS_XY], CT_NBINS_D);
    }
  } 
#endif // GPU_OMP_FULL

  return 0;
}


/* ------------------------------------- Reset the collapse times ------------------------------------- */

int reset_collapse_times(int ismooth){
	
#ifdef DEBUG
	if (counter)
	{
#ifdef BILINEAR_SPLINE
		flag = 2;  // Use bilinear spline interpolation
#endif
		double aa      = ave_diff / (double)counter;
		double std_dev = sqrt(var_diff / (double)counter) - aa * aa;

		/* Print debug information */
		printf("%d   %2d  %10d   %3d %3d %5.2f %5.2f   %20g  %20g DEBUG\n", flag, ismooth, counter, CT_NBINS_D, CT_NBINS_XY, CT_EXPO, CT_SQUEEZE, aa, std_dev);
	}

	/* Reset variables for the next iteration */
	ave_diff = 0.0;
	var_diff = 0.0;
	counter  = 0;
#endif //DEBUG
	return 0;
}

/* ------------------------------------- Comparison routine ------------------------------------- */

int compare_search(const void *A, const void *B){
	if (*(double*)A >= *(double*)B && *(double*)A < *((double*)B+1))
		return 0;
	if (*(double*)A < *(double*)B )
		return -1;
	return 1;
}

/* ------------------------------------- Interpolation of collapse times ------------------------------------- */

inline double interpolate_collapse_time(const int ismooth, const double l1, const double l2, const double l3){
 
	/* Variable declarations and assignment */
	double ampl   = sqrt(Smoothing.Variance[ismooth]);
	double d      = (l1 + l2 + l3) / ampl;
	double x      = (l1 - l2) / ampl;
	double y      = (l2 - l3) / ampl;
	double my_bin = CT_RANGE_X /(double)(CT_NBINS_XY);
	int ix        = (int)(x / my_bin);
	int iy        = (int)(y / my_bin);

	ix = (ix >= CT_NBINS_XY - 1) ? CT_NBINS_XY - 2 : (ix < 0) ? 0 : ix;
	iy = (iy >= CT_NBINS_XY - 1) ? CT_NBINS_XY - 2 : (iy < 0) ? 0 : iy;


#ifdef BILINEAR_SPLINE
    
	/* Calculate the offsets */
	double dx = x / my_bin - ix;
	double dy = y / my_bin - iy;

#if defined(GPU_OMP_FULL)

  // GPU Interpolation
  double val00 = custom_cubic_spline_eval(CT_Spline[ix * CT_NBINS_XY + iy], d, CT_NBINS_D);
  double val10 = custom_cubic_spline_eval(CT_Spline[(ix + 1) * CT_NBINS_XY + iy], d, CT_NBINS_D);
  double val01 = custom_cubic_spline_eval(CT_Spline[ix * CT_NBINS_XY + (iy + 1)], d, CT_NBINS_D);
  double val11 = custom_cubic_spline_eval(CT_Spline[(ix + 1) * CT_NBINS_XY + (iy + 1)], d, CT_NBINS_D);
  
  return  (1.0 - dx) * (1.0 - dy) * val00 +
          dx * (1.0 - dy)         * val10 +
          (1.0 - dx) * dy         * val01 +
          dx * dy                 * val11;

#endif
#endif

}

/*-------------------------- Consistency check of the header of CT Table file with the current run --------------------------------------------------------*/

int check_CTtable_header(){

	/* Checks that the information contained in the header is consistent with the run 
	NB: this is done by Task 0 */

	int fail = 0;
    int    dummy;
    double fdummy;

	/* Check the CT table type based on compilation flags */
    fread(&dummy, sizeof(int), 1, CTtableFilePointer);
	
#ifdef ELL_CLASSIC 
	if (dummy != 1)
  { /* classic ellipsoidal collapse, tabulated */	
    printf("ERROR: CT table not constructed for ELL_CLASSIC, %d\n", dummy);
    fail = 1;
	}
#endif
#ifdef  ELL_SNG
#ifndef MOD_GRAV_FR
	if (dummy != 3)
  { /* standard gravity, numerical ellipsoidal collapse */
		printf("ERROR: CT table not constructed for ELL_SNG and standard gravity, %d\n", dummy);
    fail = 1;
	}
#else
	if (dummy != 4){ 
    /* f(R) gravity, numerical ellipsoidal collapse */
		printf("ERROR: CT table not constructed for ELL_SNG and MOD_GRAV_FR, %d\n", dummy);
    fail = 1;	
	}
#endif
#endif

	/* Check the CT table parameters */
	fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
	if (fabs(fdummy - params.Omega0) > 1.e-10)
  {
		printf("ERROR: CT table constructed for the wrong Omega0, %f in place of %f\n", fdummy, params.Omega0);
    fail = 1;
	}
  fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
  if (fabs(fdummy - params.OmegaLambda) > 1.e-10)
  {
    printf("ERROR: CT table constructed for the wrong OmegaLambda, %f in place of %f\n", fdummy, params.OmegaLambda);
    fail = 1;
  }
	fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
  if (fabs(fdummy - params.Hubble100) > 1.e-10)
  {
    printf("ERROR: CT table constructed for the wrong Hubble100, %f in place of %f\n", fdummy, params.Hubble100);
    fail = 1;
  }

	/* Check the CT table size and sampling */
  fread(&dummy, sizeof(int), 1, CTtableFilePointer);
  if (dummy != Ncomputations)
  {
    printf("ERROR: CT table has the wrong size, %d in place of %d\n", dummy, Ncomputations);
    fail = 1;
  }
  fread(&dummy, sizeof(int), 1, CTtableFilePointer);
  if (dummy != CT_NBINS_D)
  {
    printf("ERROR: CT table has the wrong density sampling, %d in place of %d\n", dummy, CT_NBINS_D);
    fail = 1;
  }
  fread(&dummy, sizeof(int), 1, CTtableFilePointer);
  if (dummy != CT_NBINS_XY)
  {
    printf("ERROR: CT table has the wrong x and y sampling, %d in place of %d\n", dummy, CT_NBINS_XY);
    fail = 1;
  }
    
  return fail;
}

/*-------------------------- Write the infos into the header of CT_Table_file --------------------------------------------------------*/

void write_CTtable_header(){

  int dummy;

#ifdef ELL_CLASSIC 
	/* Classic ellipsoidal collapse, tabulated */
  dummy = 1;   
  fwrite(&dummy, sizeof(int), 1, CTtableFilePointer);
#endif
#ifdef ELL_SNG
#ifndef MOD_GRAV_FR
	/* Standard gravity, numerical ellipsoidal collapse */
  dummy = 3; 
  fwrite(&dummy, sizeof(int), 1, CTtableFilePointer);
#else
	/* f(R) gravity, numerical ellipsoidal collapse */
  dummy = 4;   
  fwrite(&dummy, sizeof(int), 1, CTtableFilePointer);
#endif
#endif

	/* Write Omega0 */
	fwrite(&params.Omega0, sizeof(double), 1, CTtableFilePointer);  

	/* Write OmegaLambda */
  fwrite(&params.OmegaLambda, sizeof(double), 1, CTtableFilePointer); 

	/* Write Hubble100 */
  fwrite(&params.Hubble100, sizeof(double), 1, CTtableFilePointer);  

	/* Write Ncomputations */
  fwrite(&Ncomputations, sizeof(int), 1, CTtableFilePointer); 

	/* Write CT_NBINS_D */
  dummy = CT_NBINS_D;
  fwrite(&dummy, sizeof(int), 1, CTtableFilePointer);    

	/* Write CT_NBINS_XY */
  dummy = CT_NBINS_XY;
  fwrite(&dummy, sizeof(int), 1, CTtableFilePointer);         
}
#endif // TABULATED_CT
#endif // GPU_OMP || FULL_GPU_OMP
