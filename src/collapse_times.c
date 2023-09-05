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


#include "pinocchio.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <immintrin.h>

/*------------------------------------------------------- Macros declaration --------------------------------------------------------*/


#define SMALL ((double)1.e-20)


/*------------------------------------------------------- Functions definition -------------------------------------------------------*/


/* Ellipsoidal solution at 3rd perturbative order in two different ways */

__attribute__((always_inline)) double ell_classic               ( int, double, double, double);  
__attribute__((always_inline)) double ell_sng                   ( int, double, double, double);  

/* Collapase time calculation */

__attribute__((always_inline)) double ell                       ( int, double, double, double); 

/* Calculation of inverse collpase time */

__attribute__((always_inline)) double inverse_collapse_time     ( int, double *, double *, double *, double *, int *);

/* Interpolation of collpase time from a table containing collpase times */

#ifdef TABULATED_CT

__attribute__((always_inline)) double interpolate_collapse_time (int, double, double, double);

#endif 

/* Order inverse collpase time */

__attribute__((always_inline)) void ord                         (double *,double *,double *);


/*------------------------------------------------------- Global Variables declaration ----------------------------------------------------*/


/* Declaring the replicated variables that will be use by each single thread */
/* NOTE: Threadprivate directive specifies that variables are replicated, with each thread having its own copy. It's a declarative directive */

#if defined( _OPENMP )

/* Cpu_time for elliptical collapse */

double cputime_ell;
#pragma omp threadprivate(cputime_ell) // threadprivate directive specifies that variables are replicated, with each thread having its own copy


/* It will be use as a fail "flag" indicating that for some reason the calculation of collapse time failed */

int fails;
#pragma omp threadprivate(fails)

/* If there's no _OPENMP macro then define cputime_ell */

#else
#define cputime_ell            cputime.ell
#endif

/* Declaring smoothing radius */

int ismooth;


/*------------------------------------------------------- Functions implementation --------------------------------------------------------*/

/* Classical ellipsoidal collapse solution */

inline double ell_classic(int ismooth, double l1, double l2, double l3) {

    /*
    This routine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
    */

    /* Local variables declaration */
    double ell;
    double del = l1 + l2 + l3;   
    double det = l1 * l2 * l3;    

    /* Vanishing lambda1 eigenvalue case */
    if (fabs(l1) < SMALL) {
        ell = -0.1;
    }
    /* Not vanishing lambda1 eigenvalue case */
    else {
        double den = det / 126. + 5. * l1 * del * (del - l1) / 84.;
        /* Check 1st perturbative order conditions */
        if (fabs(den) < SMALL) {
            if (fabs(del - l1) < SMALL) {
                if (l1 > 0.0) {
                    /* Zel'dovich approximation */
                    ell = 1. / l1;
                } else {
                    ell = -.1;
                }
            } else {
                /* Check 2nd perturbative order conditions */
                double dis = 7. * l1 * (l1 + 6. * del);
                if (dis < 0.0) {
                    ell = -.1;
                } else {
                    /* 2nd order solution */
                    ell = (7. * l1 - sqrt(dis)) / (3. * l1 * (l1 - del));
                    if (ell < 0.) {
                        ell = -.1;
                    }
                }
            }
        } else {
            /* 3rd order perturbative solution. For more details about the equations implemented, see Monaco 1996a */
            
            /* Intermediate values */
            double rden = 1.0 / den;
            double a1 = 3. * l1 * (del - l1) / 14. * rden;
            double a1_2 = a1 * a1;
            double a2 = l1 * rden;
            double a3 = -1.0 * rden;

            /* The collapse time b_c will be a combination of R, Q, and D == den */
            double q = (a1_2 - 3. * a2) / 9.;
            double r = (2. * a1_2 * a1 - 9. * a1 * a2 + 27. * a3) / 54.;
            double r_2_q_3 = r * r - q * q * q;

            /* Check 3rd perturbative order conditions */

            /* ---------------- Case 1 --------------- */
            /* If R^2 - Q^2 > 0, which is valid for spherical and quasi-spherical perturbations */
            
            /* 3rd order solution */
            if (r_2_q_3 > 0) {
                double fabs_r = fabs(r);
                double sq = pow(sqrt(r_2_q_3) + fabs_r, 0.333333333333333);
                ell = -fabs_r / r * (sq + q / sq) - a1 / 3.;
                if (ell < 0.) {
                    ell = -.1;
                }
            }

            /* ---------------- Case 2 --------------- */
            /* The solution has to be chosen as the smallest non-negative one between s1, s2, and s3 */
            
            /* 3rd order solution */
            else {
                double sq = 2 * sqrt(q);
                double inv_3 = 1.0 / 3;
                double t = acos(2 * r / q / sq);
                double s1 = -sq * cos(t * inv_3) - a1 * inv_3;
                double s2 = -sq * cos((t + 2. * PI) * inv_3) - a1 * inv_3;
                double s3 = -sq * cos((t + 4. * PI) * inv_3) - a1 * inv_3;
                if (s1 < 0.) {
                    s1 = 1.e10;
                }
                if (s2 < 0.) {
                    s2 = 1.e10;
                }
                if (s3 < 0.) {
                    s3 = 1.e10;
                }
                ell = (s1 < s2 ? s1 : s2);
                ell = (s3 < ell ? s3 : ell);
                if (ell == 1.e10) {
                    ell = -.1;
                }
            }
        }
    }
    
    if (del > 0. && ell > 0.) {
        double inv_del = 1.0 / del;
        ell += -.364 * inv_del * exp(-6.5 * (l1 - l2) * inv_del - 2.8 * (l2 - l3) * inv_del);
    }

    return ell;
}

/* Ellipsoidal collapse following Nadkarni-Ghosh & Singhal (2016) */

/* We follow the dynamics of triaxial collapse in terms of eigenvalues of the deformation tensor (lambda_a), the velocity derivative tensor (lambda_v) and the gravity Hessian 
(lambda_d). The idea is: starting from BM96, where the dynamic is characterized by the evolution of the three principal axes (the physical coordinates is r=a_i(t)q(i = 1,2,3)). 
In this case the evolution of the ellipse is completely determined once six parameters are known: the three axes lengths and their velocities at some initial epoch a_init. 
An alternate description of the ellipse can be given by a set of nine dimensionless parameters lambda_a, lambda_v, lambda_d. In this description when an axis is collpasing 
lamda_a ----> 1, whereas for an expanding axes lambda_a -----> -inf.
delta = lambda_d_1 + lambda_d_2 + lambda_d_3 
The nine (dimensionless) eigenvalues completely characterize the density, velocity and shape perturbations in this model of ellipsoidal collapse. 
For the net system for the nine eigenvalues see Nadkarni-Ghosh & Singhal (2016) pag. 5 */


/* We have then to solve a system of differential equations          */
/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */

int sng_system(double t, const double y[], double f[], void *sng_par) {

    /* Local variables declaration */
	int i,j;
	double sum;

	/* Needed cosmological parameters calculation at the given z */
	double omegam = OmegaMatter(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */
	double omegal = OmegaLambda(1. / t-1.);  /* Cosmological mass density parameter as a function of redshift DIMENSIONLESS. From cosmo.c */

	/* y array contains the eigenvalues lambda_a, lambda_v and lambda_d */
	/* In this case y[6] = lambda_d_1, y[7] = lambda_d_2, y[8] = lambda_d_3 */
	double delta = y[6] + y[7] + y[8];

    /* Calculating r.h.s sums of lambda_d_i evolution */
	/* Loop from 0 -----> 3 beacause lambda_i has 3 component (it's the same for lambda_a and lambda_v) */
	
	/* ---------------- Sum  i != j --------------- */
	 for (i = 0; i < 3; i++) {
        sum = 0.;
        for (j = 0; j < 3; j++) {
            if (i == j || y[i] == y[j]) {
                continue;
            } else {
                sum += (y[j + 6] - y[i + 6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i + 3]) -
                       (1. - y[j]) * (1. - y[j]) * (1. + y[j + 3])) /
                       ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
            }
        }

		/* ---------------- Gathering r.h.s of lambda_a_i equation --------------- */
        f[i] = (y[i + 3] * (y[i] - 1.0)) / t;

		  /* ---------------- Gathering r.h.s of lambda_v_i equation --------------- */
        /* NOTE: here is the part where models of modified gravity can be introduced */
        f[i + 3] = (0.5 * (y[i + 3] * (omegam - 2.0 * omegal - 2.0) 
	
#ifdef MOD_GRAV_FR
                       - 3.0 * omegam * y[i + 6] * (1. + ForceModification(*(double *)sng_par, t, delta))
#else
                       - 3.0 * omegam * y[i + 6]
#endif
                       - 2.0 * y[i + 3] * y[i + 3])) / t;

		/* ---------------- Gathering r.h.s of lambda_d_i equation --------------- */
        f[i + 6] = ((5. / 6. + y[i + 6]) *
                    ((3. + y[3] + y[4] + y[5]) - (1. + delta) / (2.5 + delta) * (y[3] + y[4] + y[5])) -
                    (2.5 + delta) * (1. + y[i + 3]) + sum) / t;
    }

    return GSL_SUCCESS;
}


/* Solving the ODE system for the ellipsoidal collapse following SNG */

inline double  ell_sng(int ismooth, double l1, double l2, double l3) {

	/* Needed variables for the integrations step */
	double ode_param, hh = 1.e-6;  // hh = initial step-size, ode_param = arbitrary parameters of the system == Smoothing_radius in our case

    /* Setting integration time */
	double amin = 1.e-5, amax = 5.0;
	double mya = amin;
    
	/* Step type: Explicit embedded Runge-Kutta-Fehlberg (4, 5) method */
	const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;

	/* Newly allocated instance of a stepping function of type T for a system of dim = 9 */
	gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9); 

	/* The control function examines the proposed change to the solution produced by a stepping function and attempts to determine the optimal step-size for a user-specified level of error.
	The standard control object is a four parameter heuristic based on absolute and relative errors eps_abs and eps_rel,and scaling factors a_y and a_dydt 
	for the system state y(t) and derivatives y'(t) respectively. */
	gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);

	/* The evolution function combines the results of a stepping function and control function to reliably advance the solution forward one step using an acceptable step-size.
	This function returns a pointer to a newly allocated instance of an evolution function for a system of dim = 9 */
	gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);

	/* A system of equations is defined using the gsl_odeiv2_system datatype 
	This data type defines a general ODE system with arbitrary parameters. In this case we need the r.h.s of the system that will be solved: sng_system
	The vector of derivatives elements: jac
	The dimension of the system of equations: 9 
	A pointer to the arbitrary parameters of the system: (void*)&ode_param */
	gsl_odeiv2_system   ode_sys   = {sng_system, jac, 9, (void*)&ode_param};

/* GrowingMode is linear growing mode, interpolated on the grid. See cosmo.c for details (double GrowingMode(double z, double k))*/
#ifdef SCALE_DEPENDENT

	double D_in = GrowingMode(1./amin-1.,params.k_for_GM/Smoothing.Radius[ismooth]*params.InterPartDist); 

#else

	double D_in =  GrowingMode(1./amin-1.,1./Smoothing.Radius[ismooth]);

#endif

	double y[9] = {l1*D_in, l2*D_in, l3*D_in,
		          l1*D_in/(l1*D_in - 1.), l2*D_in/(l2*D_in - 1.), l3*D_in/(l3*D_in - 1.),
		          l1*D_in, l2*D_in, l3*D_in};
    
/* Assignment of ode_param in two different cases */
#ifdef MOD_GRAV_FR

	if (ismooth<Smoothing.Nsmooth-1)
		{
		ode_param = Smoothing.Radius[ismooth];
		}
	else
		{
		ode_param = Smoothing.Radius[ismooth-1];
		}
#endif

    /*---------------- Integration step --------------- */
	double olda = mya; 
	double oldlam = y[0]; 

	while (mya < amax) {
        /* gsl_odeiv2_evolve_apply advances the system from time t and position y using the stepping function selected before i.e. rkf45. 
		The new time and position are stored in t and y on output*/
        int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &mya, amax, &hh, y);

        /* Correct integration check */
        if (status != GSL_SUCCESS) {
            printf("ERROR on task %d: integration of cosmological quantities failed\n", ThisTask);
            fflush(stdout);
            return -1;
        }

        /* ---------------- Update ellipsoid axis --------------- */
        if (y[0] >= 0.99999) {

            return olda + (1. - oldlam) * (mya - olda) / (y[0] - oldlam);

        }
    }

    /* In this case the ellipsoid does not collapse */
    return 0;
}

/* --------------------------------------------------- Modified gravity model --------------------------------------------------------*/

#ifdef MOD_GRAV_FR

double ForceModification(double size, double a, double delta) {

    double ff        = 4. * params.OmegaLambda / params.Omega0;
    double thickness = FR0 / params.Omega0 / pow(H_over_c * size, 2.0) *
                       pow(a, 7.) * pow((1. + delta), -1. / 3.) *
                       (pow((1.0 + ff) / (1.0 + ff * pow(a, 3.)), 2.0) -
                       pow((1.0 + ff) / (1.0 + delta + ff * pow(a, 3.)), 2.0));

    double F3 = (thickness * (3. + thickness * (-3. + thickness)));
    if (F3 < 0.) {
        F3 = 0.;
    }
    return (F3 < 1. ? F3 / 3. : 1. / 3);
}

#endif


/*---------------------------------------- Calculation of b_c == growing mode at collapse time ---------------------------------------*/

inline double ell(int ismooth, double l1, double l2, double l3) {

#ifdef ELL_CLASSIC

    double bc = ell_classic(ismooth, l1, l2, l3); 

    if (bc > 0.0) {
        return 1. + InverseGrowingMode(bc, ismooth);
    } else {
        return 0.0;
    }
#endif

#ifdef ELL_SNG

    double bc = ell_sng(ismooth, l1, l2, l3);

    if (bc > 0.0) {
        return 1. / bc;
    } else {
        return 0.0;
    }
#endif
}

/* ------------------------------------  Computation of collapse time i.e. F = 1 + z_collapse and variance ------------------------------------ */

inline int compute_collapse_times(int ismooth) {    

	 /*---------------------DEBUG---------------------------------------*/

    // Open a file for writing collapse times

    FILE *collapseTimeFile = fopen("collapse_times.txt", "w");
    if (collapseTimeFile == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

    /*----------------------------------------------------------------*/

    /* Local average and variance declaration */
    double local_average, local_variance;

    /* Initialization */
    local_variance = 0.0;
    local_average = 0.0;

	/* Initialize Fmax, Rmax and velocity if it is the first smoothing 
	The initialization is done in a parallel way 
	Each thread initializes its corresponding array elements
	This can improve the initialization performance when executed on systems with multiple processors or cores */

	if (!ismooth)

		/* Loop on all particles */ 
		#pragma omp parallel for 
		for (int i = 0; i < MyGrids[0].total_local_size; i++){

			/*----------- Common initialization ----------- */		
			products[i].Fmax   = -10.0;
			products[i].Rmax   = -1;

			/*----------- Zel'dovich case ----------- */
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
			
	/* Declaration and initialization of variables when the openmp option is ON	
	Inside the parallel region, each thread initializes its own local variables i.e.
	mylocal_average, mylocal_variance, cputime_invcoll, and cputime_ell 
	These variables are private to each thread and are not shared among threads, ensuring data consistency.*/
     
#if defined( _OPENMP )

	double invcoll_update = 0, ell_update = 0;

	#pragma omp parallel
	{

	/* Thread-specific variables declaration */	
	double  mylocal_average, mylocal_variance;
    double  cputime_invcoll;
		
	/* Support variables */	
	int     fails = 0;
	int     pid = omp_get_thread_num();
		
	/* Initializitation */
	mylocal_average     = 0;
	mylocal_variance    = 0;
	cputime_invcoll     = 0;
	cputime_ell         = 0;

#else

	/* If openmp is OFF 
    This approach is used to avoid conditional compilation directives within the code and maintain consistency in variable
	names across the parallel and serial versions of the code. 
	When OpenMP is disabled, the program uses the same variable names, but they refer to the original variables directly */

#define mylocal_average   local_average
#define mylocal_variance  local_variance
#define cputime_invcoll   cputime.invcoll
#define cputime_ell       cputime.ell
	
#endif

/*----------------- Calculation of variance, avarage and collapse time ----------------*/

	double tmp = MPI_Wtime();

/* When you use a parallel region, OpenMP will automatically wait for all threads to finish before execution */
/* If you do not need synchronization after the loop, you can disable it with nowait */
#ifdef _OPENMP
#pragma omp for nowait
#endif

	/* Loop on all particle */
	for (int index = 0; index < MyGrids[0].total_local_size; index++){

		/* Computation of second derivatives of the potential i.e. the gravity Hessian */
		double diff_ten[6]; 
		for (int i = 0; i < 6; i++){
			diff_ten[i] = second_derivatives[0][i][index];
		}

		/* Computation of the variance of the linear density field */
		double delta      = diff_ten[0] + diff_ten[1] + diff_ten[2];
		mylocal_average  += delta;
		mylocal_variance += delta*delta;

		/* Computation of the collapse time */
		double lambda1,lambda2,lambda3;
		int    fail;
		double Fnew = inverse_collapse_time(ismooth, diff_ten, &lambda1, &lambda2, &lambda3, &fail); 

		/*---------------------------DEBUG----------------------------------*/

		fprintf(collapseTimeFile, "%d %d %f %f %f %g\n", index, ismooth, lambda1, lambda2, lambda3, Fnew);

		/*-----------------------------------------------------------------*/

		/* Fail check during computation of the collapse time */
		if (fail){
			printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
			fflush(stdout);

/* Workflow control */
#if !defined( _OPENMP )		    		    
			return 1;
#else
			fails = 1;
#endif
		}
    
		/* Updating collapse time */
		if (products[index].Fmax < Fnew){
				products[index].Fmax = Fnew;
				products[index].Rmax = ismooth;
		}	      
	}

	/* CPU collapse time */
	cputime_invcoll += MPI_Wtime() - tmp;
	

	/* ------------------- Updates local variance and local average using OpenMP atomic operations  ------------------*/
	/* local_variance and local_average are shared variables that are updated by each thread (mylocal_variance and mylocal_average)*/

#if defined( _OPENMP )

	#pragma omp atomic
	local_variance += mylocal_variance;

	#pragma omp atomic
	local_average += mylocal_average;  

	/* ------------------------------------ Updates total CPU collapse time -----------------------------------*/
	/* cputime_invcoll and cputime_ell are local variables that are updated by each thread and then accumulated using atomic operations */

	#pragma omp atomic
	invcoll_update += cputime_invcoll;

	#pragma omp atomic	
	ell_update     += cputime_ell;

#endif	

#if defined( _OPENMP ) 
	}
#endif

	int all_fails = 0; 

	/* -------------------------- Updates cpu collapse time for a single threads -----------------------------*/

#if defined (_OPENMP)

	#pragma omp atomic
	cputime.invcoll += invcoll_update/internal.nthreads_omp; 

	#pragma omp atomic	
	cputime.ell     += ell_update/internal.nthreads_omp;

	/* fails is a local variable that indicates the number of failures in the computation. 
	Each thread updates fails, and then the values are accumulated using atomic operations */

	#pragma omp atomic
	all_fails       += fails;

#endif
    
	/* Fail check during computation of the inverse collapse time */
	/* If there were failures, an error message is printed and the function returns 1 */
	if (all_fails){
		printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
		fflush(stdout);
		return 1;
	}

	/* Calculating obtained variance and avarage */		
	double global_variance = 0.0;
	double global_average  = 0.0;
    
	/* Calculates the global variance and average by reducing the values of local_variance and local_average across all MPI tasks */
	MPI_Reduce(&local_variance, &global_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&local_average , &global_average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

	/* If the current task is the root task (task 0), it divides the global variance and average by the total number of data points (MyGrids[0].Ntotal) 
	to obtain the average values per data point */
	if (!ThisTask){
		global_variance /= (double)MyGrids[0].Ntotal;
		global_average  /= (double)MyGrids[0].Ntotal;
	}
	
	/* Broadcasts the value of global_variance from the root task to all other tasks using MPI_Bcast */
	MPI_Bcast(&global_variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	/* Stores the true variance*/
	Smoothing.TrueVariance[ismooth] = global_variance;

	return 0;
}


/* ------------------------------------  Diagonalization of the potential Hessian ------------------------------------ */


inline double inverse_collapse_time(int ismooth, double * restrict deformation_tensor, double * restrict x1, double * restrict x2, double * restrict x3, int * restrict fail) {
    
	/* Local variables declaration */
	/* mu1, mu2 and mu3 are the principal invariants of the 3x3 tensor of second derivatives */
	double mu1, mu2, mu3;
	double dtensor[6] = { deformation_tensor[0], deformation_tensor[1] ,
					 	  deformation_tensor[2], deformation_tensor[3] ,
					 	  deformation_tensor[4], deformation_tensor[5] };
	*fail = 0;

	/*  Performs diagonalization of the tensor by calculating the values of mu1, mu2, and mu3 */
	mu1 = dtensor[0] + dtensor[1] + dtensor[2];
	double mu1_2 = mu1 * mu1;
	mu2 = 0.5 * mu1_2;

    /* By define add[3] inside and outside the scope, the code ensures that any previous variable with the same name add is not affected
	and a new variable add is created within the block. This is useful when the same variable name needs to be used again in the code without
	conflicting with any previous variable of the same name */
	{
		double add[3];
		add[0] = dtensor[0]*dtensor[0];
		add[1] = dtensor[1]*dtensor[1];
		add[2] = dtensor[2]*dtensor[2];
		mu2 -= 0.5 * (add[0] + add[1] + add[2]);
	}

	/* Redefinition outside the block to calculate the squares of dtensor[3], dtensor[4], and dtensor[5].
	 This allows separate storage for these intermediate calculations, preventing any unintended interference or confusion between the two sets of calculations.*/
	double add[3];
	add[0] = dtensor[3]*dtensor[3];
	add[1] = dtensor[4]*dtensor[4];
	add[2] = dtensor[5]*dtensor[5];
	mu2 -= add[0] + add[1] + add[2];

    /* mu3 calculation */
	mu3 =  dtensor[0]*dtensor[1]*dtensor[2] +
		2.*dtensor[3]*dtensor[4]*dtensor[5] -
		dtensor[0]*add[2] -
		dtensor[1]*add[1] -
		dtensor[2]*add[0];

	/* Check if the tensor is already diagonal */
	double q;
	q = ( mu1_2 - 3.0*mu2 ) /9.0;

	if (q == 0.){
		/* In this case the tensor is already diagonal */
		*x1 = dtensor[0];
		*x2 = dtensor[1];
		*x3 = dtensor[2];
	} else {  
		/* The solution has to be chosen as the smallest non-negative one between x1, x2 and x3 */
		double r = -(2.*mu1_2*mu1 - 9.0*mu1*mu2 + 27.0*mu3)/54.;
        
		/* Fail check */
		if (q*q*q < r*r || q<0.0){
			return -10.0;
		}

        /* Calculating x1, x2, x3 solution in the same way as it done in ell_classic */
		double sq = 2 * sqrt(q);
		double t = acos(2*r/q/sq);
		double inv_3 = 1.0/3.0;
		*x1 = -sq*cos(t * inv_3) + mu1 * inv_3;
		*x2 = -sq*cos((t + 2.*PI) * inv_3 )+ mu1 * inv_3;
		*x3 = -sq*cos((t + 4.*PI) * inv_3 )+ mu1 * inv_3;
		
	}

	/* Ordering and inverse collapse time */
	ord(x1, x2, x3);
	
	/* -------------------------------- Tabulated collapse time case -------------------------------- */

#ifdef TABULATED_CT
	double t = interpolate_collapse_time(ismooth,*x1,*x2,*x3);

#ifdef DEBUG
	double t2 = ell(ismooth,*x1,*x2,*x3);
	if (t2 > 0.9){
		ave_diff += t - t2;
		var_diff += pow(t - t2, 2.0);
		counter += 1;
#ifdef PRINTJUNK
		double ampl = sqrt(Smoothing.Variance[ismooth]);
		fprintf(JUNK,"%d  %f %f %f   %f %f %f   %20g %20g %20g\n",ismooth, (*x1 + *x2 + *x3)/ampl, (*x1 - *x2)/ampl, (*x2 - *x3)/ampl, *x1, *x2, *x3, t, t2, t - t2);
#endif
	}
#endif

#else 
	/* Final computation of collpse time*/
	double t = ell(ismooth,*x1,*x2,*x3);
	
#endif

	return t;
}

/* ------------------------------------ Tabluated collapse times case ------------------------------------ */

#ifdef TABULATED_CT
#define CT_NBINS_XY (50) 
#define CT_NBINS_D (100) 
#define CT_SQUEEZE (1.2) 
#define CT_EXPO   (1.75)  
#define CT_RANGE_D (7.0)  
#define CT_RANGE_X (3.5)
#define CT_DELTA0 (-1.0)

/* Variables declaration */
int Ncomputations, start, length;
double *CT_table = 0x0;
double bin_x;
double *delta_vector; 
static gsl_interp_accel *accel = 0x0;
gsl_spline ***CT_Spline;

/* Variable declaration in peculiar cases */
#ifdef TRILINEAR
#ifdef HISTO
int *CT_histo = 0x0;
#endif
#endif

#ifdef DEBUG
double ave_diff = 0.0, var_diff = 0.0;
int counter = 0;

#ifdef PRINTJUNK
FILE *JUNK;
#endif
#endif

/* Variable and function declarations */
FILE *CTtableFilePointer;
int check_CTtable_header();
void write_CTtable_header();

/* ------------------------------------- Initialization of the collapse times ------------------------------------- */

int initialize_collapse_times(int ismooth, int onlycompute){	

	/* Variable declarations */
	double l1, l2, l3, del, ampl, x, y, interval, ref_interval, deltaf;

	/* Variables used for loop indices, counters, and flags */
	int	id, ix, iy, i, j, dummy, fail; 

	/* File name */
	char fname[LBLENGTH];

	/* The bin size for delta varies like (delta - CT_DELTA0)^CT_EXPO, 
	so it is smallest around CT_DELTA0 if CT_EXPO>1. The bin size 
	is never smaller than CT_MIN_INTERVAL.

	/*-------------------------- FIRST SMOOTHING RADIUS CASE --------------------------*/ 
	
	if (!ismooth){
		
		/* Dynamically allocation of memory for the delta_vector */
		delta_vector = (double*)malloc(CT_NBINS_D * sizeof(double));

		/* Sets the sampling in density delta */
		/* CT_EXPO value determines the sampling method*/

		if (CT_EXPO == 1){
			/* In this case the spacing is even */
			/* The interval between the elements is 2 times the range of the delta vector divided by the number of bins in the delta vector */
			interval = 2. * CT_RANGE_D / (double)CT_NBINS_D;
			for (int id = 0; id < CT_NBINS_D; id++){
				delta_vector[id] = id * interval - CT_RANGE_D;
			}
		} else {
			/* In this case the sampling is finer around CT_DELTA0 */		
			deltaf = pow(CT_SQUEEZE / CT_EXPO, 1. / (CT_EXPO - 1.));
			if (CT_EXPO == 2){
				ref_interval = ((log((CT_RANGE_D - CT_DELTA0) / deltaf) + log((CT_RANGE_D + CT_DELTA0) / deltaf))/ CT_EXPO + 
							   2. * deltaf / CT_SQUEEZE ) / (CT_NBINS_D - 2.0);
			} else {
				ref_interval = ((pow(CT_RANGE_D - CT_DELTA0, 2. - CT_EXPO) + pow(CT_RANGE_D + CT_DELTA0, 2. - CT_EXPO) 
					           - 2. *pow(deltaf, 2. - CT_EXPO) ) / CT_EXPO / (2. - CT_EXPO) + 2. *deltaf / CT_SQUEEZE ) / (CT_NBINS_D - 2.0);
			}
	 	
			del = -CT_RANGE_D;
			int id = 0;

			do {

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
		if (!ThisTask){
			printf("[%s] Grid for interpolating collapse times: CT_NBINS_D=%d, CT_NBINS_XY=%d\n",fdate(),CT_NBINS_D,CT_NBINS_XY);
		}

		/* Computations are divided among tasks */
		Ncomputations = CT_NBINS_D * CT_NBINS_XY * CT_NBINS_XY;
		bin_x = CT_RANGE_X /(double)(CT_NBINS_XY);

		/* Allocate memory for spline interpolations 
		gsl_interp_accel type is a structure defined in the GSL library that provides an acceleration mechanism for interpolation 
		The function gsl_interp_accel_alloc() dynamically allocates memory for a gsl_interp_accel object and returns a pointer to the allocated memory
		Allocate memory for gsl_interp_accel object */
		accel = gsl_interp_accel_alloc();

		/* Allocate memory for CT_Spline array 
		gsl_spline*** CT_Spline represents a three-dimensional array of gsl_spline objects. The dimensions of this array are CT_NBINS_XY x CT_NBINS_XY x CT_NBINS_D */
      	CT_Spline = (gsl_spline***)calloc(CT_NBINS_XY, sizeof(gsl_spline**));

		for (int i = 0; i < CT_NBINS_XY; ++i){
			CT_Spline[i] = (gsl_spline**)calloc(CT_NBINS_XY, sizeof(gsl_spline*));	
			for (int j = 0; j < CT_NBINS_XY; j++){	
				CT_Spline[i][j] = gsl_spline_alloc(gsl_interp_cspline, CT_NBINS_D);
			}
		}

		/* Allocate memory for CT_table array */
		CT_table = (double*)calloc(Ncomputations, sizeof(double));

#ifdef TRILINEAR
#ifdef HISTO
		/* Allocate memory for CT_histo array */
		CT_histo = (int *)calloc(Ncomputations, sizeof(int));
#endif
#endif

#ifdef DEBUG
#ifdef PRINTJUNK
		/* Create a debug file for writing debug output */
		char filename[50];
		sprintf(filename, "debug.task%d.txt", ThisTask);
		JUNK = fopen(filename, "w");
#endif
#endif

	}

	/*--------------------------  NOT FIRST SMOOTHING RADIUS CASE --------------------------*/
	
	/* Check if collapse times should be read from a file and not in the "onlycompute" mode */
	if (strcmp(params.CTtableFile,"none") && !onlycompute){
		/* In the case where ismooth is false (i.e., the current smoothing iteration is not being performed), Task 0 reads the file */
		if (!ismooth){
			if (!ThisTask){
				/* Performs consistency checks using check_CTtable_header */
				CTtableFilePointer = fopen(params.CTtableFile , "r");
				fail = check_CTtable_header(CTtableFilePointer);
			}

			/* Broadcast the fail status to all tasks to ensure consistency across the parallel execution */
			MPI_Bcast(&fail, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
			if (fail)
				return 1;
		}
		if (!ThisTask){
			/* Read collapse times from the file */
			fread(&dummy, sizeof(int), 1, CTtableFilePointer);
			fread(CT_table, sizeof(double), Ncomputations, CTtableFilePointer);
		}
			
		/* Broadcast CT_table array to all tasks */
		MPI_Bcast(CT_table, Ncomputations, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		/* If the current task is task 0 and ismooth is true for the last smoothing iteration */
		if (!ThisTask && ismooth == Smoothing.Nsmooth - 1)
		fclose(CTtableFilePointer);
	} else {
		/* The collapse time table is computed for each smoothing radius */
		ampl = sqrt(Smoothing.Variance[ismooth]);

		/* Create a map of computations for all tasks by dividing the total number of computations (Ncomputations) among the tasks */
		int *counts   = (int*)malloc(NTasks*sizeof(int));
		int *displs   = (int*)malloc(NTasks*sizeof(int));
		int bunch     = Ncomputations / NTasks;
		int remainder = Ncomputations % NTasks;

		for (int i = 0, ix = 0; i < NTasks; i++){
			counts[i] = bunch + (i < remainder);
			displs[i] = bunch*i + ((remainder > i) ? i : remainder );
		}

		/* Computations for this task */
		for (int i = displs[ThisTask]; i < displs[ThisTask] + counts[ThisTask]; ++i){

			int id = i % CT_NBINS_D;
			ix = (i / CT_NBINS_D) % CT_NBINS_XY;
			iy = i / CT_NBINS_D / CT_NBINS_XY;
			x   = ix * bin_x;
			y   = iy * bin_x;
			l1  = (delta_vector[id] + 2.*x +    y) / 3.0*ampl;
			l2  = (delta_vector[id] -    x +    y) / 3.0*ampl;
			l3  = (delta_vector[id] -    x - 2.*y) / 3.0*ampl;
			CT_table[i]= ell(ismooth,l1,l2,l3);
		}
			
        /* Synchronize computed CT_table across all tasks */
		MPI_Allgatherv( MPI_IN_PLACE, counts[ThisTask], MPI_DOUBLE, CT_table, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD );
            
		/* Free memory */
		free(displs);
		free(counts);

		/* Write the table on a file */			
		if (!ThisTask){
			if (onlycompute){
				strcpy(fname, params.CTtableFile);
			} else {
				sprintf(fname,"pinocchio.%s.CTtable.out",params.RunFlag);
			}
			if (!ismooth){
				/* Create a new file and write header if ismooth is false */
				CTtableFilePointer = fopen(fname,"w");
				write_CTtable_header(CTtableFilePointer);
			} else {
				/* Append to the existing file if ismooth is true */
				CTtableFilePointer = fopen(fname,"a");
			}
      
#ifdef ASCII
            /* Write table contents in ASCII format */
			fprintf(CTtableFilePointer,"################\n# smoothing %2d #\n################\n",ismooth);

			for (i = 0; i < Ncomputations; i++){
				id = i % CT_NBINS_D;
				ix = (i / CT_NBINS_D) % CT_NBINS_XY;
				iy = i / CT_NBINS_D / CT_NBINS_XY;
				x   = ix * bin_x;
				y   = iy * bin_x;
				l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
				l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
				l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

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

	/* Initialize the splines */
	for (int i = 0; i < CT_NBINS_XY; ++i){
		for (int j = 0; j < CT_NBINS_XY; ++j){
			gsl_spline_init(CT_Spline[i][j], delta_vector, &CT_table[i * CT_NBINS_D + j * CT_NBINS_D * CT_NBINS_XY], CT_NBINS_D);
		}
	}

	return 0;
}

/* ------------------------------------- Reset the collapse times ------------------------------------- */

int reset_collapse_times(int ismooth){
	
#ifdef DEBUG
	if (counter)
		{
#ifdef TRILINEAR
		flag = 1;  // Use trilinear interpolation
#endif
#ifdef BILINEAR_SPLINE
		flag = 2;  // Use bilinear spline interpolation
#endif
#ifdef ALL_SPLINE
		flag = 3;  // Use all spline interpolation
#endif
		double aa = ave_diff / (double)counter;
		double std_dev = sqrt(var_diff / (double)counter) - aa * aa;

		/* Print debug information */
		printf("%d   %2d  %10d   %3d %3d %5.2f %5.2f   %20g  %20g DEBUG\n",
			   flag, ismooth, counter, CT_NBINS_D, CT_NBINS_XY, CT_EXPO, CT_SQUEEZE,
			   aa, std_dev);
	}

	/* Reset variables for the next iteration */
	ave_diff = 0.0;
	var_diff = 0.0;
	counter = 0;
#endif

#ifdef TRILINEAR
#ifdef HISTO

	/* Write the histogram of the table */
	/* Create and open the file for writing the histogram */
	char fname[50];
	FILE *fd;

    /* Variable declarations */
	double x, y, l1, l2, l3;
	int i, id, ix, iy;

	/* Set name and open output file */
	sprintf(fname, "popolazione.txt");
	fd = fopen(fname, (ismooth ? "a" : "w"));

	/* Write the header for the smoothing iteration */
	fprintf(fd, "################\n# smoothing %2d #\n################\n", ismooth);

	/* Compute amplitude based on variance */
	double ampl = sqrt(Smoothing.Variance[ismooth]);

	/* Iterate over the computations and write the histogram values */
	for (i = 0; i < Ncomputations; ++i){

		id = i % CT_NBINS_D;
		ix = (i / CT_NBINS_D) % CT_NBINS_XY;
		iy = i / CT_NBINS_D / CT_NBINS_XY;

		x = ix * bin_x;
		y = iy * bin_x;
		l1 = (delta_vector[id] + 2.0 * x + y) / 3.0 * ampl;
		l2 = (delta_vector[id] - x + y) / 3.0 * ampl;
		l3 = (delta_vector[id] - x - 2.0 * y) / 3.0 * ampl;

		fprintf(fd, " %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20d\n",
				ismooth, id, ix, iy, delta_vector[id], x, y, l1, l2, l3, CT_histo[i]);
	}

	/* Close the file */
	fclose(fd);

#endif
#endif

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

inline double interpolate_collapse_time(int ismooth, double l1, double l2, double l3){
	
	/* Variable declarations and assignment */
	double ampl = sqrt(Smoothing.Variance[ismooth]);
	double d    = (l1 + l2 + l3) / ampl;
	double x    = (l1 - l2) / ampl;
	double y    = (l2 - l3) / ampl;
	int ix      = (int)(x / bin_x);
	int iy      = (int)(y / bin_x);

	/* In the unlikely event the point is beyond the limits */
	ix = (ix >= CT_NBINS_XY - 1) ? CT_NBINS_XY - 2 : (ix < 0) ? 0 : ix;
    iy = (iy >= CT_NBINS_XY - 1) ? CT_NBINS_XY - 2 : (iy < 0) ? 0 : iy;

#ifdef ALL_SPLINE
    
	/* Perform bicubic interpolation using a 4x4 grid of control points */
	int ixx, iyy, ixstart, iystart, N=4;
	double xls[4], yls[4], zls[16];

	const gsl_interp2d_type *T = gsl_interp2d_bicubic;
	gsl_interp_accel *xacc = gsl_interp_accel_alloc();
	gsl_interp_accel *yacc = gsl_interp_accel_alloc();

	/* Determine the starting indices for the control point grid */
	ixstart = (ix == 0) ? 0 : (ix >= CT_NBINS_XY - 2) ? CT_NBINS_XY - 4 : ix - 1;
	iystart = (iy == 0) ? 0 : (iy >= CT_NBINS_XY - 2) ? CT_NBINS_XY - 4 : iy - 1;

	/* Create a 4x4 grid of control points for interpolation */
	for (ixx=0; ixx<4; ixx++){
		xls[ixx]=(ixx+ixstart)*bin_x;
		yls[ixx]=(ixx+iystart)*bin_x;
	}
	for (ixx=0; ixx<4; ixx++)
		for (iyy=0; iyy<4; iyy++)
			zls[ixx+iyy*4]=my_spline_eval(CT_Spline[ixx+ixstart][iyy+iystart], d, accel);

	/* Initialize the bicubic interpolation */
	gsl_spline2d *LittleSpline = gsl_spline2d_alloc(T, N, N);
	gsl_spline2d_init(LittleSpline, xls, yls, zls, N, N);

	/* Evaluate the interpolated value using bicubic interpolation */
	double a = gsl_spline2d_eval(LittleSpline, x, y, xacc, yacc);
	gsl_spline2d_free(LittleSpline);

	return a;

#endif

	/* Trilinear interpolation based on the values of d, x, and y. The result is calculated using eight neighboring control points and their corresponding weights */
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
    
	/* Calculates the interpolation coefficients dd, dx, and dy */
	double dd = (d - delta_vector[id]) / (delta_vector[id + 1] - delta_vector[id]);
	double dx = x / bin_x - ix;
	double dy = y / bin_x - iy;
    
	/* Perform trilinear interpolation using the eight neighboring control points*/
	return (((1. - dd) * (1. - dx) * (1. - dy) * CT_table[id + (ix) * CT_NBINS_D + (iy) * CT_NBINS_D * CT_NBINS_XY]) +
           ((dd) * (1. - dx) * (1. - dy) * CT_table[(id + 1) + (ix) * CT_NBINS_D + (iy) * CT_NBINS_D * CT_NBINS_XY]) +
    	   ((1. - dd) * (dx) * (1. - dy) * CT_table[id + (ix + 1) * CT_NBINS_D + (iy) * CT_NBINS_D * CT_NBINS_XY])   +
           ((dd) * (dx) * (1. - dy) * CT_table[(id + 1) + (ix + 1) * CT_NBINS_D + (iy) * CT_NBINS_D * CT_NBINS_XY])  +
           ((1. - dd) * (1. - dx) * (dy) * CT_table[id + (ix) * CT_NBINS_D + (iy + 1) * CT_NBINS_D * CT_NBINS_XY])   +
           ((dd) * (1. - dx) * (dy) * CT_table[(id + 1) + (ix) * CT_NBINS_D + (iy + 1) * CT_NBINS_D * CT_NBINS_XY])  +
           ((1. - dd) * (dx) * (dy) * CT_table[id + (ix + 1) * CT_NBINS_D + (iy + 1) * CT_NBINS_D * CT_NBINS_XY])    +
           ((dd) * (dx) * (dy) * CT_table[(id + 1) + (ix + 1) * CT_NBINS_D + (iy + 1) * CT_NBINS_D * CT_NBINS_XY]))  ;

#endif
#ifdef BILINEAR_SPLINE
    
	/* Calculate the offsets */
	double dx = x / bin_x-ix;
	double dy = y / bin_x-iy;

    /* Perform bilinear interpolation using the four neighboring control points */
	return ((1. - dx) * (1. - dy) * my_spline_eval(CT_Spline[ix][iy], d, accel) +
		    (dx) * (1. - dy) * my_spline_eval(CT_Spline[ix + 1][iy], d, accel)  +
		    (1. - dx) * (dy) * my_spline_eval(CT_Spline[ix][iy + 1], d, accel)  +
		    (dx) * (dy) * my_spline_eval(CT_Spline[ix + 1][iy + 1], d, accel))  ;
#endif
}

/*-------------------------- Consistency check of the header of CT Table file with the current run --------------------------------------------------------*/

int check_CTtable_header(){

	/* Checks that the information contained in the header is consistent with the run 
	NB: this is done by Task 0 */

	int fail = 0;
    int dummy;
    double fdummy;

	/* Check the CT table type based on compilation flags */
    fread(&dummy, sizeof(int), 1, CTtableFilePointer);
	
#ifdef ELL_CLASSIC 
	if (dummy != 1){   /* classic ellipsoidal collapse, tabulated */	
        printf("ERROR: CT table not constructed for ELL_CLASSIC, %d\n", dummy);
        fail = 1;
	}
#endif
#ifdef ELL_SNG
#ifndef MOD_GRAV_FR
	if (dummy != 3){   /* standard gravity, numerical ellipsoidal collapse */
		printf("ERROR: CT table not constructed for ELL_SNG and standard gravity, %d\n", dummy);
        fail = 1;
	}
#else
	if (dummy != 4){   /* f(R) gravity, numerical ellipsoidal collapse */
		printf("ERROR: CT table not constructed for ELL_SNG and MOD_GRAV_FR, %d\n", dummy);
        fail = 1;	
	}
#endif
#endif

	/* Check the CT table parameters */
	fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
	if (fabs(fdummy - params.Omega0) > 1.e-10){
		printf("ERROR: CT table constructed for the wrong Omega0, %f in place of %f\n", fdummy, params.Omega0);
        fail = 1;
	}
    fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
    if (fabs(fdummy - params.OmegaLambda) > 1.e-10){
        printf("ERROR: CT table constructed for the wrong OmegaLambda, %f in place of %f\n", fdummy, params.OmegaLambda);
        fail = 1;
    }
	fread(&fdummy, sizeof(double), 1, CTtableFilePointer);
    if (fabs(fdummy - params.Hubble100) > 1.e-10){
        printf("ERROR: CT table constructed for the wrong Hubble100, %f in place of %f\n", fdummy, params.Hubble100);
        fail = 1;
    }

	/* Check the CT table size and sampling */
    fread(&dummy, sizeof(int), 1, CTtableFilePointer);
    if (dummy != Ncomputations){
        printf("ERROR: CT table has the wrong size, %d in place of %d\n", dummy, Ncomputations);
        fail = 1;
    }
    fread(&dummy, sizeof(int), 1, CTtableFilePointer);
    if (dummy != CT_NBINS_D){
        printf("ERROR: CT table has the wrong density sampling, %d in place of %d\n", dummy, CT_NBINS_D);
        fail = 1;
    }
    fread(&dummy, sizeof(int), 1, CTtableFilePointer);
    if (dummy != CT_NBINS_XY){
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
#endif 

/* ------------------------------------  Tabulated collapse time : END ------------------------------------ */

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

/* Orders a,b,c in decreasing order a>b>c */
inline void ord(double * restrict a, double * restrict b, double * restrict c){

	double lo , hi;
	hi = max( max(*a, *b), *c );
	lo = min( min(*a, *b), *c );
	*b = *a + *b+ *c - lo - hi;
	*a = hi;
	*c = lo;
}

