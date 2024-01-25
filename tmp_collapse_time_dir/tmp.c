
inline int compute_collapse_times(int ismooth) {    

	 /*---------------------DEBUG---------------------------------------*/

#ifdef DEBUG

    // Open a file for writing collapse times

    FILE *collapseTimeFile = fopen("collapse_times.txt", "w");
    if (collapseTimeFile == NULL) {
        printf("Error opening file for writing.\n");
        return 1;
    }

#endif

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
		#pragma acc kernels
		{ 
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
		}	
	/* Declaration and initialization of variables when the openmp option is ON	
	Inside the parallel region, each thread initializes its own local variables i.e.
	mylocal_average, mylocal_variance, cputime_invcoll, and cputime_ell 
	These variables are private to each thread and are not shared among threads, ensuring data consistency.*/
     
#if defined( _OPENACC )

	double invcoll_update = 0, ell_update = 0;

	#pragma acc kernels
	{

	/* Thread-specific variables declaration */	
	double  mylocal_average, mylocal_variance;
    double  cputime_invcoll;
		
	/* Support variables */	
	int     fails = 0;
	// int     pid = omp_get_thread_num();
		
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
// #ifdef _OPENMP
// #pragma omp for nowait
// #endif

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

#ifdef DEBUG
		
		fprintf(collapseTimeFile, "%d %d %f %f %f %g\n", index, ismooth, lambda1, lambda2, lambda3, Fnew);

#endif

		/*-----------------------------------------------------------------*/

		/* Fail check during computation of the collapse time */
		if (fail){
			printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
			// fflush(stdout);

/* Workflow control */
#if !defined( _OPENACC )		    		    
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

// #if defined( _OPENMP )

	#pragma acc atomic
	local_variance += mylocal_variance;

	#pragma acc atomic
	local_average += mylocal_average;  

	/* ------------------------------------ Updates total CPU collapse time -----------------------------------*/
	/* cputime_invcoll and cputime_ell are local variables that are updated by each thread and then accumulated using atomic operations */

	#pragma acc atomic
	invcoll_update += cputime_invcoll;

	#pragma acc atomic	
	ell_update     += cputime_ell;

// #endif	

#if defined( _OPENACC ) 
	}
#endif

	int all_fails = 0; 

	/* -------------------------- Updates cpu collapse time for a single threads -----------------------------*/

// #if defined (_OPENMP)

	#pragma acc atomic
	cputime.invcoll += invcoll_update/1; 

	#pragma omp atomic	
	cputime.ell     += ell_update/1;

	/* fails is a local variable that indicates the number of failures in the computation. 
	Each thread updates fails, and then the values are accumulated using atomic operations */

	#pragma acc atomic
	all_fails       += fails;

// #endif
    
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