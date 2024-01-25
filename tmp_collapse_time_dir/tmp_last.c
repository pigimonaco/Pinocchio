inline int compute_collapse_times(int ismooth) {    

/*--------------------- GPU CALCULATION OF COLLAPSE TIME ----------------------------*/

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


	if (!ismooth)

    /* OpenMP target directive to offload the initialization part to the GPU */
	/* Transfering data to GPU before initializing it. In this way, the data will be already there when needed in the offloaded region */
	#pragma omp target data map(tofrom: products[0:MyGrids[0].total_local_size])
	{
    	#pragma omp target teams distribute parallel for
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

/* Declaration and initialization of variables when the OpenMP option is ON
Inside the parallel region, each thread initializes its own local variables, i.e., mylocal_average, mylocal_variance, cputime_invcoll, and cputime_ell
These variables are private to each thread and are not shared among threads, ensuring data consistency. */

#if defined( _OPENMP )

	double invcoll_update = 0, ell_update = 0;
	int num_teams, team_size;

	#pragma omp target teams map(tofrom: local_average, local_variance, products[0:MyGrids[0].total_local_size]) map(to:second_derivatives[0][6][MyGrids[0].total_local_size]) map(from: num_teams, team_size)
	{

	/* Thread-specific variables declaration */	
	double  mylocal_average, mylocal_variance;
    double  cputime_invcoll;
		
	/* Support variables */	
	int     fails = 0;

	num_teams = omp_get_num_teams();
    // team_size = omp_get_team_size();

	/* Initializitation */
	mylocal_average     = 0;
	mylocal_variance    = 0;
	cputime_invcoll     = 0;
	cputime_ell         = 0;

#else

/* If OpenMP is OFF 
This approach is used to avoid conditional compilation directives within the code and maintain consistency in variable
names across the parallel and serial versions of the code. 
When OpenMP is disabled, the program uses the same variable names, but they refer to the original variables directly */

#define mylocal_average   local_average
#define mylocal_variance  local_variance
#define cputime_invcoll   cputime.invcoll
#define cputime_ell       cputime.ell

#endif

	/*----------------- Calculation of variance, average, and collapse time ----------------*/

	double tmp = omp_get_wtime();

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

    		// fprintf(collapseTimeFile, "%d %d %f %f %f %g\n", index, ismooth, lambda1, lambda2, lambda3, Fnew);

    		/*-----------------------------------------------------------------*/

    		/* Fail check during computation of the collapse time */
    		if (fail){
        		printf("ERROR on task %d: failure in inverse_collapse_time\n", ThisTask);
        		// fflush(stdout);
  

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
	cputime.invcoll += omp_get_wtime() - tmp;

	/* Update CPU variables with GPU results */

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
	cputime.invcoll += invcoll_update/num_teams; 

	#pragma omp atomic	
	cputime.ell     += ell_update/num_teams;

	/* fails is a local variable that indicates the number of failures in the computation. 
	Each thread updates fails, and then the values are accumulated using atomic operations */

	#pragma omp atomic
	all_fails       += fails;

#endif

/*-------------------- END OF GPU CALCULATION -----------------------------------------------*/

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
