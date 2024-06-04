#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif



int main() 
{
  
  int num_devices = omp_get_num_devices();
  printf("Number of available devices %d\n", num_devices);
   
  double test_time = 0;
  int nteams;
  int fails = 0 ;
  int nthreads;
  int grid_size = 10;
  
  /* Local average and variance declaration */
	double local_average  = 0;
  double local_variance = 0;

#ifdef _OPENMP

  double invcoll_update = 0;
  
  
  #pragma omp target teams map(tofrom: nteams, nthreads, invcoll_update, local_variance, local_average) map(to:grid_size)  num_teams(2)  reduction(+:fails)
  {
    
    /* Thread-specific variables declaration */	
    double cputime_invcoll = 0;
		double mylocal_average  = 0;
    double mylocal_variance = 0;

#else 

#define cputime_invcoll  test_time

#endif
    
    double tmp = omp_get_wtime();

    if (omp_is_initial_device()) {
        printf("Running on host\n");    
      } else {
        nteams= omp_get_num_teams(); 
        nthreads= omp_get_num_threads();
        
        for (int index = 0; index < grid_size; index ++){
          double diff_ten[6];
          for (int i = 0; i < 6 ; i++){
            diff_ten[i] = 1;
          }
        
        double delta      = diff_ten[0] + diff_ten[1] + diff_ten[2];
    	  mylocal_average  += delta;
    	  mylocal_variance += delta*delta;

/* Workflow control */
#if !defined( _OPENMP )		    		    
    		return 1;
#else
       	fails = 1;
#endif

        }
        // printf("Running on device with %d teams in total and %d threads in each team\n",nteams,nthreads);
      }

      cputime_invcoll += omp_get_wtime() - tmp;

#if defined( _OPENMP )

      #pragma omp atomic
	    invcoll_update += cputime_invcoll;
      
      #pragma omp atomic
	    local_variance += mylocal_variance;

	    #pragma omp atomic
	    local_average += mylocal_average;  

#endif	
      
      printf("mylocal_average:  %lf\n", mylocal_average);
      printf("mylocal_variance: %lf\n", mylocal_variance);
      printf("cputime_invcoll:  %lf\n", cputime_invcoll);
      printf("fails:            %d\n",  fails);

#if defined( _OPENMP ) 
	}
#endif

  int all_fails = 0;

#if defined (_OPENMP)

	#pragma omp atomic
	test_time += invcoll_update/nteams; 

  // #pragma omp atomic
	all_fails += fails;

#endif
   
  printf("local_average:  %lf\n", local_average);
  printf("local_variance: %lf\n", local_variance);
  printf("test_time:      %lf\n", test_time);
  printf("all_fails:      %d\n",  all_fails);
  
  //  printf("Running on device with %d teams in total and %d threads in each team\n",nteams,nthreads);

}