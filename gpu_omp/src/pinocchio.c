/* ######HEADER###### */

#include "pinocchio.h"
#include "def_splines.h"

void abort_code(void);
void write_cputimes(void);

int main(int argc, char **argv, char **envp)
{

#ifdef USE_GPERFTOOLS
  ProfilerStart("pinocchio_gprofile.log");
#endif
  
  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level); // Hybrid MPI and OPENMP parallel
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask); // gives you your rank i.e. ID
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks); // size of your pool 

#ifdef _OPENMP
  /* initialization of OpemMP */
#pragma omp parallel
  {
#pragma omp master
    internal.nthreads_omp = omp_get_num_threads();
  }
#endif

  /* GPU initialization */
#ifdef GPU_OMP
  if (initialization_gpu_omp())
    abort_code();
#endif // GPU_OMP
  
  /* timing of the code */
  cputime.total=MPI_Wtime();
  greetings();

  /* checks that the parameter file is given in the command line */
  if (argc<2)
    {
      if (!ThisTask)
	printf("Usage: pinocchio.x parameterfile\n");
      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif
      return 0;
    }

  /* initialization */
  memset(&params, 0, sizeof(param_data));
  strcpy(params.ParameterFile,argv[1]);
  if (initialization())
    abort_code();


  /*****************************************/
  /*********** Special behaviour ***********/
  /*****************************************/
  /* On request, it writes the density field */

  if (params.WriteDensity || (argc>=3 && atoi(argv[2])==1) )
    {
#ifdef SNAPSHOT
      /* called as "pinocchio.x parameterfile 1" it writes the density field
	 in configuration space and exits */
      if (argc>=3 && atoi(argv[2])==1 && !ThisTask)
	printf("In this configuration pinocchio only writes the linear density field\n");

      for (int ThisGrid=0; ThisGrid<Ngrids; ThisGrid++)
	{
	  write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	  double time=reverse_transform(ThisGrid);
	  if (!ThisTask)
	    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);
	  write_from_rvector(ThisGrid, density[ThisGrid]);
	  if (write_density(ThisGrid))
	    abort_code();
	}
#else
      if (!ThisTask)
	printf("Please recompile the code with SNAPSHOT directive to write the density\n");
#endif
      if (!ThisTask)
	{
	  write_cputimes();
	  printf("Pinocchio done!\n");
	}

      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif
      return 0;
    }

  /* called as "pinocchio.x parameterfile 3" it computes and writes displacement fields, then exit */

  if (argc>=3 && atoi(argv[2])==3)
    {
      if (!ThisTask)
	      {
	        printf("In this configuration pinocchio will only produce a GADGET snapshot\n");
	        printf("at the first redshift specified by the %s file (z=%f)\n", params.OutputList,outputs.z[0]);
        }

#ifdef THREE_LPT
      if (!ThisTask)
	      {
	        printf("*********************************************************\n");
	        printf("Writing of LPT snapshot presently does not work with 3LPT\n");
	        printf("Please recompile without the THREE_LPT directive\n");
	        printf("*********************************************************\n");
	      }

#else

      if (compute_displacements())
          abort_code();

      if (write_LPT_snapshot(outputs.z[0]))
          abort_code();

#endif

      if (!ThisTask) 
        printf("Pinocchio done!\n");
    
      MPI_Finalize();

 #ifdef USE_GPERFTOOLS 
      ProfilerStop();

 #endif 

      return 0; 
    } 


  /* called as "pinocchio.x parameterfile 2" it computes and writes collapse time table, then exit */
  if (argc>=3 && atoi(argv[2])==3)
    {
#ifdef TABULATED_CT
      if (!ThisTask)
	{
	  printf("In this configuration pinocchio will only compute a table of collapse times\n");
	}

      /* CYCLE ON SMOOTHING RADII */
      for (int ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
	{
	  double cputmp=MPI_Wtime();

	  if (!ThisTask)
	    printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n", 
		   fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
		   sqrt(Smoothing.Variance[ismooth]) );

	  if (initialize_collapse_times(ismooth,1))
	    return 1;

	  if (!ThisTask)
	    printf("[%s] Collapse times computed, cpu time =%f s\n",fdate(),cputmp);

	}

      if (!ThisTask)
	printf("Pinocchio done!\n");
      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif

      return 0;

#else
      if (!ThisTask)
	{
	  printf("Sorry but you have to compile the code with TABULATED_CT to compute CT table\n");
	}
      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif

      return 0;
#endif
    }

  /******************************************/
  /*********** Standard behaviour ***********/
  /******************************************/


  if (params.ReadProductsFromDumps)
    {
      /* on request, read products from dump files and skip the first part */
      if (read_dumps())
	abort_code();
    }
  else
    {
      /* computation of collapse times and displacements */
      if (compute_fmax())
	abort_code();

      /* on request, dump products to files for skipping fmax */
      if (params.DumpProducts)
	if (dump_products())
	  abort_code();
    }

  /* fragmentation of the collapsed medium */
  if (fragment_driver())
    abort_code();

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* output detailed cpu times */
  cputime.total=MPI_Wtime()-cputime.total;
  if (!ThisTask)
    write_cputimes();

  /* done */
  if (!ThisTask)
    printf("Pinocchio done!\n");
  MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif

  return 0;
}


void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}


void write_cputimes()
{
  printf("Total               : %14.6f\n", cputime.total);
  printf("Initialization      : %14.6f (%5.2f%%)\n", cputime.init, 100.*cputime.init/cputime.total);
  printf("  Density in PS     : %14.6f (%5.2f%%)\n", cputime.dens, 100.*cputime.dens/cputime.total);
  printf("fmax                : %14.6f (%5.2f%%)\n", cputime.fmax, 100.*cputime.fmax /cputime.total);
#ifdef TWO_LPT
  printf("  LPT               : %14.6f (%5.2f%%)\n", cputime.lpt,  100.*cputime.lpt  /cputime.total);
#endif
  printf("  Derivatives       : %14.6f (%5.2f%%)\n", cputime.deriv,  100.*cputime.deriv  /cputime.total);
  printf("    Mem transfer    : %14.6f (%5.2f%%)\n", cputime.mem_transf, 100.*cputime.mem_transf  /cputime.total);
  printf("    FFTs            : %14.6f (%5.2f%%)\n", cputime.fft,  100.*cputime.fft  /cputime.total);
  printf("  Collapse times    : %14.6f (%5.2f%%)\n", cputime.coll, 100.*cputime.coll /cputime.total);
#ifdef GPU_OMP
  printf("    GPU computation : %14.6f (%5.2f%%)\n", gputime.computation.collapse_times, 100.*gputime.computation.collapse_times/cputime.total);
  printf("    GPU mem transfer: %14.6f (%5.2f%%)\n", gputime.memory_transfer.collapse_times, 100.*gputime.memory_transfer.collapse_times/cputime.total);
#endif // GPU_OMP  
  printf("  Velocities        : %14.6f (%5.2f%%)\n", cputime.vel,  100.*cputime.vel  /cputime.total);
  printf("Fragmentation       : %14.6f (%5.2f%%)\n", cputime.frag, 100.*cputime.frag /cputime.total);
  printf("  Redistribution    : %14.6f (%5.2f%%)\n", cputime.distr,100.*cputime.distr/cputime.total);
  printf("  Sorting           : %14.6f (%5.2f%%)\n", cputime.sort, 100.*cputime.sort /cputime.total);
#ifdef PLC
  printf("  Groups total      : %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
  printf("  Groups PLC        : %14.6f (%5.2f%%)\n", cputime.plc,100.*cputime.plc/cputime.total);
#else
  printf("  Groups            : %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
#endif
  printf("Total I/O           : %14.6f (%5.2f%%)\n", cputime.io,   100.*cputime.io   /cputime.total);
}

