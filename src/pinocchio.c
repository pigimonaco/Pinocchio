/*****************************************************************
 *                        PINOCCHIO  V4.1                        *
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
#include "def_splines.h"

void abort_code(void);
void write_cputimes(void);
void greetings(void);

int main(int argc, char **argv, char **envp)
{
  double time;
  int ThisGrid;


#ifdef USE_GPERFTOOLS
  ProfilerStart("pinocchio_gprofile.log");
#endif
  
  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);

  /* timing of the code */
  cputime.total=MPI_Wtime();
  greetings();

#ifdef _OPENMP
  /* initialization of OpemMP */
#pragma omp parallel
  {
#pragma omp master
    internal_data.nthreads_omp = omp_get_num_threads();
  }
#endif

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
    for (ThisGrid=0; ThisGrid<Ngrids; ThisGrid++)
      {
	write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	time=reverse_transform(ThisGrid);
	if (!ThisTask)
	  printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);
	write_from_rvector(ThisGrid, density[ThisGrid]);
	if (write_density(ThisGrid))
	  abort_code();
      }
    
  /* called as "pinocchio.x parameterfile 1" it writes the density field of Grid 0
     in configuration space and exits */
  if (argc>=3 && atoi(argv[2])==1)
    {
      if (!ThisTask)
	{
	  printf("In this configuration pinocchio only writes the linear density field\n");
	  write_cputimes();
	  printf("Pinocchio done!\n");
	}

      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif
      return 0;
    }

/*   /\* called as "pinocchio.x parameterfile 2" it computes and writes displacement fields, then exit *\/ */
/*   if (argc>=3 && atoi(argv[2])==2) */
/*     { */
/*       if (!ThisTask) */
/* 	{ */
/* 	  printf("In this configuration pinocchio will only produce a GADGET snapshot\n"); */
/* 	  printf("at the first redshift specified by the %s file (z=%f)\n", */
/* 		 params.OutputList,outputs.z[0]); */
/* 	} */


/* #ifdef THREE_LPT */
/*       if (!ThisTask) */
/* 	{ */

/* 	  // METTERE 3LPT!!!! */

/* 	  printf("*********************************************************\n"); */
/* 	  printf("Writing of LPT snapshot presently does not work with 3LPT\n"); */
/* 	  printf("Please recompile without the THREE_LPT directive\n"); */
/* 	  printf("*********************************************************\n"); */
/* 	} */
/* #else */

/*       if (compute_displacements()) */
/*         abort_code(); */

/*       /\* if (write_LPT_snapshot(outputs.z[0])) *\/  // RIMETTERE!!! */
/*       /\*   abort_code(); *\/ */

/* #endif */
	
/*       if (!ThisTask) */
/* 	printf("Pinocchio done!\n"); */
/*       MPI_Finalize(); */
/* #ifdef USE_GPERFTOOLS */
/*       ProfilerStop(); */
/* #endif */

/*       return 0; */
/*     } */

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

      /* on request, write products on a snapshot */
      if (params.WriteProducts)
	if (write_products())
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
  printf("Total:            %14.6f\n", cputime.total);
  printf("Initialization:   %14.6f (%5.2f%%)\n", cputime.init, 100.*cputime.init/cputime.total);
  printf("  Density in PS:  %14.6f (%5.2f%%)\n", cputime.dens, 100.*cputime.dens/cputime.total);
  printf("fmax:             %14.6f (%5.2f%%)\n", cputime.fmax, 100.*cputime.fmax /cputime.total);
#ifdef TWO_LPT
  printf("  LPT:            %14.6f (%5.2f%%)\n", cputime.lpt,  100.*cputime.lpt  /cputime.total);
#endif
  printf("  Derivatives:    %14.6f (%5.2f%%)\n", cputime.deriv,  100.*cputime.deriv  /cputime.total);
  printf("    Mem transfer: %14.6f (%5.2f%%)\n", cputime.mem_transf, 100.*cputime.mem_transf  /cputime.total);
  printf("    FFTs:         %14.6f (%5.2f%%)\n", cputime.fft,  100.*cputime.fft  /cputime.total);
  printf("  Collapse times: %14.6f (%5.2f%%)\n", cputime.coll, 100.*cputime.coll /cputime.total);
  printf("    inv.collapse: %14.6f (%5.2f%%)\n", cputime.invcoll, 100.*cputime.invcoll /cputime.total);
  printf("    ellipsoid:    %14.6f (%5.2f%%)\n", cputime.ell, 100.*cputime.ell /cputime.total);
  printf("  Velocities:     %14.6f (%5.2f%%)\n", cputime.vel,  100.*cputime.vel  /cputime.total);
  printf("Fragmentation:    %14.6f (%5.2f%%)\n", cputime.frag, 100.*cputime.frag /cputime.total);
  printf("  Redistribution: %14.6f (%5.2f%%)\n", cputime.distr,100.*cputime.distr/cputime.total);
  printf("  Sorting:        %14.6f (%5.2f%%)\n", cputime.sort, 100.*cputime.sort /cputime.total);
#ifdef PLC
  printf("  Groups total:   %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
  printf("  Groups PLC:     %14.6f (%5.2f%%)\n", cputime.plc,100.*cputime.plc/cputime.total);
#else
  printf("  Groups:         %14.6f (%5.2f%%)\n", cputime.group,100.*cputime.group/cputime.total);
#endif
  printf("Total I/O:        %14.6f (%5.2f%%)\n", cputime.io,   100.*cputime.io   /cputime.total);
}

void greetings(void)
{
  /* This is a list of messages to declare the most relevant precompiler directives in the stdout */

  if (!ThisTask)
    {
      printf("[%s] This is pinocchio V4.XX, running on %d MPI tasks\n\n",fdate(),NTasks);
#ifdef _OPENMP
      printf( "Using %d OpenMP threads\n", internal_data.nthreads_omp );
#endif

#ifdef USE_FFT_THREADS
      printf( "Using threaded-FFTs\n");
#endif
#ifdef TWO_LPT
#ifndef THREE_LPT
      printf("This version uses 2LPT displacements\n");
#else
      printf("This version uses 3LPT displacements\n");
#endif
#else
      printf("This version uses Zeldovich displacements\n");
#endif
#ifdef ROTATE_BOX
      printf("The output will be rotated to reproduce N-GenIC and 2LPTic orientation\n");
#endif
#ifdef NORADIATION
      printf("Radiation is not included in the Friedmann equations\n");
#else
      printf("Radiation is included in the Friedmann equations\n");
#endif
#ifdef TIMELESS_SNAPSHOT
      printf("Production of the timeless snapshot has been activated\n");
#endif
#ifdef LIGHT_OUTPUT
      printf("Catalogs will be written in the light version\n");
#endif

#ifdef TABULATED_CT
#ifdef ELL_CLASSIC
      printf("Ellipsoidal collapse will be tabulated as Monaco (1995)\n");
#endif
#ifdef ELL_SNG
      printf("Numerical integration of ellipsoidal collapse will be tabulated\n");
#endif
#else
      printf("Ellipsoidal collapse will be computed as Monaco (1995)\n");
#endif

#ifdef WHITENOISE
      printf("Initial conditions will be read from a white noise file\n");
#endif

#ifdef SCALE_DEPENDENT
      printf("This version of the code works with scale-dependent growing modes;\n");
#ifdef MOD_GRAV_FR
      printf("Scales will range from %10g to %10g 1/Mpc, in %d steps\n",
	     0.0, pow(10.,LOGKMIN+(NkBINS-1)*DELTALOGK), NkBINS);
      printf("Gravity will be given by Hu-Sawicki f(R) with f_R0=%7g\n",FR0);
#endif
#ifdef READ_PK_TABLE
      printf("Scales will range from %10g to %10g 1/Mpc, in %d steps\n",
	     pow(10.,LOGKMIN), pow(10.,LOGKMIN+(NkBINS-1)*DELTALOGK), NkBINS);
      printf("Scale-dependent growth rates will be worked out from CAMB P(k) files\n");
#ifdef ONLY_MATTER_POWER
      printf("The power spectrum will include only dark matter + baryon fluctuations, excluding neutrinos (if present)\n");
#else
      printf("The power spectrum will include TOTAL matter fluctuations, including neutrinos (if present)\n");
#endif
#endif
#endif

#ifdef NO_RANDOM_MODULES
      printf("Initial conditions will be generated with non-random modules of the Fourier modes\n");
#endif

      printf("\n");

    }
}




  /* //LEVARE */
  /* write_in_cvector(0, kdensity[0]); */

  /* dump_cvector((double*)kdensity[0], params.use_transposed_fft*internal.tasks_subdivision_dim, */
  /* 	       MyGrids[0].GSglobal[_x_], */
  /* 	       MyGrids[0].GSlocal_k, */
  /* 	       MyGrids[0].GSstart_k, "kdensity.check", 0); */


  /* time=reverse_transform(0); */
  /* write_from_rvector(0, density[0]); */

  /* dump_rvector((double*)rvector_fft[0], */
  /* 	       MyGrids[0].GSglobal[_x_], */
  /* 	       MyGrids[0].GSlocal, */
  /* 	       MyGrids[0].GSstart, "density.check", 0); */

  /* fflush(stdout); */
  /* MPI_Barrier(MPI_COMM_WORLD); */
