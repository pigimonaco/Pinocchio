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

#ifdef PRECISE_TIMING

#define SET_WTIME cputime.partial = MPI_Wtime();
#define ASSIGN_WTIME(INIT, ACC) do { double ttt= MPI_Wtime(); cputime.ACC = ttt - cputime.INIT; } while(0)
#define ACCUMULATE_WTIME(INIT, ACC) do { double ttt= MPI_Wtime(); cputime.ACC += ttt - cputime.INIT; } while(0)
#else
#define SET_WTIME
#define ASSIGN_WTIME(INIT, ACC)
#define ACCUMULATE_WTIME(INIT, ACC)
#endif

int          init_cosmology(void);
int          set_smoothing(void);
int          generate_densities(void);
int          set_subboxes(void);
int          set_plc(void);
unsigned int gcd(unsigned int, unsigned int);
int          set_fft_decomposition(void);
/* int get_task_decomposition(int, int**); */
/* int integer_factorization(UL, UL **); */


int initialization()
{
  
  // timing
  cputime.init = MPI_Wtime();

  // this is for gsl integration
  workspace = gsl_integration_workspace_alloc(NWINT);
  
  // this is the initialization of the random number generator
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  // reading of parameters from file and few other parameter initializations
  if (set_parameters())
    return 1;

#ifdef USE_FFT_THREADS
  if ( internal.nthreads > 1 )
    dprintf(VMSG, 0, "Using %d threads for FFTs\n", internal.nthreads );
#endif
  
  /* Initialize pfft */
  pfft_init();

  /* Inititalize fftw */
#ifdef USE_FFT_THREADS
    fftw_init_threads();
#endif

  fftw_mpi_init();    

  
  if(set_fft_decomposition())
    return 1;

  if(ThisTask == 0)
    dprintf(VMSG, ThisTask, "cube subdivision [%d dim]: %d x %d x %d = %d processes\n",
	    internal.tasks_subdivision_dim,
	    internal.tasks_subdivision_3D[0],
	    internal.tasks_subdivision_3D[1],
	    internal.tasks_subdivision_3D[2],
	    internal.tasks_subdivision_3D[0] *
	    internal.tasks_subdivision_3D[1] *
	    internal.tasks_subdivision_3D[2]);
  
  if( pfft_create_procmesh(internal.tasks_subdivision_dim, MPI_COMM_WORLD, internal.tasks_subdivision_3D, &FFT_Comm) )
    {
      int all = 1;
      for(int iii = 0; iii < internal.tasks_subdivision_dim; iii++)
  	all *= internal.tasks_subdivision_3D[iii];
      
      pfft_fprintf(MPI_COMM_WORLD, stderr, "Error while creating communicator and mesh with %d processes\n", all);
      return 1;
    }

#ifdef SCALE_DEPENDENT_GROWTH
  // reads the P(k) tables from CAMB
  if (read_power_table_from_CAMB())
    return 1;
#endif

  // call to cosmo.c: initialization of cosmological functions
  if (initialize_cosmology())
    return 1;

  // set the properties of grids and initialize FFTW quantities, including vectors
  if (set_grids())
    return 1;

  // computes the smoothing radii
  if (set_smoothing())
    return 1;

  // now it re-initializes the variance with a top-hat filter
  WindowFunctionType=2;
  if (initialize_MassVariance())
    return 1;

#ifdef SCALE_DEPENDENT_GROWTH
  // initialize scale-dependent growing modes and fomegas
  if (initialize_ScaleDependentGrowth())
    return 1;
#endif

  SET_WTIME;
  // computes the number of sub-boxes for fragmentation
  if (set_subboxes())
    {
      if (!ThisTask)
	printf("Pinocchio done!\n");
      MPI_Finalize();
      exit (0);
    }
  ASSIGN_WTIME(partial, set_subboxes);

  SET_WTIME;
  // initializes quantities needed for the on-the-fly reconstruction of PLC
  if (set_plc())
    return 1;
  ASSIGN_WTIME(partial, set_plc);

  // this barrier is set to have correct stdout in case the run cannot start
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);


  SET_WTIME;
  // allocations of memory for fmax and memory tests 
  if (allocate_main_memory())
    return 1;
  ASSIGN_WTIME(partial, memory_allocation);

  
  SET_WTIME;
  // initialization of fft plans
  if (initialize_fft())
    return 1;
  ASSIGN_WTIME(partial, fft_initialization);

  
  // generation of initial density field
  if (generate_densities())
    return 1;

  cputime.init = MPI_Wtime()-cputime.init;

  if (!ThisTask)
    {      
      dprintf(VMSG, ThisTask, "[%s] initialization done, initialization cpu time = %14.6f\n", fdate(), cputime.init);
      dprintf(VMSG, ThisTask, "\t\t set subboxes time = %14.6f\n"
	      "\t\t set plc time = %14.6f\n"
	      "\t\t memory allocation time = %14.6f\n"
	      "\t\t fft initialization time = %14.6f\n"
	      "\t\t density generation time = %14.6f\n",
	      cputime.set_subboxes, cputime.set_plc, cputime.memory_allocation,
	      cputime.fft_initialization, cputime.dens);
    }

  return 0;
}


int set_parameters()
{
  int i;

  // set default internal parameters
  internal.verbose_level                   = VDIAG;
  internal.dump_seedplane                  = 0;
  internal.dump_kdensity                   = 0;
  internal.large_plane                     = 1;
  internal.mimic_original_seedtable        = 0;
  internal.dump_vectors                    = 0;
  internal.constrain_task_decomposition[0] = 0;
  internal.constrain_task_decomposition[1] = 0;
  internal.constrain_task_decomposition[2] = 0;
  internal.tasks_subdivision_3D[0]         = 0;
  internal.tasks_subdivision_3D[1]         = 0;
  internal.tasks_subdivision_3D[2]         = 0;
  internal.nthreads                        = 1;
  
  if(read_parameter_file())
    return 1;

  if (params.BoxInH100)
    {
      params.BoxSize_h100  = params.BoxSize;
      params.BoxSize_htrue = params.BoxSize/params.Hubble100;
    }
  else
    {
      params.BoxSize_h100  = params.BoxSize*params.Hubble100;
      params.BoxSize_htrue = params.BoxSize;
    }
  params.InterPartDist = params.BoxSize_htrue / params.GridSize[0];

  params.ParticleMass = 2.775499745e11 *
    params.Hubble100 * params.Hubble100 * params.Omega0 *
    pow(params.InterPartDist, 3.0);
  
  strcpy(params.DataDir,"Data/");

  if (!params.NumFiles)
    params.NumFiles=1;

  /* The number of files must be a divisor of the number of tasks */
  if (NTasks%params.NumFiles != 0)
    {
      while (NTasks%params.NumFiles != 0)
	params.NumFiles--;

      if (!ThisTask)
	printf("Warning: NumFiles must be a divisor of NTasks, it has been fixed to %d\n",
	       params.NumFiles);
    }


  /* inverse collapse times for the required outputs */
  for (i=0; i<outputs.n; i++)
    /* if F is the inverse growing mode: */
    /* outputs.F[i]=1./GrowingMode(outputs.z[i]); */
    outputs.F[i]=1.+outputs.z[i];
  outputs.Flast=outputs.F[outputs.n-1];

  
  if (!ThisTask)
    {
      dprintf(VMSG, 0, "Flag for this run: %s\n\n",params.RunFlag);
      dprintf(VMSG, 0, "PARAMETER VALUES from file %s:\n",params.ParameterFile);
      dprintf(VMSG, 0, "Omega0                      %f\n",params.Omega0);
      dprintf(VMSG, 0, "OmegaLambda                 %f\n",params.OmegaLambda);    
      dprintf(VMSG, 0, "OmegaBaryon                 %f\n",params.OmegaBaryon);
      if (strcmp(params.TabulatedEoSfile,"no"))
	{
	  dprintf(VMSG, 0, "Dark Energy EoS will be read from file %s\n",params.TabulatedEoSfile);
	}
      else
	{
	  dprintf(VMSG, 0, "DE EoS parameters           %f %f\n",params.DEw0,params.DEwa);
	}

      dprintf(VMSG, 0, "Hubble100                   %f\n",params.Hubble100);
      dprintf(VMSG, 0, "Sigma8                      %f\n",params.Sigma8);
      dprintf(VMSG, 0, "PrimordialIndex             %f\n",params.PrimordialIndex);
      dprintf(VMSG, 0, "RandomSeed                  %d\n",params.RandomSeed);
      dprintf(VMSG, 0, "OutputList                  %s\n",params.OutputList);
      dprintf(VMSG, 0, "Number of outputs           %d\n",outputs.n);
      dprintf(VMSG, 0, "Output redshifts           ");
      for (i=0; i<outputs.n; i++)
	dprintf(VMSG, 0, " %f ",outputs.z[i]);
      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "GridSize                    %d %d %d\n",params.GridSize[0],params.GridSize[1],params.GridSize[2]);
      dprintf(VMSG, 0, "BoxSize (true Mpc)          %f\n",params.BoxSize_htrue);
      dprintf(VMSG, 0, "BoxSize (Mpc/h)             %f\n",params.BoxSize_h100);
      dprintf(VMSG, 0, "Particle Mass (true Msun)   %g\n",params.ParticleMass);
      dprintf(VMSG, 0, "Particle Mass (Msun/h)      %g\n",params.ParticleMass*params.Hubble100);
      dprintf(VMSG, 0, "Inter-part dist (true Mpc)  %f\n",params.InterPartDist);
      dprintf(VMSG, 0, "Inter-part dist (Mpc/h)     %f\n",params.InterPartDist*params.Hubble100);
      dprintf(VMSG, 0, "MinHaloMass (particles)     %d\n",params.MinHaloMass);
      dprintf(VMSG, 0, "BoundaryLayerFactor         %f\n",params.BoundaryLayerFactor);
      dprintf(VMSG, 0, "MaxMem per task (Mb)        %d\n",params.MaxMem);
      dprintf(VMSG, 0, "MaxMem per particle (b)     %f\n",params.MaxMemPerParticle);
      dprintf(VMSG, 0, "CatalogInAscii              %d\n",params.CatalogInAscii);
      dprintf(VMSG, 0, "NumFiles                    %d\n",params.NumFiles);
      dprintf(VMSG, 0, "DoNotWriteCatalogs          %d\n",params.DoNotWriteCatalogs);
      dprintf(VMSG, 0, "DoNotWriteHistories         %d\n",params.DoNotWriteHistories);
      dprintf(VMSG, 0, "WriteSnapshot               %d\n",params.WriteSnapshot);
      dprintf(VMSG, 0, "OutputInH100                %d\n",params.OutputInH100);
      dprintf(VMSG, 0, "WriteFmax                   %d\n",params.WriteFmax);
      dprintf(VMSG, 0, "WriteVmax                   %d\n",params.WriteVmax);
      dprintf(VMSG, 0, "WriteRmax                   %d\n",params.WriteRmax);
      switch(params.AnalyticMassFunction)
	{
	case 0:
	  dprintf(VMSG, 0, "Using Press & Schechter (1974) for the analytic mass function\n");
	  break;
	case 1:
	  dprintf(VMSG, 0, "Using Sheth & Tormen (2001) for the analytic mass function\n");
	  break;
	case 2:
	  dprintf(VMSG, 0, "Using Jenkins et al. (2001) for the analytic mass function\n");
	  break;
	case 3:
	  dprintf(VMSG, 0, "Using Warren et al. (2006) for the analytic mass function\n");
	  break;
	case 4:
	  dprintf(VMSG, 0, "Using Reed et al. (2007) for the analytic mass function\n");
	  break;
	case 5:
	  dprintf(VMSG, 0, "Using Crocce et al. (2010) for the analytic mass function\n");
	  break;
	case 6:
	  dprintf(VMSG, 0, "Using Tinker et al. (2008) for the analytic mass function\n");
	  break;
	case 7:
	  dprintf(VMSG, 0, "Using Courtin et al. (2010) for the analytic mass function\n");
	  break;
	case 8:
	  dprintf(VMSG, 0, "Using Angulo et al. (2012) for the analytic mass function\n");
	  break;
	case 9:
	  dprintf(VMSG, 0, "Using Watson et al. (2013) for the analytic mass function\n");
	  break;
	case 10:
	  dprintf(VMSG, 0, "Using Crocce et al. (2010) with forced universality for the analytic mass function\n");
	  break;
	default:
	  dprintf(VMSG, 0, "Unknown value for AnalyticMassFunction, Using Watson et al. (2013)\n");
	  params.AnalyticMassFunction=9;
	  break;
	}
      dprintf(VMSG, 0, "\n");

      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "GENIC parameters:\n");
      dprintf(VMSG, 0, "InputSpectrum_UnitLength_in_cm %f\n",params.InputSpectrum_UnitLength_in_cm);
      dprintf(VMSG, 0, "FileWithInputSpectrum          %s\n",params.FileWithInputSpectrum);
      dprintf(VMSG, 0, "WDM_PartMass_in_kev            %f\n",params.WDM_PartMass_in_kev);
      dprintf(VMSG, 0, "\n");
    }

  /* Task 0 may have changed the value of this parameter */
  MPI_Bcast(&params.AnalyticMassFunction, 1, MPI_INT, 0, MPI_COMM_WORLD);

  return 0;
}


#define NSIGMA ((double)6.0)
#define STEP_VAR ((double)0.15)

int set_smoothing()
{
  int ismooth;
  double var_min, var_max, rmin;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=-1;   /* here we want to use the standard, scale-independent growing mode */
#endif
  
  var_min    = 1.686 / NSIGMA / GrowingMode( outputs.zlast );
  var_min   *= var_min;
  /* rmax       = Radius(var_min); */
  rmin       = params.InterPartDist / 6.0;
  var_max    = MassVariance( rmin );
  Smoothing.Nsmooth = (log10(var_max) - log10(var_min)) / STEP_VAR + 2;

  if ( Smoothing.Nsmooth <= 0 )
    {
      if ( !ThisTask )
	dprintf(VERR, 0, "I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
      Smoothing.Nsmooth = 1;
    }

  if ( !ThisTask )
    dprintf(VMSG, 0, "Min variance: %f12.6, max variance: %f12.6, number of smoothing radii: %d\n",
	   var_min,var_max,Smoothing.Nsmooth);

  Smoothing.Radius      =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.Variance    =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.TrueVariance=(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  if ( Smoothing.Radius == 0x0 || Smoothing.Variance == 0x0 || Smoothing.TrueVariance == 0x0)
    {
      dprintf(VXERR, 0, "ERROR on task %d: allocation of Smoothing failed\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  for ( ismooth = 0; ismooth < Smoothing.Nsmooth - 1; ismooth++ )
    {
      Smoothing.Variance[ismooth] = pow(10., log10(var_min) + STEP_VAR * ismooth);
      Smoothing.Radius[ismooth]   = Radius(Smoothing.Variance[ismooth]);
    }
  Smoothing.Radius[ismooth]   = 0.0;
  Smoothing.Variance[ismooth] = var_max;

  if (!ThisTask)
    for ( ismooth = 0; ismooth < Smoothing.Nsmooth; ismooth++ )
      dprintf(VMSG, 0, "           %2d)  Radius=%10f, Variance=%10f\n",ismooth+1,Smoothing.Radius[ismooth],Smoothing.Variance[ismooth]);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


int generate_densities()
{

  cputime.dens = MPI_Wtime();

  if (!ThisTask)
    dprintf(VMSG, 0, "[%s] Generating density in Fourier space\n",fdate());

#ifdef WHITENOISE

  if (Ngrids>1)
    {
      if (!ThisTask)
	dprintf(VMSG, 0, "Sorry, this works only with a single grid\n");
      return 1;
    }

  if (read_white_noise())
    return 1;

#else

  int igrid;
  for (igrid=0; igrid<Ngrids; igrid++)
    if (GenIC_large(igrid))
      return 1;

#endif

  cputime.dens = MPI_Wtime()-cputime.dens;
    if (!ThisTask)
      dprintf(VMSG, 0, "[%s] Done generating density in Fourier space, cputime = %f s\n",fdate(), cputime.dens);

  return 0;
}


int set_grids()
{
  /* initialization of fftw quantities on grids (one for the moment) */

  Ngrids = 1;

  if( ( MyGrids = (grid_data*) malloc(Ngrids * sizeof(grid_data)) ) == NULL )
    {
      dprintf(VXERR, ThisTask,
	      "unable to allocate MyGrids (%zd bytes)\n", Ngrids * sizeof(grid_data));
      MPI_Finalize();
#ifdef USE_GPERFTOOLS
      ProfilerStop();
#endif
    }


  MyGrids[0].GSglobal[_x_] = params.GridSize[0];
  MyGrids[0].GSglobal[_y_] = params.GridSize[1];
  MyGrids[0].GSglobal[_z_] = params.GridSize[2];

  MyGrids[0].BoxSize        = params.BoxSize_htrue;
  MyGrids[0].lower_k_cutoff = 0.;
  MyGrids[0].upper_k_cutoff = NYQUIST * PI;
  
  /* allocates pointers */
  cvector_fft        = (pfft_complex**) malloc(Ngrids * sizeof(fftw_complex*));
  rvector_fft        = (double**)       malloc(Ngrids * sizeof(double*));

  kdensity           = (double**)       malloc(Ngrids * sizeof(double*));
  density            = (double**)       malloc(Ngrids * sizeof(double*));
  first_derivatives  = (double***)      malloc(Ngrids * sizeof(double**));
  second_derivatives = (double***)      malloc(Ngrids * sizeof(double**));
  VEL_for_displ      = (double**)       malloc(3 * sizeof(double*));
#ifdef TWO_LPT
  VEL2_for_displ     = (double**)       malloc(3 * sizeof(double*));
#endif
  
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    {
      first_derivatives[igrid]  = (double**) malloc(3 * sizeof(double*));
      second_derivatives[igrid] = (double**) malloc(6 * sizeof(double*));
    }
  // moved to GenIC
  //  seedtable=(unsigned int**)malloc(Ngrids * sizeof(unsigned int*));
 
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    if ( set_one_grid(igrid) )
      return 1;

  return 0;
}



#ifdef PLC

double myz;
int cone_and_cube_intersect(double *, double *, double *, double *, double , double *, double *, int *);
double maxF(double *, double *, double *, double *, double *);

int set_plc(void)
{
  int NAll,ir,jr,kr,ic,this,intersection,axis;
  double Largest_r,Smallest_r,smallestr,largestr,x[3],l[3],z,d,mod;
  FILE *fout;
  char filename[BLENGTH];

  /* ordering of coordinates to accomodate for rotation caused by fft ordering */
#ifdef ROTATE_BOX
  static int rot[3] = {1,2,0};
#else
  static int rot[3] = {0,1,2};
#endif

  if (params.StartingzForPLC<0.)
    {
      plc.Nreplications = 0;
      plc.Fstart=plc.Fstop = -1.;
      plc.Nmax=0;
      if (!ThisTask)
	printf("Negative value of StartingzForPLC, no Past Light Cone output will be given\n\n");

      return 0;
    }


  /* here we define the vertex and axis direction of the cone */
  if (params.PLCProvideConeData)
    {
      /* in this case data are provided in the parameter file */
      plc.center[rot[0]] = params.PLCCenter[0] / params.BoxSize * params.GridSize[0];
      plc.center[rot[1]] = params.PLCCenter[1] / params.BoxSize * params.GridSize[0];
      plc.center[rot[2]] = params.PLCCenter[2] / params.BoxSize * params.GridSize[0];
      plc.zvers[rot[0]]  = params.PLCAxis[0];
      plc.zvers[rot[1]]  = params.PLCAxis[1];
      plc.zvers[rot[2]]  = params.PLCAxis[2];
    }
  else
    {
      /* in this case the center is randomly placed and the direction points toward the main diagonal */
      gsl_rng_set(random_generator, params.RandomSeed);
      plc.center[0] = gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_x_];
      plc.center[1] = gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_y_];
      plc.center[2] = gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_z_];
      plc.zvers[0] = 1.0;
      plc.zvers[1] = 1.0;
      plc.zvers[2] = 1.0;
    }
  
  /* normalization of the cone axis direction */
  mod = sqrt( plc.zvers[0]*plc.zvers[0] + plc.zvers[1]*plc.zvers[1] + plc.zvers[2]*plc.zvers[2] );
  for (ir=0;ir<3;ir++)
    plc.zvers[ir] /= mod;

  /* here we define a system where zvers is the z-axis */
  if (plc.zvers[2] == 1.0)
    {
      /* if the zversor corresponds to the z axis use the existing axes */
      plc.xvers[0] = 1.0;
      plc.xvers[1] = 0.0;
      plc.xvers[2] = 0.0;
      plc.yvers[0] = 0.0;
      plc.yvers[1] = 1.0;
      plc.yvers[2] = 0.0;
    }
  else
    {
      /* x axis will be oriented as the cross product of zvers and the z axis */
      mod = sqrt( plc.zvers[0]*plc.zvers[0]+plc.zvers[1]*plc.zvers[1] );
      plc.xvers[0] =  plc.zvers[1] / mod;
      plc.xvers[1] = -plc.zvers[0] / mod;
      plc.xvers[2] =  0.0;
      /* y axis will be the cross product of z and x */
      plc.yvers[0] = plc.zvers[1]*plc.xvers[2] - plc.zvers[2]*plc.xvers[1];
      plc.yvers[1] = plc.zvers[2]*plc.xvers[0] - plc.zvers[0]*plc.xvers[2];
      plc.yvers[2] = plc.zvers[0]*plc.xvers[1] - plc.zvers[1]*plc.xvers[0];
    }

  /* initialization to compute the number of realizations */
  NAll = (int)(ProperDistance(params.StartingzForPLC)/MyGrids[0].BoxSize) + 2;
  plc.Fstart = 1.+params.StartingzForPLC;
  plc.Fstop  = 1.+params.LastzForPLC;
  /* if F is the inverse collapse time: */
  /* plc.Fstart = 1./GrowingMode(params.StartingzForPLC); */
  /* plc.Fstop  = 1./GrowingMode(params.LastzForPLC); */

  Largest_r=ProperDistance(params.StartingzForPLC)/params.InterPartDist;
  Smallest_r=ProperDistance(params.LastzForPLC)/params.InterPartDist;

  l[0]=(double)(MyGrids[0].GSglobal[_x_]);
  l[1]=(double)(MyGrids[0].GSglobal[_y_]);
  l[2]=(double)(MyGrids[0].GSglobal[_z_]);

  /* first, it counts the number of replications needed */
  plc.Nreplications=0;
  for (ir=-NAll; ir<=NAll; ir++)
    for (jr=-NAll; jr<=NAll; jr++)
      for (kr=-NAll; kr<=NAll; kr++)
	{
	  x[0]=ir*l[0];
	  x[1]=jr*l[1];
	  x[2]=kr*l[2];

	  intersection = cone_and_cube_intersect(x, l, plc.center, plc.zvers, params.PLCAperture, &smallestr, &largestr, &axis);

	  if (intersection && !(smallestr>Largest_r || largestr<Smallest_r))
	    plc.Nreplications++;
	}

  /* second, it allocates the needed memory for the replications */
  plc.repls=malloc(plc.Nreplications * sizeof(replication_data));

  if (!ThisTask)
    {
      sprintf(filename,"pinocchio.%s.geometry.out",params.RunFlag);
      fout=fopen(filename,"w");
      fprintf(fout, "# N. replications: %d (out of %d checked)\n",plc.Nreplications,(2*NAll+1)*(2*NAll+1)*(2*NAll+1));
      fprintf(fout, "# distance range: %10.6f %10.6f\n",Smallest_r,Largest_r);
      fprintf(fout, "# V   = %10.6f %10.6f %10.6f\n",plc.center[0],plc.center[1],plc.center[2]);
      fprintf(fout, "# D   = %10.6f %10.6f %10.6f\n",plc.zvers[0],plc.zvers[1],plc.zvers[2]);
      fprintf(fout, "# L   = %10.6f %10.6f %10.6f\n",l[0],l[1],l[2]);
      fprintf(fout, "# A   = %10.6f\n",params.PLCAperture);
      fprintf(fout, "# IPD = %10.6f\n",params.InterPartDist);
      fprintf(fout, "#\n");
    }

  /* third, it stores information on replications */
  this=0;
  for (ir=-NAll; ir<=NAll; ir++)
    for (jr=-NAll; jr<=NAll; jr++)
      for (kr=-NAll; kr<=NAll; kr++)
	{
	  x[0]=ir*l[0];
	  x[1]=jr*l[1];
	  x[2]=kr*l[2];

	  intersection = cone_and_cube_intersect(x, l, plc.center, plc.zvers, params.PLCAperture, &smallestr, &largestr, &axis);

	  if (intersection && !(smallestr>Largest_r || largestr<Smallest_r))
	    {
	      if (!ThisTask)
		fprintf(fout," %3d  %3d %3d %3d   %10.6f %10.6f   %d  %d\n",this,ir,jr,kr,smallestr,largestr,intersection,axis);
	      plc.repls[this].i=ir;
	      plc.repls[this].j=jr;
	      plc.repls[this].k=kr;
	      plc.repls[this].F1=-largestr;
	      plc.repls[this++].F2=-smallestr;
	    }	  
	}
  if (!ThisTask)
    fclose(fout);

  /* fourth it transforms distances to redshifts */
  for (z=100.; z>=0.0; z-=0.01)
    {
      d=ProperDistance(z)/params.InterPartDist;
      for (this=0; this<plc.Nreplications; this++)
	{
	  if (plc.repls[this].F1<=0.0 && d<-plc.repls[this].F1)
	    plc.repls[this].F1=z+0.01+1.0;
	  if (plc.repls[this].F2<=0.0 && d<-plc.repls[this].F2)
	    plc.repls[this].F2=z-0.01+1.0;
	}
    }
  for (this=0; this<plc.Nreplications; this++)
    {
      if (plc.repls[this].F1<=0.0)
	plc.repls[this].F1=1.0;
      if (plc.repls[this].F2<=0.0)
	plc.repls[this].F2=1.0;
    }

  plc.Nmax = subbox.Npart / 10;

  if (!ThisTask)
    {
      dprintf(VMSG, 0, "The Past Light Cone will be reconstruct from z=%f to z=%f\n", params.StartingzForPLC,params.LastzForPLC);
      if (params.PLCProvideConeData)
      	dprintf(VMSG, 0, "Cone data have been provided in the parameter file\n");
      else
       	dprintf(VMSG, 0, "Cone data have been decided by the code\n");
      dprintf(VMSG, 0, "Past Light Cone will be centred on point [%f,%f,%f] (true Mpc)\n",
	     plc.center[0]*params.InterPartDist,plc.center[1]*params.InterPartDist,plc.center[2]*params.InterPartDist);
      dprintf(VMSG, 0, "The cone vertex will be pointed toward [%f,%f,%f]\n",plc.zvers[0],plc.zvers[1],plc.zvers[2]);
      dprintf(VMSG, 0, "It will have an aperture of %f degrees\n",params.PLCAperture);
#ifdef ROTATE_BOX
      if (params.PLCProvideConeData)
	printf("(NB: rotation has been applied to the provided coordinates)\n");
#endif
      dprintf(VMSG, 0, "The proper distance at the starting redshift, z=%f, is: %f Mpc\n",
	     params.StartingzForPLC, Largest_r*params.InterPartDist);
      dprintf(VMSG, 0, "The proper distance at the stopping redshift, z=%f, is: %f Mpc\n",
	     params.LastzForPLC, Smallest_r*params.InterPartDist);
      dprintf(VMSG, 0, "The reconstruction will be done for %f < z < %f\n",params.LastzForPLC,params.StartingzForPLC);
      dprintf(VMSG, 0, "The corresponding F values are: Fstart=%f, Fstop=%f\n",plc.Fstart,plc.Fstop);
      dprintf(VMSG, 0, "The box will be replicated %d times to construct the PLC\n",plc.Nreplications);
      for (ic=0; ic<plc.Nreplications; ic++)
	printf("   Replication %2d: shift (%2d,%2d,%2d), from F=%f to F=%f\n",
	       ic,plc.repls[ic].i,plc.repls[ic].j,plc.repls[ic].k,
	       plc.repls[ic].F1,plc.repls[ic].F2);
      dprintf(VMSG, 0, "Task 0 will use plc.Nmax=%d\n",plc.Nmax);
      dprintf(VMSG, 0, "\n");
    }

  return 0;

}


double maxF(double *P, double *V, double *U, double *D, double *L)
{
  /* determines the smallest angle between a segment and a cone direction;
     the cone vertex is in a point vec(V), its direction is versor(D);
     the segment starts from point vec(P) and goes along the direction versor(U), for a length L;
     the code returns the cos of the smallest angle between the cone axis and
     the line that joins vec(V) with a point of the segment.
  */

  double dP=sqrt((P[0]-V[0])*(P[0]-V[0])+(P[1]-V[1])*(P[1]-V[1])+(P[2]-V[2])*(P[2]-V[2]));
  if (dP==0.0)
    return 1.0;
  double cosDU=(D[0]*U[0]+D[1]*U[1]+D[2]*U[2]);
  double cosDP=(D[0]*(P[0]-V[0])+D[1]*(P[1]-V[1])+D[2]*(P[2]-V[2]))/dP;
  double cosUP=(U[0]*(P[0]-V[0])+U[1]*(P[1]-V[1])+U[2]*(P[2]-V[2]))/dP;

  if (cosDP-cosDU*cosUP==0.0)
    return 0.0;
  double tmax=(cosDU-cosDP*cosUP)/(cosDP-cosDU*cosUP);
  if (tmax<0)
    tmax=0.0;
  else if (tmax>*L/dP)
    tmax=*L/dP;

  return (cosDP+tmax*cosDU)/sqrt(1.0+tmax*tmax+2.*tmax*cosUP);
}

int cone_and_cube_intersect(double *Oc, double *L, double *V, double *D, double theta, double *rmin, double *rmax, int *axis)
{
  int i, j, k, ivec[3], dim, dim1, dim2;
  double r, x, F, Fmax, costh, proj, U[3], P[3];

  /* This routine returns 1 if the cone with vertex V, axis direction D and semi-aperture theta (deg)
     intersects the cube starting from point Oc and with edges of lenght L aligned with the axes.
     It returns 0 if the two solids do not intersect.
     It also computes the smallest and largest distances of the cube from the cone vertex V.
   */


  /* Computation of rmin and rmax */
  *rmin = 1.e32;
  *rmax = 0.0;

  /* max distance from vertices */
  for ( i = 0; i < 2; i++ )
    for ( j = 0; j < 2; j++ )
      for ( k = 0; k < 2; k++ )
	{
	  r = sqrt( pow(Oc[0] + i*L[0] - V[0], 2.0) +
		    pow(Oc[1] + j*L[1] - V[1], 2.0) +
		    pow(Oc[2] + k*L[2] - V[2], 2.0) );
	  if ( r > *rmax )
	    *rmax = r;
	}


  /* min distance from cube faces */
  /* and intersection of axis with cube faces */  
  *axis = 0;
  for ( dim = 0; dim < 3; dim++ )  /* three dimension (normals to cube faces) */
    for ( i = 0; i < 2; i++ )      /* two faces per dimension */
      {
	proj = Oc[dim] - V[dim] + i*L[dim];
	dim1 = (dim+1) % 3;
	dim2 = (dim+2) % 3;
	
	/* minimum distance */
	r = proj*proj;                          /* the normal component always contributes */
	if ( V[dim1] < Oc[dim1] )               /* only of their projection is outside the face */
	  r += pow( V[dim1]-Oc[dim1], 2.0);
	
	else if ( V[dim1] > Oc[dim1] + L[dim1] )
	  r += pow( V[dim1]-Oc[dim1]-L[dim1], 2.0 );
	
	if ( V[dim2] < Oc[dim2] )
	  r += pow( V[dim2]-Oc[dim2], 2.0 );
	
	else if ( V[dim2] > Oc[dim2]+L[dim2] )
	  r += pow( V[dim2]-Oc[dim2]-L[dim2], 2.0 );
	
	r = sqrt(r);
	
	if ( r < *rmin )
	  *rmin = r;

	/* axis intersection */
	if ( (x=proj/D[dim]) > 0.0 &&
	     V[dim1] + x*D[dim1] >= Oc[dim1] &&
	     V[dim1] + x*D[dim1] <= Oc[dim1] + L[dim1] &&
	     V[dim2] + x*D[dim2] >= Oc[dim2] &&
	     V[dim2] + x*D[dim2] <= Oc[dim2] + L[dim2] )
	  *axis+=1<<(dim+i*3);
      }

  /* step 1: if the whole sky is required, only rmin and rmax are needed */
  if (theta>=180.)
    return 1;

  /* step 2: if the vertex V is inside the cube and axis>0 then they intersect */
  if (( *axis &&
	V[0]>=Oc[0] && V[0]<=Oc[0]+L[0] &&
	V[1]>=Oc[1] && V[1]<=Oc[1]+L[1] &&
	V[2]>=Oc[2] && V[2]<=Oc[2]+L[2] ) )
    {
      *rmin = 0.0;  /* in this case rmin is unrelated to the cube boundary */
      return 2;
    }

  /* step 3: if the axis intersects one face then there is an intersection */
  if ( *axis )
    return 3;

  /* step4: compute maximum of ** F = (P-V) dot D /|P-V| - cos theta ** 
     for each cube edge */
  Fmax=-10.0;
  costh = cos( theta / 180. * PI );
  for ( i = 0; i < 2; i++ )
    for (j = 0; j < 2 ; j++ )
      for ( k = 0; k < 2; k++ )
	{
	  ivec[0] = i;
	  ivec[1] = j;
	  ivec[2] = k;
	  
	  for ( dim = 0; dim < 3; dim++)
	    if ( !ivec[dim] )
	      {
		U[dim]       = 1.0;
		U[(dim+1)%3] = 0.0;
		U[(dim+2)%3] = 0.0;
		P[0]         = Oc[0] + ivec[0]*L[0];
		P[1]         = Oc[1] + ivec[1]*L[1];
		P[2]         = Oc[2] + ivec[2]*L[2];
		F            = maxF( P, V, U, D, L+dim ) - costh;
		if ( F > Fmax )
		  Fmax=F;
	      }
	}

  /* if the nearest vertex is inside the cone then exit */
  if (Fmax>0)
    return 4;

  /* at this point the cone and the cube do not intersect */
  return 0;

}


#else

int set_plc()
{
  if (!ThisTask)
    dprintf(VMSG, 0, "PLC flag at compilation was not set, no Past Light Cone output will be given\n\n");

  return 0;
}

#endif


/* division in sub-boxes */
int set_subboxes()
{

  int i1,j1,k1, i2,j2,k2,NS2, NN,NN1,N1,N2,N3,ssafe;
  double tdis,size,this,BytesPerParticle,FmaxBPP,FragG_BPP,FragP_BPP,
    TotalNP,TotalNP_pertask,ratio,smallest, MemPerTask;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=-1;
#endif

  /* typical displacement at zlast */
  tdis = GrowingMode(outputs.zlast) * sqrt( DisplVariance(params.InterPartDist) );

  /* mass of the largest halo expected in the box */
  {
    params.Largest = 1.e18;
    
    double cc = 1.0 / (params.BoxSize_htrue * params.BoxSize_htrue * params.BoxSize_htrue);    
    double aa = AnalyticMassFunction(params.Largest, outputs.zlast);
    
    while ( aa*params.Largest < cc)
      {
	params.Largest *= 0.99; 
	aa = AnalyticMassFunction(params.Largest, outputs.zlast);
      }
  }
    
  /* boundary layer */
  size                = SizeForMass(params.Largest);
  subbox.SafetyBorder = params.BoundaryLayerFactor * size;
  subbox.safe         = (int)(subbox.SafetyBorder / params.InterPartDist) + 1;

  if (!ThisTask)
    {
      dprintf(VMSG, ThisTask, "\n");
      dprintf(VMSG, ThisTask, "Determination of the boundary layer\n");
      dprintf(VMSG, ThisTask, "   growing mode at z=%f: %f\n",
	      outputs.zlast, GrowingMode( outputs.zlast ));
      dprintf(VMSG, ThisTask, "   largest halo expected in this box at z=%f: %e Msun\n",
	      outputs.zlast, params.Largest);
      dprintf(VMSG, ThisTask, "   its Lagrangian size: %f Mpc\n",
	      size);
      dprintf(VMSG, ThisTask, "   typical displacement: %f \n",
	      tdis);
      dprintf(VMSG, ThisTask, "   the boundary layer will be %f, a factor of %f with respect to the typical displacement\n",
	      subbox.SafetyBorder, subbox.SafetyBorder/tdis);
    }

  /* finds the optimal number of sub-boxes to use for the fragmentation */
  ssafe           = 2.0 * subbox.safe;
  FmaxBPP         = (double)sizeof(product_data) +
    10.0*(double)sizeof(double) + 
    (double)sizeof(int) * (double)NTasks /
    (double)MyGrids[0].GSglobal[_z_];
  
  TotalNP         = (double)( MyGrids[0].GSglobal[_x_] * MyGrids[0].GSglobal[_y_] * MyGrids[0].GSglobal[_z_] );
  TotalNP_pertask = TotalNP / (double)NTasks;

  FragP_BPP       = (double)sizeof(product_data);
  FragG_BPP       = 3.0*(double)sizeof(int)+(double)(sizeof(group_data) +sizeof(histories_data))/10.0;
#ifdef PLC
  FragG_BPP      += (double)sizeof(plcgroup_data)/10.;
#endif

  smallest = 1.e10;
  NSlices  = 0;

  do
    {
      ++NSlices;

      BytesPerParticle = 1e10;
      
      for ( int k = 1; k <= NTasks; k++ )
	for ( int j = 1; j <= NTasks/k; j++ )
	  for ( int i = 1; i <= NTasks/k/j; i++ )
	    /* the three indices must be exact divisors of the three grid lengths */
	    if ( i*j*k == NTasks )
	      {
		/* number of particles in the sub-box */
		N1 = find_length(MyGrids[0].GSglobal[_x_], i, 0);
		N2 = find_length(MyGrids[0].GSglobal[_y_], j, 0);
		N3 = find_length(MyGrids[0].GSglobal[_z_], k*NSlices, 0);
		
		if ( N1 < ssafe || N2 < ssafe || N3 < ssafe )
		  continue;
		
		NN = (N1 + (i==1? 0 : ssafe)) 
		  *  (N2 + (j==1? 0 : ssafe))
		  *  (N3 + (k*NSlices==1? 0 : ssafe));

		ratio = (double)NN/TotalNP_pertask;
		
		if ( NSlices > 1 )
		  this = (double)sizeof(product_data) + ratio * (FragP_BPP + FragG_BPP);
		else		    
		  this = ( (double)sizeof(product_data) > ratio * FragG_BPP ?
			   (double)sizeof(product_data) : ratio * FragG_BPP) +
		    ratio * FragP_BPP;
		if ( this < FmaxBPP )
		  this = FmaxBPP;
		
		if ( this < smallest )
		  {
		    smallest = this;
		    i2 = i; j2 = j; k2 = k; NS2 = NSlices;
		  }

		if ( this < BytesPerParticle )
		  {
		    BytesPerParticle = this;
		    NN1 = NN;
		    i1  = i;
		    j1  = j;
		    k1  = k;
		  }
	      }
      if ( BytesPerParticle > 1000.0 )
	break;
    }
  while ( BytesPerParticle > params.MaxMemPerParticle );


  if ( BytesPerParticle > 1000.0 )
    {
      if ( !ThisTask )
	{
	  dprintf( VXERR, ThisTask, "ERROR: no possible division of sub-boxes found up to Nslices=%d\n", 
		 NSlices );
	  dprintf( VXERR, ThisTask, "lowest possible value of memory per particle is %f ", smallest );
	  dprintf( VXERR, ThisTask, "found on a subdivision %d-%d-%d on %d slices\n",      i2, j2, k2, NS2 );
	  dprintf( VXERR, ThisTask, "please decrease BoundaryLayerFactor or increase MaxMemPerParticle\n" );
	  fflush( stdout );
	}
      return 1;
    }

  subbox.nbox_x  = i1;
  subbox.nbox_y  = j1;
  subbox.nbox_z_thisslice = k1;
  subbox.nbox_z_allslices = k1*NSlices;

  subbox.safe_x  = ( subbox.nbox_x > 1 ? subbox.safe : 0 );
  subbox.safe_y  = ( subbox.nbox_y > 1 ? subbox.safe : 0 );
  subbox.safe_z  = ( subbox.nbox_z_allslices > 1 ? subbox.safe : 0 );

  subbox.pbc_x   = ( subbox.nbox_x == 1 );
  subbox.pbc_y   = ( subbox.nbox_y == 1 );
  subbox.pbc_z   = ( subbox.nbox_z_allslices ==1 );

  /* this will be mybox for the first slice */
  NN             = subbox.nbox_y*subbox.nbox_z_thisslice;
  subbox.mybox_x = ThisTask / NN;
  NN1            = ThisTask - subbox.mybox_x*NN;
  subbox.mybox_y = NN1 / subbox.nbox_z_thisslice;
  subbox.mybox_z = NN1 - subbox.mybox_y*subbox.nbox_z_thisslice;

  subbox.Lgrid_x = find_length( MyGrids[0].GSglobal[_x_], subbox.nbox_x, subbox.mybox_x );
  subbox.Lgrid_y = find_length( MyGrids[0].GSglobal[_y_], subbox.nbox_y, subbox.mybox_y );
  subbox.Lgrid_z = find_length( MyGrids[0].GSglobal[_z_], subbox.nbox_z_allslices, subbox.mybox_z );

  subbox.Lgwbl_x = subbox.Lgrid_x + 2*subbox.safe_x; 
  subbox.Lgwbl_y = subbox.Lgrid_y + 2*subbox.safe_y;
  subbox.Lgwbl_z = subbox.Lgrid_z + 2*subbox.safe_z;

  subbox.Npart   = subbox.Lgwbl_x * subbox.Lgwbl_y * subbox.Lgwbl_z;

  subbox.start_x = find_start(MyGrids[0].GSglobal[_x_] ,subbox.nbox_x,subbox.mybox_x);
  subbox.start_y = find_start(MyGrids[0].GSglobal[_y_],subbox.nbox_y,subbox.mybox_y);
  subbox.start_z = find_start(MyGrids[0].GSglobal[_z_],subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.stabl_x = subbox.start_x - subbox.safe_x;
  subbox.stabl_y = subbox.start_y - subbox.safe_y;
  subbox.stabl_z = subbox.start_z - subbox.safe_z;

  subbox.overhead=(double)subbox.Npart/(double)(subbox.Lgrid_x * subbox.Lgrid_y * subbox.Lgrid_z);

  MemPerTask  = BytesPerParticle * TotalNP_pertask / 1024. / 1024. / 1024.;

  /* NSlices>1 is incompatible with WriteSnapshot */
  if (NSlices>1 && params.WriteSnapshot)
    {
      params.WriteSnapshot=0;
      if (!ThisTask)
	printf("Sorry, but snapshots cannot be written if fragmentation is done in slices\n");
    }

  /* initialization of quantities required by compute_mf */
   if (params.OutputInH100)
    mf.hfactor       = params.Hubble100;
  else
    mf.hfactor       = 1.0;
  mf.hfactor4        = pow(mf.hfactor,4.);
  mf.vol             = (double)MyGrids[0].GSglobal[_x_]*
    (double)MyGrids[0].GSglobal[_y_] *
    (double)MyGrids[0].GSglobal[_z_] *
    pow(params.InterPartDist,3.0);
  mf.mmin            = log10(params.MinHaloMass * params.ParticleMass) - 0.001*DELTAM;
  mf.mmax            = log10(params.Largest) +3.0*DELTAM;
  mf.NBIN            = (int)( (mf.mmax - mf.mmin) / DELTAM) +1;
  mf.ninbin          = (int*)calloc(mf.NBIN,sizeof(int));
  mf.ninbin_local    = (int*)calloc(mf.NBIN,sizeof(int));
  mf.massinbin       = (double*)calloc(mf.NBIN,sizeof(double));
  mf.massinbin_local = (double*)calloc(mf.NBIN,sizeof(double));

  /* messages */
  if (!ThisTask)
    {
      dprintf( VMSG, 0, "\n" );
      dprintf( VMSG, 0, "FRAGMENTATION:\n" );
      if ( NSlices > 1 )
	dprintf( VMSG, 0, "The box will be fragmented in %d slices\n",   NSlices );
      else
	dprintf( VMSG, 0, "The box will be fragmented in one go\n");
      dprintf( VMSG, 0, "Number of sub-boxes per dimension: %d %d %d\n", subbox.nbox_x,subbox.nbox_y,subbox.nbox_z_allslices );
      dprintf( VMSG, 0, "Boundary layer (true Mpc):         %f\n",       subbox.SafetyBorder );
      dprintf( VMSG, 0, "Boundary layer (gridpoints):       %d\n",       subbox.safe );
      dprintf( VMSG, 0, "Core 0 will work on a grid:        %d %d %d\n", subbox.Lgwbl_x,subbox.Lgwbl_y,subbox.Lgwbl_z );
      dprintf( VMSG, 0, "Number of particles for core 0:    %d\n",       subbox.Npart );
      dprintf( VMSG, 0, "The resolved box will be:          %d %d %d\n", subbox.Lgrid_x,subbox.Lgrid_y,subbox.Lgrid_z );
      dprintf( VMSG, 0, "Periodic boundary conditions:      %d %d %d\n", subbox.pbc_x,subbox.pbc_y,subbox.pbc_z );
      dprintf( VMSG, 0, "Required bytes per fft particle:   %f\n",       BytesPerParticle );
      dprintf( VMSG, 0, "The overhead for fragmentation is: %f\n",       subbox.overhead );
      dprintf( VMSG, 0, "Required memory per task:          %4.0fMb - Maxmem=%dMb\n", MemPerTask*1024.,params.MaxMem );
      dprintf( VMSG, 0, "\nThe mass function will be computed from Log M = %f to Log M = %f (%d bins)\n",
	       mf.mmin, mf.mmax, mf.NBIN );
      dprintf( VMSG, 0, "\n" );
    }

  if ( MemPerTask > params.MaxMem/1024.0 )
    {
      if ( !ThisTask )
	dprintf(VXERR, 0, "ERROR: your requirements overshoot the available memory per MPI task\n");
      return 1;
    }

  return 0;
}


int find_start(int L,int n,int ibox)
{
  int LL,MM;

  if (n==1)
    return 0;
  else
    {
      LL=L/n;
      MM=L%n;
      if (ibox==0)
        return 0;
      else if (ibox<=MM)
        return ibox*(LL+1);
      else
        return ibox*LL+MM;
    }

}


int find_length(int L, int n, int ibox)
{
  /* finds the length of a subbox, given the grid length L, 
     the number of subboxes n and the subbox id ibox */

  if (n==1) 
    return L;
  else
    {
      int LL = L/n;
      int MM = L%n;
      if (ibox < MM)
        return LL+1;
      else
	return LL;
    }
}




unsigned int gcd(unsigned int u, unsigned int v)
// this version of greatest common divisor taken
// from Daniel Lemire's blog
// lemire.me/blog/2013/12/26/fastest-way-to-compute-the-greatest-common-divisor/
{
    if (u == 0) return v;
    if (v == 0) return u;
    int shift = __builtin_ctz(u | v);
    u   >>= __builtin_ctz( u );
    do {
        v >>= __builtin_ctz( v );
        if (u > v) {
            unsigned int t = v;
            v = u;
            u = t;
        }  
        v = v - u;
    } while (v != 0);
    return u << shift;
}


int set_fft_decomposition(void)
{

  /* initialize task mesh for pfft */
  /* it's up to you to decide HOW to subdivide work in 3D, and then to store it in task_subdivision_3D*/

  int decomposition_done = 0;
  
  if(  internal.constrain_task_decomposition[0] +
       internal.constrain_task_decomposition[1] +
       internal.constrain_task_decomposition[2] > 0)

    // some constraints about how to decompose fft are set in parameter file
    {
      // --- check trivial errors
      // just ot be sure about trivial typos, check that none is < 0
      if (  internal.constrain_task_decomposition[0] < 0 ||
	    internal.constrain_task_decomposition[1] < 0 ||
	    internal.constrain_task_decomposition[2] < 0 )
	{
	  dprintf(VXERR, 0, "you can't constraint FFt decomposition with negative values\n");
	  return 1;
	}
      
      if ( internal.constrain_task_decomposition[0] == 0 )
	{	
	  dprintf(VXERR, 0, "you can't constraint FFt decomposition leaving first dimension to 0\n");
	  return 1;
	}
      // -------------------------
      
      // set first dimension
      internal.tasks_subdivision_3D[0] = internal.constrain_task_decomposition[0];
      internal.tasks_subdivision_3D[1] = internal.tasks_subdivision_3D[2] = 1;
      decomposition_done = 1;
      
      if( internal.constrain_task_decomposition[0] == NTasks )
	{
	  // all tasks are in dimension 1
	  internal.tasks_subdivision_dim = 1;
	  decomposition_done = 3;
	}
      
      else if( internal.constrain_task_decomposition[1] > 0 )
	{
	  // --- check trivial errors
	  if(internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1] > NTasks )
	    {
	      dprintf(VXERR, 0, "you specified a wrong fft decomposition: Dim0 x Dim1 = %d > %d tasks\n",
		      internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1], NTasks);
	      return 1;
	    }
	  // ------------------------
	  
	  internal.tasks_subdivision_3D[1] = internal.constrain_task_decomposition[1];
	  
	  if(internal.constrain_task_decomposition[0]*internal.constrain_task_decomposition[1] == NTasks)
	    {
	      // all tasks are in dimension 1 and 2
	      internal.tasks_subdivision_dim = 2;
	      internal.tasks_subdivision_3D[2] = 1;
	      decomposition_done = 3;
	    }
	  else
	    {
	      internal.tasks_subdivision_dim = 3;
	      decomposition_done = 3;
	      
	      internal.tasks_subdivision_3D[2] = NTasks /(internal.tasks_subdivision_3D[0] * internal.tasks_subdivision_3D[1]);
	      
	      if(internal.constrain_task_decomposition[2] > 0 &&
		 internal.tasks_subdivision_3D[2] != internal.constrain_task_decomposition[2])
		{
		  dprintf(VXERR, 0, "you specified a wrong fft decomposition: dim2 should be %d instead of %d\n", internal.tasks_subdivision_3D[2], internal.constrain_task_decomposition[2]);
		  return 1;
		}
	    }
	}
      else // constrain_task_decomposition[1] > 0
	decomposition_done = 1;
      
      // close if constrain_task_decomposition[1] > 0
      
    } // close constrain_task_decomposition initial if
  

  if(decomposition_done < 3)
    {
      // decomposition is still to be made or completed
      // we try to use as less dimensions as possible, maximizing contiguity

      
      if(decomposition_done == 1)
	// only the first dimension has been specified in the param file,
	// but with a number of tasks smaller than NTasks
	{	  
	  internal.tasks_subdivision_3D[1] = NTasks / internal.tasks_subdivision_3D[0];
	  internal.tasks_subdivision_3D[2] = 1;
	  internal.tasks_subdivision_dim = 2;
	  
	  if(NTasks % internal.tasks_subdivision_3D[0])
	    {
	      dprintf(VXERR, 0, "you specified a wrong fft decomposition\n");
	      return 1;
	    }  	  
	}
      else
	// no dimension has been constrained in the param file
	{
	  // NOTE : no non-cubic grids, no multi-grids
	  
	  int Ngrid = params.GridSize[0];

	  if( NTasks <= Ngrid )
	    // prefer 1D decomposition for the most obvious case
	    // that minimize communications in FFTs
	    {
	      internal.tasks_subdivision_dim = 1;
	      internal.tasks_subdivision_3D[0] = NTasks;
	      internal.tasks_subdivision_3D[2] = internal.tasks_subdivision_3D[1] = 1;
	      return 0;
	    }

	  int Ngrid2 = Ngrid * Ngrid / (DECOMPOSITION_LIMIT_FACTOR_2D * DECOMPOSITION_LIMIT_FACTOR_2D);

	  if( NTasks <= Ngrid2)
	    // check whether exact 2D pencil decomposition
	    {
	      unsigned GCD   = gcd( Ngrid, NTasks );
	      unsigned GCD_2 = gcd( Ngrid, (NTasks / GCD) );

	      if( GCD * GCD_2 != NTasks)
		// no exact decomposition is possible,
		// revert to 1D decomposition
		{
		  internal.tasks_subdivision_dim = 1;
		  internal.tasks_subdivision_3D[0] = NTasks;
		  internal.tasks_subdivision_3D[2] = internal.tasks_subdivision_3D[1] = 1;
		  return 0;
		}
	      else
		{
		  internal.tasks_subdivision_dim = 2;
		  internal.tasks_subdivision_3D[0] = GCD;
		  internal.tasks_subdivision_3D[1] = GCD_2;
		  return 0;
		}

	    }  // close if( NTasks < Ngrid2)

	  else
	    // try 3d decomposition
	    {
	      int Ngrid_limit = Ngrid / DECOMPOSITION_LIMIT_FACTOR_2D;
	      
	      unsigned GCD   = gcd( Ngrid_limit, NTasks );
	      unsigned GCD_2 = gcd( Ngrid_limit, (NTasks / GCD) );
	      
	      if( NTasks % (GCD * GCD_2) )
		{
		  dprintf(VXERR, 0, "3D decomposition is not possible\n");
		  return 1;
		}

	      internal.tasks_subdivision_dim = 3;
	      internal.tasks_subdivision_3D[0] = GCD;
	      internal.tasks_subdivision_3D[1] = GCD_2;

	      internal.tasks_subdivision_3D[2] = NTasks / (GCD * GCD_2);
	      return 0;
	    }
	  
	}
    }

  return 0;
  
  /* internal.tasks_subdivision_dim = 1; */
  
  /* if(NTasks >= 8) */
  /*   { */
  /*     int SQlimit = (int)floor(sqrt(NTasks) + 0.5); */
  /*     int Qlimit = (int)(floor(pow(NTasks, 1.0/3.0) + 0.5)); */
      
  /*     internal.tasks_subdivision_3D[0] = factor(NTasks, Qlimit); */
  /*     internal.tasks_subdivision_3D[1] = factor(NTasks / internal.tasks_subdivision_3D[0], SQlimit); */
  /*     internal.tasks_subdivision_3D[2] = NTasks / internal.tasks_subdivision_3D[0] / internal.tasks_subdivision_3D[1]; */
      
  /*     internal.tasks_subdivision_dim = 3 - (internal.tasks_subdivision_3D[2] == 1) - (internal.tasks_subdivision_3D[1] == 1); */
  /*   } */
  /* else if (NTasks < 4) */
  /*   { */
  /*     internal.tasks_subdivision_dim = 1; */
  /*     internal.tasks_subdivision_3D[0] = NTasks;   */
  /*   } */
  /* else */
  /*   {       */
  /*     internal.tasks_subdivision_dim = 2; */

  /*     if (NTasks % 2 == 0) */
  /* 	{ */
  /* 	  internal.tasks_subdivision_3D[0] = NTasks / 2; */
  /* 	  internal.tasks_subdivision_3D[1] = NTasks / internal.tasks_subdivision_3D[0]; */
  /* 	} */
  /*     else */
  /* 	{ */
  /* 	  // divides in stripes along x or y */
  /* 	  // randomizes the choice of the axis */
  /* 	  int    dir; */
  /* 	  double thistime = MPI_Wtime(); */
  /* 	  if(thistime < 1) */
  /* 	  	thistime = 1.0 / thistime; */
  /* 	  dir = (int)thistime % 2; */
	  
  /* 	  internal.tasks_subdivision_3D[(dir+1)%2] = 1; */
  /* 	  internal.tasks_subdivision_3D[dir] = NTasks;	   */
  /* 	} */
  /*   } */
}



/* int factor(int A, int limit) */
/* { */
/*   int i; */
/*   int factor = 1; */

/*   if( A % limit == 0) */
/*     return limit; */
    
/*   for (i = limit-1; i >= 2; ) */
/*     if( A % i == 0 && */
/* 	factor * i < limit) */
/*       { */
/* 	factor *= i; */
/* 	A /= i; */
/*       } */
/*     else */
/*       i--; */
  
/*   return factor; */
/* } */



/* int large_factorization(UL A, int limit) */
/* { */
/*   UL i; */
/*   UL factor = 1; */

/*   if( A % limit == 0) */
/*     return limit; */
    
/*   for (i = limit-1; i >= 2; ) */
/*     if( A % i == 0 && */
/* 	factor * i < limit) */
/*       { */
/* 	factor *= i; */
/* 	A /= i; */
/*       } */
/*     else */
/*       i--; */
  
/*   return factor; */
/* } */

/* int integer_factorization(UL A, UL **factors) */
/* { */
/*   int fcount; */
/*   UL size; */
/*   UL i, prev; */
/*   UL *myfactors; */

/*   if(A < 2) */
/*     { */
/*       *factors = (UL*)malloc(sizeof(UL) * 2); */
/*       (*factors)[0] = 1; */
/*       (*factors)[1] = 1; */
/*     return 1; */
/*     } */

/*   // guess: there will be no more than sqrt(A) factors */
/*   size = (UL)floor(sqrt((long double)A) + 1); */
/*   myfactors = (UL*)calloc(size * 2, sizeof(UL)); */

/*   for(fcount = -1, prev = 0, i = 2; i <= A;) */
/*     if(A % i == 0) */
/*       { */
/* 	if(prev < i) */
/* 	  { */
/* 	    fcount++; */
/* 	    prev = i; */
/* 	  } */
/* 	myfactors[fcount] = i; */
/* 	myfactors[size + fcount]++; */
/* 	A /= i;	 */
/*       } */
/*     else */
/*       i++; */

/*   fcount += 1; */
/*   *factors = (UL*)malloc(sizeof(UL) * 2 * fcount); */
/*   for(i = 0; i < fcount; i++) */
/*     { */
/*       (*factors)[i] = myfactors[i]; */
/*       (*factors)[fcount + i] = myfactors[size + i]; */
/*     } */

/*   free(myfactors); */
/*   return fcount; */
/* } */


/* int get_task_decomposition(int ThisGrid, int **N) */
/* { */
/*   unsigned int NT = NTasks; */
/*   int number_of_factors = 1; */
/*   int decomposition_dimensions; */
/*   UL G = MyGrids[ThisGrid].GSglobal[_x_]; */
/*   UL GxG = G*G; */
/*   UL remainder; */

/*   if(NT <= G) */
/*     // if the number of tasks used is less than the number of grid points, */
/*     // use 1-D decomposition (usual FFTW slab decomposition) */
/*     { */
/*       decomposition_dimensions = 1; */
/*       remainder = G % NT; */
/*     } */
/*   else */
/*     { */
/*       if(NT <= GxG) */
/* 	// if the number of tasks requested is less than the square of grid points, */
/* 	// use a 2-D decomposition (pencil decomposition) */
/* 	{ */
/* 	  decomposition_dimensions = 2; */
/* 	  number_of_factors = 2; */
/* 	} */
/*       else */
/* 	// otherwise, use a 3-D decomposition */
/* 	{ */
/* 	  decomposition_dimensions = 3; */
/* 	  number_of_factors = 4; */
/* 	} */

/*       remainder = (GxG) % NT; */
/*     } */
  
/*   if(remainder) */
/*     { */
/*       for(; NT >= 2; NT--) */
/* 	if(GxG % NT == 0) */
/* 	  break; */

/*       unsigned long long int NTG;       */
/*       for(NTG = NTasks; NT <= GxG; NT++) */
/* 	if(GxG % NTG == 0) */
/* 	  break; */

/*       if(ThisTask == 0) */
/* 	printf("The number %d of tasks that you requested does not factorize exactly with the grid number %llu for Grid %d.\n" */
/* 	       "At the moment, in the pfft implementation, that is not allowed.\n" */
/* 	       "We suggest you to use either %u or %llu tasks insted\n", */
/* 	       NTasks, G, ThisGrid, */
/* 	       NT, NTG); */
/*       return 0; */
/*     } */

/*   (*N) = (int*)malloc(sizeof(int) * number_of_factors); */


/*   /\*\/  */
/*     - */
/*     - if the decomposition dimensions are less than 3, i.e. we want to decompose in 1 or 2 factors */
/*     - the number of tasks, that is quite immediate to be achieved. */
/*     - */
/*     /\*\/ */
  
/*   if(number_of_factors < 4) */
/*     { */
/*       if(number_of_factors == 1) */
/* 	(*N)[0] = NT; */
/*       else */
/* 	{ */
/* 	  UL SQlimit = (UL)floor(sqrt((double)G)+0.5); */
/* 	  (*N)[0] = large_factorization(G, SQlimit); */
/* 	  (*N)[1] = G / (*N)[0]; */
/* 	} */
/*       return decomposition_dimensions; */
/*     } */

/*   /\*\/ */
/*     - */
/*     - if the decomposition dimensions are 3, i.e. the factors are 4, it is needed to work it out */
/*     - */
/*     /\*\/ */


/*   UL P_NF; */
/*   UL *P_factors, *factors; */
/*   UL product_common_divisors; */
/*   int fcount, ccount, i, j, k; */
  


/*   // factorize the number of tasks */
/*   P_NF = integer_factorization((UL)NT, &P_factors); */
  

/*   // find common divisors between GxG and NT */
/*   product_common_divisors = 1; */
/*   for(ccount = 0, i = P_NF-1; i >= 0; i--) */
/*     if(GxG % P_factors[i] == 0) */
/*       // P_factors[i] is a common divisor */
/*       { */
/* 	ccount++; */
/* 	fcount = 1; */
/* 	// now find to what power it is a common divisor */
/* 	while(fcount < P_factors[P_NF + i]) */
/* 	{ */
/* 	  GxG /= P_factors[i]; */
/* 	  product_common_divisors *= P_factors[i]; */
/* 	  if(GxG % P_factors[i] == 0) */
/* 	    fcount++; */
/* 	  else */
/* 	    P_factors[P_NF + i] = 0; */
/* 	} */
/* 	P_factors[P_NF + i] = fcount; */
/*       } */
/*     else */
/*       // P_factor[i] is not a commond divisor */
/*       P_factors[P_NF + i] = 0; */

/*   // at this moment: */
/*   // - ccount is the number of unique divisors */
/*   // - fcount is the number of total divisors (i.e. for 2^3, ccount = 1, fcount = 3) */
  
/*   // eliminate non-common divisors */
/*   for(i = 0; i < P_NF; i++) */
/*     if(P_factors[P_NF + i] == 0) */
/*       for(j = i+1; j < P_NF && P_factors[P_NF+j] == 0; j++) */
/* 	{ */
/* 	  P_factors[i] = P_factors[j]; */
/* 	  P_factors[P_NF + i] = P_factors[P_NF + j]; */
/* 	  P_factors[P_NF + j] = 0; */
/* 	} */
/*   // make the powers be correct for common divisors */
/*   for(fcount = i = 0; i < ccount; i++) */
/*     { */
/*       P_factors[ccount + i] = P_factors[P_NF + i]; */
/*       fcount += P_factors[ccount + i]; */
/*     } */

/*   // if we do not have enough single divisors, we can't factorize correctly */
/*   if(fcount < number_of_factors) */
/*     { */
/*       printf("decomposition in %d factors is impossible\n", number_of_factors); */
/*       return 0; */
/*     } */

/*   // find the final optimal factorization */
/*   else */
/*     { */
      
/*       factors = (UL*)malloc(sizeof(UL) * number_of_factors); */
  
/*       if(fcount == number_of_factors) */
/* 	for(k = 0, i = ccount; i >= 0; i--) */
/* 	  for(j = 0; j < P_factors[ccount + i]; j++) */
/* 	    factors[k++] = P_factors[i]; */
/*       else */
/* 	{  */
/* 	  int max, bit; */
/* 	  UL *ffactors; */
	  
/* 	  ffactors = (UL*)malloc(sizeof(UL) * fcount); */

/* 	  // instead of recording the factorization as couples [factor, power] */
/* 	  // write all the single factors on the array ffactors */
/* 	  for(k = i = 0; i < ccount; i++) */
/* 	    for(j = 0; j < P_factors[ccount + i]; j++) */
/* 	      ffactors[k++] = P_factors[i]; */
	  
/* 	  max = fcount; */
/* 	  bit = 0; */

/* 	  // Now we "fold" the ffactors as many times as possible, in order to */
/* 	  // end up with a list of number_of_factors factors that is as equilibrate */
/* 	  // as possible. */
/* 	  // To achieve that, we multiply tha smallest factor by the largest, the */
/* 	  // second-smaller by the second-larger and so on. */
/* 	  // In practice, we take the first 0,..,number_of_factors-1 entries in */
/* 	  // the ffactors array as the destination of final factors. */
/* 	  // ( note that at the begin, ffactors[0] is the smallest factor and by */
/* 	  // construction ffactors[i] <= ffactors[i+1] ). */
/* 	  // */
/* 	  while(max / number_of_factors > 1) */
/* 	    {	       */
/* 	      if(!(bit & 1)) */
/* 		// reverse order for multiplication */
/* 		for(i = 0; i < number_of_factors; i++) */
/* 		  ffactors[i] *= ffactors[max-1-i]; */
/* 	      else */
/* 		// same order of multiplication */
/* 		for(i = 0; i < number_of_factors; i++) */
/* 		  ffactors[i] *= ffactors[max-number_of_factors+i]; */
/* 	      // we reverse the order so that to "equilibrate" at maximum the results */
/* 	      bit++; */
/* 	      max -= number_of_factors; */
/* 	    } */

/* 	  // some factors are left */
/* 	  if(max > number_of_factors) */
/* 	    { */
/* 	      if(!(bit & 1)) */
/* 		for(i = 0; i < max-number_of_factors; i++) */
/* 		  ffactors[i] *= ffactors[max-1-i]; */
/* 	      else */
/* 		for(i = 0; i < max - number_of_factors; i++) */
/* 		  ffactors[i] *= ffactors[number_of_factors+i]; */
/* 	      max = number_of_factors; */
/* 	    } */
	  
/* 	  // copy the results */
/* 	  for(i = 0; i < number_of_factors; i++) */
/* 	    (*N)[i] = ffactors[i]; */
	  
/* 	  free(ffactors); */
	  
/*  	} */

/*       free(factors); */
/*     } */

/*  return decomposition_dimensions; */
/* } */
