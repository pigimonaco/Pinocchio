/* ######HEADER###### */

#include "pinocchio.h"

//#define BL_GRANDISSIMA

#ifdef PRECISE_TIMING  // SI PUO` TOGLIERE?

#define SET_WTIME cputime.partial = MPI_Wtime();
#define ASSIGN_WTIME(INIT, ACC) do { double ttt= MPI_Wtime(); cputime.ACC = ttt - cputime.INIT; } while(0)
#define ACCUMULATE_WTIME(INIT, ACC) do { double ttt= MPI_Wtime(); cputime.ACC += ttt - cputime.INIT; } while(0)
#else
#define SET_WTIME
#define ASSIGN_WTIME(INIT, ACC)
#define ACCUMULATE_WTIME(INIT, ACC)
#endif

int initialize_fft(void);
int init_cosmology(void);
int generate_densities(void);
int set_plc(void);
#ifdef SCALE_DEPENDENT
int set_scaledep_GM(void);
#endif
unsigned int gcd(unsigned int, unsigned int);
int set_fft_decomposition(void);

int initialization()
{

  /* timing */
  cputime.init=MPI_Wtime();

  /* this is for gsl integration */
  workspace = gsl_integration_workspace_alloc(NWINT);
  /* this is the initialization of the random number generator */
  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  /* reading of parameters from file and few other parameter initializations */
  if (set_parameters())
    return 1;

  /* call to cosmo.c: initialization of cosmological functions */
  if (initialize_cosmology())
    return 1;

  /* initialize pfft and fftw functions */
  if (initialize_fft())
    return 1;

  /* set the properties of grids and initialize FFTW quantities, including vectors */
  if (set_grids())
    return 1;

  /* computes the smoothing radii */
  if (set_smoothing())
    return 1;

  /* now it re-initializes the variance with a top-hat filter */
  WindowFunctionType=2;
  if (initialize_MassVariance())
    return 1;

  /* initializes quantities needed for the on-the-fly reconstruction of PLC */
  SET_WTIME;
  if (set_plc())
    return 1;
  ASSIGN_WTIME(partial, set_plc);

  /* computes the number of sub-boxes for fragmentation */
  SET_WTIME;
  if (set_subboxes())
    return 1;
  ASSIGN_WTIME(partial, set_subboxes);

#ifdef SCALE_DEPENDENT
  /* computes the growth rates for displacements */
  if (set_scaledep_GM())
    return 1;
#endif
  
  /* checks that parameters and directives are coherent */
  if (check_parameters_and_directives())
    return 1;

  /* estimates the size of output file */
  if (estimate_file_size())
    return 1;

  /* this barrier is set to have correct stdout in case the run cannot start */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* allocations of memory for fmax and memory tests */
  SET_WTIME;
  if (allocate_main_memory())
    return 1;
  ASSIGN_WTIME(partial, memory_allocation);

  /* initialization of fft plans */
  SET_WTIME;
  if (compute_fft_plans())
    return 1;
  ASSIGN_WTIME(partial, fft_initialization);

  /* generation of initial density field */
  if (!params.ReadProductsFromDumps) /* no generation if products are read from dumps */
    if (generate_densities())
      return 1;

  cputime.init=MPI_Wtime()-cputime.init;

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


int initialize_fft(void)
{

#ifdef USE_FFT_THREADS
  //if ( internal.nthreads_fft < 0 )
  internal.nthreads_fft = internal.nthreads_omp;
  if ( internal.nthreads_fft > 1 )
    dprintf(VMSG, 0, "Using %d threads for FFTs\n", internal.nthreads_fft );
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

  if (!ThisTask)
    dprintf(VMSG, ThisTask, "cube subdivision [%d dim]: %d x %d x %d = %d processes\n",
	    internal.tasks_subdivision_dim,
	    internal.tasks_subdivision_3D[0],
	    internal.tasks_subdivision_3D[1],
	    internal.tasks_subdivision_3D[2],
	    internal.tasks_subdivision_3D[0] *
	    internal.tasks_subdivision_3D[1] *
	    internal.tasks_subdivision_3D[2]);
  
  if ( pfft_create_procmesh(internal.tasks_subdivision_dim, MPI_COMM_WORLD, internal.tasks_subdivision_3D, &FFT_Comm) )
    {
      int all = 1;
      for(int iii = 0; iii < internal.tasks_subdivision_dim; iii++)
  	all *= internal.tasks_subdivision_3D[iii];
      
      pfft_fprintf(MPI_COMM_WORLD, stderr, "Error while creating communicator and mesh with %d processes\n", all);
      return 1;
    }

  return 0;
}


int set_parameters()
{
  int i;

  /* set default internal parameters */
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

  if(read_parameter_file())
    return 1;

  /* the smallest legitimate value of MinHaloMass is 1 */
  if (params.MinHaloMass<=0)
    params.MinHaloMass=1;

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
  params.InterPartDist = params.BoxSize_htrue/params.GridSize[0];

  params.ParticleMass = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0 
    * pow(params.InterPartDist,3.);
  strcpy(params.DumpDir,"DumpProducts/");

  /* The Nyquist wavenumber is used in generic calls of scale-dependent growth rates */
  params.k_for_GM = PI/params.InterPartDist;

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
      dprintf(VMSG, 0, "MinHaloMass (Msun/h)        %g\n",params.MinHaloMass*params.ParticleMass*params.Hubble100);
      dprintf(VMSG, 0, "BoundaryLayerFactor         %f\n",params.BoundaryLayerFactor);
      dprintf(VMSG, 0, "MaxMem per task (Mb)        %d\n",params.MaxMem);
      dprintf(VMSG, 0, "MaxMem per particle (b)     %f\n",params.MaxMemPerParticle);
      dprintf(VMSG, 0, "CatalogInAscii              %d\n",params.CatalogInAscii);
      dprintf(VMSG, 0, "NumFiles                    %d\n",params.NumFiles);
      dprintf(VMSG, 0, "DoNotWriteCatalogs          %d\n",params.DoNotWriteCatalogs);
      dprintf(VMSG, 0, "DoNotWriteHistories         %d\n",params.DoNotWriteHistories);
      dprintf(VMSG, 0, "WriteTimelessSnapshot       %d\n",params.WriteTimelessSnapshot);
      dprintf(VMSG, 0, "OutputInH100                %d\n",params.OutputInH100);
      dprintf(VMSG, 0, "WriteDensity                %d\n",params.WriteDensity);
      dprintf(VMSG, 0, "WriteProducts               %d\n",params.WriteProducts);
      dprintf(VMSG, 0, "DumpProducts                %d\n",params.DumpProducts);
      dprintf(VMSG, 0, "ReadProductsFromDumps       %d\n",params.ReadProductsFromDumps);
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

int set_smoothing()
{
  int ismooth;
  double var_min, var_max, rmin;

  var_min    = pow(1.686/NSIGMA / GrowingMode(outputs.zlast,params.k_for_GM),2.0);
  rmin       = params.InterPartDist/6.;
  var_max    = MassVariance(rmin);
  Smoothing.Nsmooth = (log10(var_max)-log10(var_min))/STEP_VAR+2;

  if (Smoothing.Nsmooth<=0)
    {
      if (!ThisTask)
	dprintf(VERR, 0, "I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
      Smoothing.Nsmooth=1;
    }

  if (!ThisTask)
    {
      printf("\nSMOOTHING RADII\n");
      printf("Min variance: %f12.6, max variance: %f12.6, number of smoothing radii: %d\n",
	     var_min,var_max,Smoothing.Nsmooth);
    }
  Smoothing.Radius      =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.Variance    =(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.TrueVariance=(double*)malloc(Smoothing.Nsmooth * sizeof(double));
  if (Smoothing.Radius==0x0 || Smoothing.Variance==0x0 || Smoothing.TrueVariance==0x0)
    {
      printf("ERROR on task %d: allocation of Smoothing failed\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  for (ismooth=0; ismooth<Smoothing.Nsmooth-1; ismooth++)
    {
      Smoothing.Variance[ismooth] = pow(10., log10(var_min)+STEP_VAR*ismooth);
      Smoothing.Radius[ismooth]   = Radius(Smoothing.Variance[ismooth]);
    }
  Smoothing.Radius[ismooth]   = 0.0;
  Smoothing.Variance[ismooth] = var_max;

  if (!ThisTask)
    for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
      printf("           %2d)  Radius=%10f, Variance=%10f\n",ismooth+1,Smoothing.Radius[ismooth],Smoothing.Variance[ismooth]);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


int generate_densities()
{

  cputime.dens=MPI_Wtime();

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

  int igrid, dim;

  Ngrids=1;

  MyGrids=(grid_data*)malloc(Ngrids * sizeof(grid_data));

  for (dim=0; dim<3; dim++)
    MyGrids[0].GSglobal[dim] = params.GridSize[dim];
  
  MyGrids[0].Ntotal = (unsigned long long)MyGrids[0].GSglobal[_x_] * 
    (unsigned long long)MyGrids[0].GSglobal[_y_] * 
    (unsigned long long)MyGrids[0].GSglobal[_z_];

  MyGrids[0].BoxSize = params.BoxSize_htrue;
  MyGrids[0].lower_k_cutoff=0.;
  MyGrids[0].upper_k_cutoff=NYQUIST * PI;

  /* allocates pointers */
  cvector_fft=(pfft_complex**)malloc(Ngrids * sizeof(fftw_complex*));
  rvector_fft=(double**)malloc(Ngrids * sizeof(double*));

  kdensity=(double**)malloc(Ngrids * sizeof(double*));
  density=(double**)malloc(Ngrids * sizeof(double*));
  first_derivatives=(double***)malloc(Ngrids * sizeof(double**));
  second_derivatives=(double***)malloc(Ngrids * sizeof(double**));
  VEL_for_displ=(double**)malloc(3 * sizeof(double*));
#ifdef TWO_LPT
  VEL2_for_displ=(double**)malloc(3 * sizeof(double*));
#endif
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      first_derivatives[igrid]=(double**)malloc(3 * sizeof(double*));
      second_derivatives[igrid]=(double**)malloc(6 * sizeof(double*));
    }
  /* moved to GenIC */
  /* seedtable=(unsigned int**)malloc(Ngrids * sizeof(unsigned int*)); */
 
  for (igrid=0; igrid<Ngrids; igrid++)
    if (set_one_grid(igrid))
      return 1;

  /* Task 0 broadcasts its number of fft particles, that becomes the reference */
  unsigned int PPT;
  if (!ThisTask)
    PPT=MyGrids[0].total_local_size;

  MPI_Bcast(&PPT, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MyGrids[0].ParticlesPerTask=PPT;

  if ((int)(PPT * params.MaxMemPerParticle / MBYTE + 1.0) > params.MaxMem)
    {
      if (!ThisTask)
	{
	  printf("ERROR: MaxMem of %d Mb per task is insufficient to store %d bytes per particle\n",
		 params.MaxMem, (int)params.MaxMemPerParticle);
	  printf("       please increase MaxMem to at least %d\n",(int)(PPT * params.MaxMemPerParticle / MBYTE + 1.0));
	}
      return 1;
    }

  return 0;
}



#ifdef PLC
#define NSAFE 2.0

int cone_and_cube_intersect(double *, double *, double *, double *, double , double *, double *, int *);
double maxF(double *, double *, double *, double *, double *);

int set_plc(void)
{
  int NAll,ir,jr,kr,ic,this,intersection,axis;
  double Largest_r,Smallest_r,smallestr,largestr,x[3],l[3],z,d,mod,displ_variance,tdis;
  FILE *fout;
  char filename[LBLENGTH];

  /* ordering of coordinates to accomodate for rotation caused by fft ordering */
#ifdef ROTATE_BOX
  static int rot[3]={1,2,0};
#else
  static int rot[3]={0,1,2};
#endif

  if (params.StartingzForPLC<0.)
    {
      plc.Nreplications=0;
      plc.Fstart=plc.Fstop=-1.;
      plc.Nmax=0;
      plc.Nexpected=0;
      if (!ThisTask)
	printf("Negative value of StartingzForPLC, no Past Light Cone output will be given\n\n");

      return 0;
    }


  /* here we define the vertex and axis direction of the cone */
  if (params.PLCProvideConeData)
    {
      /* in this case data are provided in the parameter file */
      plc.center[rot[0]]=params.PLCCenter[0]/params.BoxSize*params.GridSize[0];
      plc.center[rot[1]]=params.PLCCenter[1]/params.BoxSize*params.GridSize[0];
      plc.center[rot[2]]=params.PLCCenter[2]/params.BoxSize*params.GridSize[0];
      plc.zvers[rot[0]]=params.PLCAxis[0];
      plc.zvers[rot[1]]=params.PLCAxis[1];
      plc.zvers[rot[2]]=params.PLCAxis[2];
    }
  else
    {
      /* in this case the center is randomly placed and the direction points toward the main diagonal */
      gsl_rng_set(random_generator, params.RandomSeed);
      plc.center[0]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_x_];
      plc.center[1]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_y_];
      plc.center[2]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal[_z_];
      double mytheta=acos(2*gsl_rng_uniform(random_generator)-1);
      double myphi=gsl_rng_uniform(random_generator)*2.0*PI;
      plc.zvers[0]=sin(mytheta)*cos(myphi);
      plc.zvers[1]=sin(mytheta)*sin(myphi);
      plc.zvers[2]=cos(mytheta);
    }
  
  /* normalization of the cone axis direction */
  mod=sqrt(plc.zvers[0]*plc.zvers[0]+plc.zvers[1]*plc.zvers[1]+plc.zvers[2]*plc.zvers[2]);
  for (ir=0;ir<3;ir++)
    plc.zvers[ir]/=mod;

  /* here we define a system where zvers is the z-axis */
  if (plc.zvers[2]==1.0)
    {
      /* if the zversor corresponds to the z axis use the existing axes */
      plc.xvers[0]=1.0;
      plc.xvers[1]=0.0;
      plc.xvers[2]=0.0;
      plc.yvers[0]=0.0;
      plc.yvers[1]=1.0;
      plc.yvers[2]=0.0;
    }
  else
    {
      /* x axis will be oriented as the cross product of zvers and the z axis */
      mod=sqrt(plc.zvers[0]*plc.zvers[0]+plc.zvers[1]*plc.zvers[1]);
      plc.xvers[0]= plc.zvers[1]/mod;
      plc.xvers[1]=-plc.zvers[0]/mod;
      plc.xvers[2]= 0.0;
      /* y axis will be the cross product of z and x */
      plc.yvers[0]=plc.zvers[1]*plc.xvers[2]-plc.zvers[2]*plc.xvers[1];
      plc.yvers[1]=plc.zvers[2]*plc.xvers[0]-plc.zvers[0]*plc.xvers[2];
      plc.yvers[2]=plc.zvers[0]*plc.xvers[1]-plc.zvers[1]*plc.xvers[0];
    }

  /* initialization to compute the number of realizations */
  NAll=(int)(ComovingDistance(params.StartingzForPLC)/MyGrids[0].BoxSize)+2;
  plc.Fstart = 1.+params.StartingzForPLC;
  plc.Fstop  = 1.+params.LastzForPLC;

  Largest_r=ComovingDistance(params.StartingzForPLC)/params.InterPartDist;
  Smallest_r=ComovingDistance(params.LastzForPLC)/params.InterPartDist;

  /* this is needed to compute the typical displacement */
  displ_variance = sqrt( DisplVariance(params.InterPartDist) ) / params.InterPartDist;
  Smallest_r -= NSAFE * GrowingMode(params.LastzForPLC,params.k_for_GM) * displ_variance;
  Smallest_r = (Smallest_r>0 ? Smallest_r : 0.);
  Largest_r += NSAFE * GrowingMode(params.StartingzForPLC,params.k_for_GM) * displ_variance;

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
      /* NSAFE times the typical displacement at z */
      tdis = NSAFE * GrowingMode(z,params.k_for_GM) * displ_variance;

      /* comoving distance at redshift z */
      d = ComovingDistance(z)/params.InterPartDist;
      for (this=0; this<plc.Nreplications; this++)
	{
	  if (plc.repls[this].F1<=0.0 && d<-plc.repls[this].F1+tdis)
	    plc.repls[this].F1=z+0.01+1.0;
	  if (plc.repls[this].F2<=0.0 && d<-plc.repls[this].F2-tdis)
	    plc.repls[this].F2=z+1.0;
	}
    }
  for (this=0; this<plc.Nreplications; this++)
    {
      if (plc.repls[this].F1<=0.0)
	plc.repls[this].F1=1.0;
      if (plc.repls[this].F2<=0.0)
	plc.repls[this].F2=1.0;
    }

  plc.Nmax = (int)(MyGrids[0].ParticlesPerTask/6 * params.PredPeakFactor);

  /* n(z) for the light cone */
  plc.delta_z=0.05;  /* NB: this is hard-coded... */
  plc.nzbins=(int)((params.StartingzForPLC-params.LastzForPLC)/plc.delta_z+0.1);
  plc.nz=(double*)calloc(plc.nzbins , sizeof(double));

  if (!ThisTask)
    {
      printf("\nThe Past Light Cone will be reconstructed from z=%f to z=%f\n",
	     params.StartingzForPLC,params.LastzForPLC);
      if (params.PLCProvideConeData)
	printf("Cone data have been provided in the parameter file\n");
      else
	printf("Cone data have been decided by the code\n");
      printf("Past Light Cone will be centred on point [%f,%f,%f] (true Mpc)\n",
	     plc.center[0]*params.InterPartDist,plc.center[1]*params.InterPartDist,plc.center[2]*params.InterPartDist);
      printf("The cone vertex will be pointed toward [%f,%f,%f]\n",plc.zvers[0],plc.zvers[1],plc.zvers[2]);
      printf("It will have an aperture of %f degrees\n",params.PLCAperture);
#ifdef ROTATE_BOX
      if (params.PLCProvideConeData)
	printf("(NB: rotation has been applied to the provided coordinates)\n");
#endif
      printf("The comoving distance at the starting redshift, z=%f, is: %f Mpc\n",
	     params.StartingzForPLC, Largest_r*params.InterPartDist);
      printf("The comoving distance at the stopping redshift, z=%f, is: %f Mpc\n",
	     params.LastzForPLC, Smallest_r*params.InterPartDist);
      printf("The reconstruction will be done for %f < z < %f\n",params.LastzForPLC,params.StartingzForPLC);
      printf("The corresponding F values are: Fstart=%f, Fstop=%f\n",plc.Fstart,plc.Fstop);
      printf("The box will be replicated %d times to construct the PLC\n",plc.Nreplications);
      for (ic=0; ic<plc.Nreplications; ic++)
	printf("   Replication %2d: shift (%2d,%2d,%2d), from F=%f to F=%f\n",
	       ic,plc.repls[ic].i,plc.repls[ic].j,plc.repls[ic].k,
	       plc.repls[ic].F1,plc.repls[ic].F2);
      printf("Task 0 will use plc.Nmax=%d\n",plc.Nmax);
      printf("The halo number density will be output in %d redshift bins\n",plc.nzbins);
      printf("\n");
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

  /* This routine returns >=1 if the cone with vertex V, axis
     direction D and semi-aperture theta (deg) intersects the cube
     (parallelepiped) starting from point Oc and with edges of lenght
     L aligned with the axes.  It returns 0 if the two solids do not
     intersect.  It also computes the smallest and largest distances
     of the cube from the cone vertex V. The bits of the integer axis
     will encode the cube faces that are intersected by the cone axis.
   */


  /* Initialization of rmin and rmax */
  *rmin=1.e32;
  *rmax=0.0;

  /* max distance is computed using the cube vertices */
  for (i=0;i<2;i++)
    for (j=0;j<2;j++)
      for (k=0;k<2;k++)
	{
	  r = sqrt(pow(Oc[0]+i*L[0]-V[0],2.0) +
		   pow(Oc[1]+j*L[1]-V[1],2.0) +
		   pow(Oc[2]+k*L[2]-V[2],2.0));
	  if (r>*rmax)
	    *rmax=r;
	}


  /* min distance from cube faces */
  /* and intersection of axis with cube faces */  
  *axis=0;
  for (dim=0; dim<3; dim++)  /* three dimension (normals to cube faces) */
    for (i=0; i<2; i++)      /* two faces per dimension */
      {
	proj=Oc[dim]-V[dim]+i*L[dim];
	dim1=(dim+1)%3;
	dim2=(dim+2)%3;
	
	/* minimum distance */
	r=proj*proj;                        /* the normal component always contributes */
	if (V[dim1]<Oc[dim1])               /* only of their projection is outside the face */
	  r+=pow(V[dim1]-Oc[dim1],2.0);
	else if (V[dim1]>=Oc[dim1]+L[dim1])
	  r+=pow(V[dim1]-Oc[dim1]-L[dim1],2.0);
	if (V[dim2]<Oc[dim2])
	  r+=pow(V[dim2]-Oc[dim2],2.0);
	else if (V[dim2]>=Oc[dim2]+L[dim2])
	  r+=pow(V[dim2]-Oc[dim2]-L[dim2],2.0);
	r=sqrt(r);
	if (r<*rmin)
	  *rmin=r;

	/* axis intersection */
	if ( (x=proj/D[dim]) > 0.0 &&
	     V[dim1] + x*D[dim1] >= Oc[dim1] &&
	     V[dim1] + x*D[dim1] < Oc[dim1] + L[dim1] &&
	     V[dim2] + x*D[dim2] >= Oc[dim2] &&
	     V[dim2] + x*D[dim2] < Oc[dim2] + L[dim2] )
	  *axis+=1<<(dim+i*3);
      }

  /* step 1: if the vertex V is inside the cube then they intersect */
  if (( V[0]>=Oc[0] && V[0]<Oc[0]+L[0] &&
	V[1]>=Oc[1] && V[1]<Oc[1]+L[1] &&
	V[2]>=Oc[2] && V[2]<Oc[2]+L[2] ) )
    {
      *rmin = 0.0;  /* in this case rmin is unrelated to the cube boundary */
      return 1;
    }

  /* step 2: if the whole sky is required, only rmin and rmax are needed */
  if (theta>=180.)
      return 2;

  /* step 3: if the axis intersects one face then there is an intersection */
  if (*axis)
    return 3;

  /* step4: compute maximum of ** F = (P-V) dot D /|P-V| - cos theta ** 
     for each cube edge */
  Fmax=-10.0;
  costh = cos( theta / 180. * PI );
  for (i=0;i<2;i++)
    for (j=0;j<2;j++)
      for (k=0;k<2;k++)
	{
	  ivec[0]=i;
	  ivec[1]=j;
	  ivec[2]=k;
	  
	  for (dim=0;dim<3;dim++)
	    if (!ivec[dim])
	      {
		U[dim]=1.0;
		U[(dim+1)%3]=0.0;
		U[(dim+2)%3]=0.0;
		P[0]=Oc[0]+ivec[0]*L[0];
		P[1]=Oc[1]+ivec[1]*L[1];
		P[2]=Oc[2]+ivec[2]*L[2];
		F=maxF(P, V, U, D, L+dim)-costh;
		if (F>Fmax)
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
    printf("PLC flag at compilation was not set, no Past Light Cone output will be given\n\n");

  return 0;
}

#endif

/* division in sub-boxes */
int set_subboxes()
{

  int i,j,k, i1=0,j1=0,k1=0, NN,NN1,N1,N2,N3;
  unsigned long long int surface,this,tt;
  double size,sizeG,cc;

  /* mass of the largest halo expected in the box */
  params.Largest=1.e18;
  cc=1./pow(params.BoxSize_htrue,3.0);
  double aa=AnalyticMassFunction(params.Largest,outputs.zlast);
  while (aa*params.Largest<cc)
    {
      params.Largest*=0.99;
      aa=AnalyticMassFunction(params.Largest,outputs.zlast);
    }
  size=SizeForMass(params.Largest);
  sizeG=size/params.InterPartDist;
  
  /*  
      The number of loadable subbox particles is equal 
      to the number allowed by the specified MaxMemPerParticle.
      The boundary layer is set to its maximum value.
  */

  /* finds the optimal number of sub-boxes to use for the fragmentation */


  surface=MyGrids[0].Ntotal;
  for (k=1; k<=NTasks; k++)
    for (j=1; j<=NTasks/k; j++)
      for (i=1; i<=NTasks/k/j; i++)
	/* the three indices must be exact divisors of the three grid lengths */
	if (i*j*k==NTasks)
	  {
	    /* number of particles in the sub-box */
	    N1 = find_length(MyGrids[0].GSglobal[_x_],i,0);
	    N2 = find_length(MyGrids[0].GSglobal[_y_],j,0);
	    N3 = find_length(MyGrids[0].GSglobal[_z_],k,0);

	    this = (unsigned long long int)(i>1? 2*(N2*N3) : 0) + 
	      (unsigned long long int)(j>1? 2*(N1*N3) : 0) +
	      (unsigned long long int)(k>1? 2*(N1*N2) : 0);
	    tt=this;
	    if (N1/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N1,2.0));
	    if (N2/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N2,2.0));
	    if (N3/2 < sizeG)
	      this+=(unsigned long long int)((double)tt*pow(2*sizeG/(double)N3,2.0));

	    if (this<surface)
	      {
		surface=this;
		i1=i; 
		j1=j; 
		k1=k; 
	      }
	  }

  subbox.nbox[_x_]=i1;
  subbox.nbox[_y_]=j1;
  subbox.nbox[_z_]=k1;

  /* mybox is the box assigned to the task */
  NN=subbox.nbox[_y_]*subbox.nbox[_z_];
  if (NN==0)
    {
      printf("ERROR: I could not find a valid subbox subdivision\n");
      printf("       subbox.nbox = [%d,%d,%d]\n",i1,j1,k1);
      printf("       please try again with a different number of tasks\n");
      return 1;
    }

  subbox.mybox[_x_]=ThisTask/NN;
  NN1=ThisTask-subbox.mybox[_x_]*NN;
  subbox.mybox[_y_]=NN1/subbox.nbox[_z_];
  subbox.mybox[_z_]=NN1-subbox.mybox[_y_]*subbox.nbox[_z_];

  subbox.Lgrid[_x_] = find_length(MyGrids[0].GSglobal[_x_],subbox.nbox[_x_],subbox.mybox[_x_]);
  subbox.Lgrid[_y_] = find_length(MyGrids[0].GSglobal[_y_],subbox.nbox[_y_],subbox.mybox[_y_]);
  subbox.Lgrid[_z_] = find_length(MyGrids[0].GSglobal[_z_],subbox.nbox[_z_],subbox.mybox[_z_]);

  subbox.pbc[_x_] = (subbox.nbox[_x_]==1);
  subbox.pbc[_y_] = (subbox.nbox[_y_]==1);
  subbox.pbc[_z_] = (subbox.nbox[_z_]==1);

/* #ifndef BL_GRANDISSIMA */

/*   subbox.safe[_x_] = (subbox.pbc[_x_] ? 0 : (find_length(MyGrids[0].GSglobal[_x_],subbox.nbox[_x_],0)-1)/2); */
/*   subbox.safe[_y_] = (subbox.pbc[_y_] ? 0 : (find_length(MyGrids[0].GSglobal[_y_],subbox.nbox[_y_],0)-1)/2); */
/*   subbox.safe[_z_] = (subbox.pbc[_z_] ? 0 : (find_length(MyGrids[0].GSglobal[_z_],subbox.nbox[_z_],0)-1)/2); */

/* #else */

  /* the boundary layer can be as large as to nearly fill the whole box,
     but the number of particles must be represented by an unsigned int */
  int BB = (int)(params.BoundaryLayerFactor*sizeG+1);
  subbox.safe[_x_] = (subbox.pbc[_x_] ? 0 : (BB > MyGrids[0].GSglobal[_x_]/2 - subbox.Lgrid[_x_]/2 - 1 ? MyGrids[0].GSglobal[_x_]/2 - subbox.Lgrid[_x_]/2 - 1 : BB));
  subbox.safe[_y_] = (subbox.pbc[_y_] ? 0 : (BB > MyGrids[0].GSglobal[_y_]/2 - subbox.Lgrid[_y_]/2 - 1 ? MyGrids[0].GSglobal[_y_]/2 - subbox.Lgrid[_y_]/2 - 1 : BB));
  subbox.safe[_z_] = (subbox.pbc[_z_] ? 0 : (BB > MyGrids[0].GSglobal[_z_]/2 - subbox.Lgrid[_z_]/2 - 1 ? MyGrids[0].GSglobal[_z_]/2 - subbox.Lgrid[_z_]/2 - 1 : BB));

//#endif

  subbox.Lgwbl[_x_] = subbox.Lgrid[_x_] + 2*subbox.safe[_x_];
  subbox.Lgwbl[_y_] = subbox.Lgrid[_y_] + 2*subbox.safe[_y_];
  subbox.Lgwbl[_z_] = subbox.Lgrid[_z_] + 2*subbox.safe[_z_];
  unsigned long long MySize = (long long)subbox.Lgwbl[_x_] * (long long)subbox.Lgwbl[_y_] * (long long)subbox.Lgwbl[_z_];
  while (MySize > (unsigned long long)1<<31)
    {
      subbox.safe[_x_] -=1;
      subbox.safe[_y_] -=1;
      subbox.safe[_z_] -=1;
      subbox.Lgwbl[_x_] = subbox.Lgrid[_x_] + 2*subbox.safe[_x_];
      subbox.Lgwbl[_y_] = subbox.Lgrid[_y_] + 2*subbox.safe[_y_];
      subbox.Lgwbl[_z_] = subbox.Lgrid[_z_] + 2*subbox.safe[_z_];
      MySize = (unsigned long long)subbox.Lgwbl[_x_] * (unsigned long long)subbox.Lgwbl[_y_] * (unsigned long long)subbox.Lgwbl[_z_];
    }
  
  subbox.start[_x_] = find_start(MyGrids[0].GSglobal[_x_],subbox.nbox[_x_],subbox.mybox[_x_]);
  subbox.start[_y_] = find_start(MyGrids[0].GSglobal[_y_],subbox.nbox[_y_],subbox.mybox[_y_]);
  subbox.start[_z_] = find_start(MyGrids[0].GSglobal[_z_],subbox.nbox[_z_],subbox.mybox[_z_]);

  subbox.stabl[_x_] = subbox.start[_x_] - subbox.safe[_x_];
  subbox.stabl[_y_] = subbox.start[_y_] - subbox.safe[_y_];
  subbox.stabl[_z_] = subbox.start[_z_] - subbox.safe[_z_];

  /* 
     Npart: total number of particles in the whole sub-volume 
     Ngood: total number of particles in the well reconstructed region
     Npredpeaks: a guess of the maximum number of peaks in the subbox
     Nalloc: number of particles for which memory has been allocated (set in organize_main_memory)
     Nstored: number of actually stored particles
  */

  subbox.Npart = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_] * subbox.Lgwbl[_z_];
  subbox.Ngood = subbox.Lgrid[_x_] * subbox.Lgrid[_y_] * subbox.Lgrid[_z_];
  /* this is a prediction of the number of peaks that will be found */
  subbox.PredNpeaks = (int)(MyGrids[0].ParticlesPerTask/6 * params.PredPeakFactor);
  subbox.Nstored = 0;
  /* this is the size of frag_map*/
  subbox.maplength = subbox.Npart/UINTLEN + (subbox.Npart%UINTLEN!=0);
  if ( (subbox.Nalloc = organize_main_memory()) == 0 )
    {
      fflush(stdout);
      if (!ThisTask)
	printf("organize_main_memory returned an invalid Nalloc, exiting\n");
      return 1;
    }

  /* messagges */
  if (!ThisTask)
    {
      printf("\n");
      printf("FRAGMENTATION:\n");
      printf("Reference number of particles:         %d\n",MyGrids[0].ParticlesPerTask);
      printf("Requested bytes per particle:          %d\n",(int)params.MaxMemPerParticle);
      printf("Number of sub-boxes per dimension:     %d %d %d\n",subbox.nbox[_x_],subbox.nbox[_y_],subbox.nbox[_z_]);
      printf("Periodic boundary conditions:          %d %d %d\n",subbox.pbc[_x_],subbox.pbc[_y_],subbox.pbc[_z_]);
      printf("Core 0 will work on a grid:            %d %d %d\n",subbox.Lgwbl[_x_],subbox.Lgwbl[_y_],subbox.Lgwbl[_z_]);
      printf("The resolved box will be:              %d %d %d\n",subbox.Lgrid[_x_],subbox.Lgrid[_y_],subbox.Lgrid[_z_]);
      printf("Boundary layer:                        %d %d %d\n",subbox.safe[_x_],subbox.safe[_y_],subbox.safe[_z_]);
      printf("Boundary layer factor:                 %f\n",params.BoundaryLayerFactor);
      printf("Number of total particles for core 0:  %d\n",subbox.Npart);
      printf("Number of good particles for core 0:   %d\n",subbox.Ngood);
      printf("Particles that core 0 will allocate:   %d\n",subbox.Nalloc);
      printf("Allowed overhead for boundary layer:   %f\n",(float)subbox.Nalloc/(float)MyGrids[0].ParticlesPerTask);
      printf("Largest halo expected in this box at z=%f: %e Msun\n",
  	     outputs.zlast, params.Largest);
      printf("   its Lagrangian size: %f Mpc (%6.2f grid points)\n",size,sizeG);
      printf("   this requires a boundary layer of %6.2f grid points \n",sizeG*params.BoundaryLayerFactor);

#ifndef BL_GRANDISSIMA 
      if ((!subbox.pbc[_x_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_x_]) || 
	  (!subbox.pbc[_y_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_y_]) || 
	  (!subbox.pbc[_z_] && params.BoundaryLayerFactor*sizeG>subbox.safe[_z_]))
	{
	  printf("WARNING: the boundary layer on some dimension is smaller than the predicted size of the largest halos\n");
	  printf("         times the BoundaryLayerFactor, the most massive halos may be inaccurate\n");
	}
#endif
    }

  /* initialization of quantities required by compute_mf */
   if (params.OutputInH100)
    mf.hfactor=params.Hubble100;
  else
    mf.hfactor=1.0;
  mf.hfactor4=pow(mf.hfactor,4.);
  mf.vol=(double)MyGrids[0].Ntotal*pow(params.InterPartDist,3.0);
  mf.mmin=log10(params.MinHaloMass*params.ParticleMass)-0.001*DELTAM;
  mf.mmax=log10(params.Largest)+3.0*DELTAM;
  mf.NBIN = (int)((mf.mmax-mf.mmin)/DELTAM) +1;
  mf.ninbin=(int*)calloc(mf.NBIN,sizeof(int));
  mf.ninbin_local=(int*)calloc(mf.NBIN,sizeof(int));
  mf.massinbin=(double*)calloc(mf.NBIN,sizeof(double));
  mf.massinbin_local=(double*)calloc(mf.NBIN,sizeof(double));

  /* messages */
  if (!ThisTask)
    {
      printf("\nThe mass function will be computed from Log M=%f to Log M=%f (%d bins)\n",
      	     mf.mmin, mf.mmax, mf.NBIN);
      printf("\n");
      fflush(stdout);      
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

  int LL,MM;

  if (n==1) 
    return L;
  else
    {
      LL=L/n;
      MM=L%n;
      if (ibox<MM)
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
}


int check_parameters_and_directives(void)
{
  

#ifndef SNAPSHOT
  if (params.WriteTimelessSnapshot || params.WriteDensity)
    {
      if (!ThisTask)
	printf("ERROR: to produce a snapshot you have to compile with SNAPSHOT directive\n");
      return 1;
    }
#endif

#ifdef SNAPSHOT

  static unsigned long long largest32 = (unsigned)1<<31;

#ifndef LONGIDS
  if ((params.WriteTimelessSnapshot || params.WriteDensity) && MyGrids[0].Ntotal > largest32)
    {
      if (!ThisTask)
	printf("ERROR: with these many particles you need to compile with LONGIDS directive\n  otherwise the snapshot IDs will be unreadable\n");
      return 1;
    }
#endif

  if (params.WriteTimelessSnapshot)
    {
      unsigned long long BlockLength = MyGrids[0].Ntotal * 12 / (unsigned long long)params.NumFiles;
      if (BlockLength > largest32)
	{
	  unsigned int NumFiles = (int)(MyGrids[0].Ntotal * 12 / largest32);
	  if ((unsigned long long)(NumFiles * 12) * largest32 < MyGrids[0].Ntotal)
	    ++NumFiles;

	  if (!ThisTask)
	    {
	      printf("ERROR: you need to write such a large snapshot with at least NumFiles=%d\n",NumFiles);
	    }
	  return 1;
	}
    }
#endif


  return 0;
}


#ifdef SCALE_DEPENDENT
#include "def_splines.h"
#define SMALLDIFF ((double)1.e-5)
#define TOLERANCE ((double)1.e-4)
#define MAXITER 20
//#define DEBUG

double ThisRadius;
double Time;

double IntegrandForSDDensVariance(double logk, void *radius)
{
  double k,w,D;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * w*w * k*k*k / (2.*PI*PI);
}

double IntegrandForSDDisplVariance(double logk, void *radius)
{
  double k,w,D;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * w*w * k / (2.*PI*PI);
}

double IntegrandForSDDispl2Variance(double logk, void *radius)
{
  double k,w,D;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_2LPT(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * w*w * k / (2.*PI*PI);
}

double IntegrandForSDDispl31Variance(double logk, void *radius)
{
  double k,w,D;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_3LPT_1(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * w*w * k / (2.*PI*PI);
}

double IntegrandForSDDispl32Variance(double logk, void *radius)
{
  double k,w,D;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_3LPT_2(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * w*w * k / (2.*PI*PI);
}


double IntegrandForSDVelVariance(double logk, void *radius)
{
  double k,w,D,fo;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode(1./Time-1.,k);
  fo=fomega(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * fo*fo * w*w * k / (2.*PI*PI);
}

double IntegrandForSDVel2Variance(double logk, void *radius)
{
  double k,w,D,fo;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_2LPT(1./Time-1.,k);
  fo=fomega_2LPT(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * fo*fo * w*w * k / (2.*PI*PI);
}

double IntegrandForSDVel31Variance(double logk, void *radius)
{
  double k,w,D,fo;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_3LPT_1(1./Time-1.,k);
  fo=fomega_3LPT_1(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * fo*fo * w*w * k / (2.*PI*PI);
}

double IntegrandForSDVel32Variance(double logk, void *radius)
{
  double k,w,D,fo;

  k=pow(10.,logk);
  w=WindowFunction(k * ThisRadius);
  D=GrowingMode_3LPT_2(1./Time-1.,k);
  fo=fomega_3LPT_2(1./Time-1.,k);
  return PowerSpectrum(k) * D*D * fo*fo * w*w * k / (2.*PI*PI);
}


int set_scaledep_GM()
{
  /* computes for each smoothing radius the k wavenumber at which the Fourier growing mode
     has the same time evolution of the rms of the density and displacement */

  double *vector;
  gsl_function Function;
  int ismooth, i, Today, Z20, iter;
  double nyquist, normv, normGM1, normGM2, normGMm, logk1, logk2, logkmid, k1, k2, kmid, 
    diff1, diff2, diffm, mindiff, result, error;

#ifdef DEBUG
  gsl_function Function2, Function31, Function32;
  FILE *fd;
#endif


#ifdef ELL_CLASSIC
  /* splines for inverse growing mode */
  SPLINE_INVGROW = (gsl_spline **)calloc(Smoothing.Nsmooth, sizeof(gsl_spline *));
  for (i=0; i<Smoothing.Nsmooth; i++)
    SPLINE_INVGROW[i] = gsl_spline_alloc (gsl_interp_cspline, NBINS);
  ACCEL_INVGROW = (gsl_interp_accel **)calloc(Smoothing.Nsmooth, sizeof(gsl_interp_accel *));
  for (i=0; i<Smoothing.Nsmooth; i++)
    ACCEL_INVGROW[i] = gsl_interp_accel_alloc();
#endif

  Smoothing.Rad_GM          = (double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.k_GM_dens       = (double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.k_GM_displ      = (double*)malloc(Smoothing.Nsmooth * sizeof(double));
  Smoothing.k_GM_vel        = (double*)malloc(Smoothing.Nsmooth * sizeof(double));

  for (i=0; i<NBINS; i++)
    if (pow(10.,SPLINE[SP_TIME]->x[i])>1.0)
      break;
  Today=i-1;

  for (i=0; i<NBINS; i++)
    if (pow(10.,SPLINE[SP_TIME]->x[i])>1./21.)
      break;
  Z20=i-1;


  /***********/
  /* density */
  /***********/
  Function.function = &IntegrandForSDDensVariance;
  vector = (double*)malloc(NBINS * sizeof(double));

  /* density (mass) variance is computed using Gaussian smoothing, 
     the integral will be truncated at the Nyquist frequency */
  WindowFunctionType=0;
  nyquist = NYQUIST * PI / params.InterPartDist;

#ifdef DEBUG
  if (!ThisTask)
    fd=fopen("accuracy_scaledependent_dens.txt","w");
#endif

  /* cycle on smoothing radii */
  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {
      /* here we use the (Gaussian) smoothing radii used in the rest of the code */
      ThisRadius=Smoothing.Radius[ismooth];

      /* integrates the growing mode and stores in vector */
      for (i=0; i<NBINS; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
  	  gsl_integration_qags (&Function, -4., nyquist, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  vector[i]=sqrt(result);
	}
      normv=vector[Today];
      for (i=0; i<NBINS; i++)
	vector[i]/=normv;

      /* bisector search of the k that gives the best approximation to the growth rate */
      iter=0;

      /* first guess at the extremes in k */
      logk1=LOGKMIN;
      logk2=LOGKMIN+(NkBINS-1)*DELTALOGK;
      k1=pow(10.,logk1);
      k2=pow(10.,logk2);

      /* differences at the extremes in k */
      normGM1=GrowingMode(0.0,k1);
      normGM2=GrowingMode(0.0,k2);
      for (diff1=diff2=0.0, i=Z20; i<=Today; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
	  diff1+=vector[i]-GrowingMode(1./Time-1.,k1)/normGM1;
	  diff2+=vector[i]-GrowingMode(1./Time-1.,k2)/normGM2;
	}
      diff1/=(double)NBINS;
      diff2/=(double)NBINS;

      /* k1 might be accurate enough */
      if (fabs(diff1) < SMALLDIFF)
	Smoothing.k_GM_dens[ismooth]=k1;

      /* k2 might be accurate enough */
      else if (fabs(diff2)<SMALLDIFF)
	Smoothing.k_GM_dens[ismooth]=k2;

      /* this is an unfortunate case: the two differences have the same sign
         and none is accurate enough. In this case the scale is set to the k that gives
	 the smallest difference, and a warning is printed out  */
      else if (diff1*diff2>0)
	{
	  if (!ThisTask)
	    printf("WARNING in scale-dependent density growth rate for smoothing radius %d (%f): accuracy not guaranteed [diff1=%g - diff2=%g]\n",
		   ismooth, Smoothing.Radius[ismooth],diff1,diff2);
	  if (fabs(diff1)<fabs(diff2))
	    Smoothing.k_GM_dens[ismooth]=logk1;
	  else
	    Smoothing.k_GM_dens[ismooth]=k2;

	}
      else
	{
	  /* smallest difference */
	  mindiff=fabs(diff1);
	  mindiff=(fabs(diff2) < mindiff ? fabs(diff2) : mindiff);

	  /* iteration for bisector search */
	  do
	    {

	      /* midpoint in log space */
	      logkmid=0.5*(logk1+logk2);
	      kmid=pow(10.,logkmid);

	      /* difference at midpoint */
	      normGMm=GrowingMode(0.0,kmid);
	      for (diffm=0.0, i=Z20; i<=Today; i++)
		{
		  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
		  diffm+=vector[i]-GrowingMode(1./Time-1.,kmid)/normGMm;
		}
	      diffm/=(double)NBINS;

	      /* update of mindiff if necessary */
	      mindiff=(fabs(diffm) < mindiff ? fabs(diffm) : mindiff);

	      /* new extreme points */
	      if (diff1*diffm>0)
		{
		  logk1=logkmid;
		  diff1=diffm;
		}
	      else
		{
		  logk2=logkmid;
		  diff2=diffm;
		}
	      ++iter;
	    }
	  while (fabs(mindiff)>SMALLDIFF && iter<=MAXITER);
    
	  Smoothing.k_GM_dens[ismooth]=kmid;
	}

#ifdef DEBUG
      printf("density, smoothing %d, iter=%d\n",ismooth,iter);
      if (!ThisTask)
	for (i=0; i<NBINS; i++)
	  {
	    Time=pow(10.,SPLINE[SP_TIME]->x[i]);
	    fprintf(fd,"  %d %g  %g  %g  %g  %g\n",
		    ismooth,Smoothing.Radius[ismooth],Smoothing.k_GM_dens[ismooth],
		    Time,vector[i],GrowingMode(1./Time-1.,Smoothing.k_GM_dens[ismooth])/GrowingMode(0.0,Smoothing.k_GM_dens[ismooth]));
	  }
#endif

#ifdef ELL_CLASSIC
      for (i=0; i<NBINS; i++)
	vector[i]=log10(vector[i]);
      gsl_spline_init(SPLINE_INVGROW[ismooth], vector, &(SPLINE[SP_TIME]->x[0]),  NBINS);
#endif

    }
#ifdef DEBUG
  if (!ThisTask)
    fclose(fd);
#endif


  /*****************/
  /* displacements */
  /*****************/
  Function.function = &IntegrandForSDDisplVariance;
#ifdef DEBUG
  Function2.function = &IntegrandForSDDispl2Variance;
  Function31.function = &IntegrandForSDDispl31Variance;
  Function32.function = &IntegrandForSDDispl32Variance;
#endif

  /* displacement variance is computed using Top-Hat smoothing, 
     the integral will be again truncated at the Nyquist frequency */
  WindowFunctionType=2;

#ifdef DEBUG
  if (!ThisTask)
    fd=fopen("accuracy_scaledependent_disp.txt","w");
#endif

  double Largest=SizeForMass(pow(10.,mf.mmax));
  /* cycle on smoothing radii */
  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {
      ThisRadius=Smoothing.Rad_GM[ismooth]=Largest*(Smoothing.Nsmooth-1-ismooth)/(double)(Smoothing.Nsmooth-1);

      /* integrates the growing mode and stores in vector */
      for (i=0; i<NBINS; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
  	  gsl_integration_qags (&Function, -4., nyquist, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  vector[i]=sqrt(result);
	}
      normv=vector[Today];
      for (i=0; i<NBINS; i++)
	vector[i]/=normv;

      /* bisector search of the k that gives the best approximation to the growth rate */
      iter=0;

      /* first guess at the extremes in k */
      logk1=LOGKMIN;
      logk2=LOGKMIN+(NkBINS-1)*DELTALOGK;
      k1=pow(10.,logk1);
      k2=pow(10.,logk2);

      /* differences at the extremes in k */
      normGM1=GrowingMode(0.0,k1);
      normGM2=GrowingMode(0.0,k2);
      for (diff1=diff2=0.0, i=Z20; i<=Today; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
	  diff1+=vector[i]-GrowingMode(1./Time-1.,k1)/normGM1;
	  diff2+=vector[i]-GrowingMode(1./Time-1.,k2)/normGM2;
	}
      diff1/=(double)NBINS;
      diff2/=(double)NBINS;

      /* k1 might be accurate enough */
      if (fabs(diff1)<SMALLDIFF)
	Smoothing.k_GM_displ[ismooth]=k1;

      /* k2 might be accurate enough */
      else if (fabs(diff2)<SMALLDIFF)
	Smoothing.k_GM_displ[ismooth]=k2;

      /* this is an unfortunate case: the two differences have the same sign
         and none is accurate enough. In this case the scale is set to the k that gives
	 the smallest difference, and a warning is printed out  */
      else if (diff1*diff2>0)
	{
	  if (!ThisTask)
	    printf("WARNING in scale-dependent density growth rate for smoothing radius %d (%f): accuracy not guaranteed [diff1=%g - diff2=%g]\n",
		   ismooth, Smoothing.Rad_GM[ismooth],diff1,diff2);
	  if (fabs(diff1)<fabs(diff2))
	    Smoothing.k_GM_displ[ismooth]=k1;

	  else
	    Smoothing.k_GM_displ[ismooth]=k2;
	}
      else
	{
	  /* smallest difference */
	  mindiff=fabs(diff1);
	  mindiff=(fabs(diff2) < mindiff ? fabs(diff2) : mindiff);

	  /* iteration for bisector search */
	  do
	    {

	      /* midpoint in log space */
	      logkmid=0.5*(logk1+logk2);
	      kmid=pow(10.,logkmid);

	      /* difference at midpoint */
	      normGMm=GrowingMode(0.0,kmid);
	      for (diffm=0.0, i=Z20; i<=Today; i++)
		{
		  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
		  diffm+=vector[i]-GrowingMode(1./Time-1.,kmid)/normGMm;
		}
	      diffm/=(double)NBINS;

	      /* update of mindiff if necessary */
	      mindiff=(fabs(diffm) < mindiff ? fabs(diffm) : mindiff);

	      /* new extreme points */
	      if (diff1*diffm>0)
		{
		  logk1=logkmid;
		  diff1=diffm;
		}
	      else
		{
		  logk2=logkmid;
		  diff2=diffm;
		}
	      ++iter;
	    }
	  while (fabs(mindiff)>SMALLDIFF && iter<=MAXITER);

	  Smoothing.k_GM_displ[ismooth]=kmid;
	}

#ifdef DEBUG
      printf("displacements, smoothing %d, iter=%d\n",ismooth,iter);
      Time=1.0;
      gsl_integration_qags (&Function2, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGM1=sqrt(result);
      gsl_integration_qags (&Function31, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGM2=sqrt(result);
      gsl_integration_qags (&Function32, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGMm=sqrt(result);
      if (!ThisTask)
	for (i=0; i<NBINS; i++)
	  {
	    Time=pow(10.,SPLINE[SP_TIME]->x[i]);
	    fprintf(fd,"  %d %g  %g  %g  %g  %g",
		    ismooth,Smoothing.Rad_GM[ismooth],Smoothing.k_GM_displ[ismooth],
		    Time,vector[i],GrowingMode(1./Time-1.,Smoothing.k_GM_displ[ismooth])/GrowingMode(0.0,Smoothing.k_GM_displ[ismooth]));
	    gsl_integration_qags (&Function2, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g",sqrt(result)/normGM1,GrowingMode_2LPT(1./Time-1.,Smoothing.k_GM_displ[ismooth])/GrowingMode_2LPT(0.0,Smoothing.k_GM_displ[ismooth]));
	    gsl_integration_qags (&Function31, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g",sqrt(result)/normGM2,GrowingMode_3LPT_1(1./Time-1.,Smoothing.k_GM_displ[ismooth])/GrowingMode_3LPT_1(0.0,Smoothing.k_GM_displ[ismooth]));
	    gsl_integration_qags (&Function32, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g\n",sqrt(result)/normGMm,GrowingMode_3LPT_2(1./Time-1.,Smoothing.k_GM_displ[ismooth])/GrowingMode_3LPT_2(0.0,Smoothing.k_GM_displ[ismooth]));
	  }
#endif

    }

#ifdef DEBUG
  if (!ThisTask)
    fclose(fd);
#endif

  /**************/
  /* velocities */
  /**************/
  Function.function = &IntegrandForSDVelVariance;
#ifdef DEBUG
  Function2.function = &IntegrandForSDVel2Variance;
  Function31.function = &IntegrandForSDVel31Variance;
  Function32.function = &IntegrandForSDVel32Variance;
#endif

  /* displacement variance is computed using Top-Hat smoothing, 
     the integral will be again truncated at the Nyquist frequency */
  WindowFunctionType=2;

#ifdef DEBUG
  if (!ThisTask)
    fd=fopen("accuracy_scaledependent_velo.txt","w");
#endif

  /* cycle on smoothing radii */
  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {
      ThisRadius=Smoothing.Rad_GM[ismooth];

      /* integrates the growing mode and stores in vector */
      for (i=0; i<NBINS; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
  	  gsl_integration_qags (&Function, -4., nyquist, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  vector[i]=sqrt(result);
	}
      normv=vector[Today];
      for (i=0; i<NBINS; i++)
	vector[i]/=normv;

      /* bisector search of the k that gives the best approximation to the growth rate */
      iter=0;

      /* first guess at the extremes in k */
      logk1=LOGKMIN;
      logk2=LOGKMIN+(NkBINS-1)*DELTALOGK;
      k1=pow(10.,logk1);
      k2=pow(10.,logk2);

      /* differences at the extremes in k */
      normGM1=GrowingMode(0.0,k1)*fomega(0.0,k1);
      normGM2=GrowingMode(0.0,k2)*fomega(0.0,k2);
      for (diff1=diff2=0.0, i=Z20; i<=Today; i++)
	{
	  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
	  diff1+=vector[i]-GrowingMode(1./Time-1.,k1)*fomega(1./Time-1,k1)/normGM1;
	  diff2+=vector[i]-GrowingMode(1./Time-1.,k2)*fomega(1./Time-1,k1)/normGM2;
	}
      diff1/=(double)NBINS;
      diff2/=(double)NBINS;

      /* k1 might be accurate enough */
      if (fabs(diff1)<SMALLDIFF)
	Smoothing.k_GM_vel[ismooth]=k1;

      /* k2 might be accurate enough */
      else if (fabs(diff2)<SMALLDIFF)
	Smoothing.k_GM_vel[ismooth]=k2;

      /* this is an unfortunate case: the two differences have the same sign
         and none is accurate enough. In this case the scale is set to the k that gives
	 the smallest difference, and a warning is printed out  */
      else if (diff1*diff2>0)
	{
	  if (!ThisTask)
	    printf("WARNING in scale-dependent density growth rate for smoothing radius %d (%f): accuracy not guaranteed [diff1=%g - diff2=%g]\n",
		   ismooth, Smoothing.Rad_GM[ismooth],diff1,diff2);
	  if (fabs(diff1)<fabs(diff2))
	    Smoothing.k_GM_vel[ismooth]=k1;
	  else
	    Smoothing.k_GM_vel[ismooth]=k2;
	}
      else
	{
	  /* smallest difference */
	  mindiff=fabs(diff1);
	  mindiff=(fabs(diff2) < mindiff ? fabs(diff2) : mindiff);

	  /* iteration for bisector search */
	  do
	    {

	      /* midpoint in log space */
	      logkmid=0.5*(logk1+logk2);
	      kmid=pow(10.,logkmid);

	      /* difference at midpoint */
	      normGMm=GrowingMode(0.0,kmid)*fomega(0.0,kmid);
	      for (diffm=0.0, i=Z20; i<=Today; i++)
		{
		  Time = pow(10.,SPLINE[SP_TIME]->x[i]);
		  diffm+=vector[i]-GrowingMode(1./Time-1.,kmid)*fomega(1./Time-1.,kmid)/normGMm;
		}
	      diffm/=(double)NBINS;

	      /* update of mindiff if necessary */
	      mindiff=(fabs(diffm) < mindiff ? fabs(diffm) : mindiff);

	      /* new extreme points */
	      if (diff1*diffm>0)
		{
		  logk1=logkmid;
		  diff1=diffm;
		}
	      else
		{
		  logk2=logkmid;
		  diff2=diffm;
		}
	      ++iter;
	    }
	  while (fabs(mindiff)>SMALLDIFF && iter<=MAXITER);

	  Smoothing.k_GM_vel[ismooth]=kmid;
	}

#ifdef DEBUG
      printf("velocities, smoothing %d, iter=%d\n",ismooth,iter);
      Time=1.0;
      gsl_integration_qags (&Function2, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGM1=sqrt(result);
      gsl_integration_qags (&Function31, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGM2=sqrt(result);
      gsl_integration_qags (&Function32, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      normGMm=sqrt(result);
      if (!ThisTask)
	for (i=0; i<NBINS; i++)
	  {
	    Time=pow(10.,SPLINE[SP_TIME]->x[i]);
	    fprintf(fd,"  %d %g  %g  %g  %g  %g",
		    ismooth,Smoothing.Rad_GM[ismooth],Smoothing.k_GM_vel[ismooth],
		    Time,vector[i],GrowingMode(1./Time-1.,Smoothing.k_GM_vel[ismooth])*fomega(1./Time-1.,Smoothing.k_GM_vel[ismooth])/
		    GrowingMode(0.,Smoothing.k_GM_vel[ismooth]) / fomega(0.,Smoothing.k_GM_vel[ismooth]));

	    gsl_integration_qags (&Function2, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g",sqrt(result)/normGM1,GrowingMode_2LPT(1./Time-1.,Smoothing.k_GM_vel[ismooth])*fomega_2LPT(1./Time-1.,Smoothing.k_GM_vel[ismooth])/
		    GrowingMode_2LPT(0.,Smoothing.k_GM_vel[ismooth]) /fomega_2LPT(0.0,Smoothing.k_GM_vel[ismooth]));

	    gsl_integration_qags (&Function31, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g",sqrt(result)/normGM2,GrowingMode_3LPT_1(1./Time-1.,Smoothing.k_GM_vel[ismooth])*fomega_3LPT_1(1./Time-1.,Smoothing.k_GM_vel[ismooth])/
		    GrowingMode_3LPT_1(0.,Smoothing.k_GM_vel[ismooth])/fomega_3LPT_1(0.,Smoothing.k_GM_vel[ismooth]));

	    gsl_integration_qags (&Function32, -4., 2., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	    fprintf(fd,"  %g  %g\n",sqrt(result)/normGMm,GrowingMode_3LPT_2(1./Time-1.,Smoothing.k_GM_vel[ismooth])*fomega_3LPT_2(1./Time-1.,Smoothing.k_GM_vel[ismooth])/
		    GrowingMode_3LPT_2(0.,Smoothing.k_GM_vel[ismooth])/fomega_3LPT_2(0.,Smoothing.k_GM_vel[ismooth]));

	  }
#endif

    }

#ifdef DEBUG
  if (!ThisTask)
    fclose(fd);
#endif

  free(vector);

  WindowFunctionType=2;
  return 0;
}


#endif


void greetings(void)
{
  /* This is a list of messages to declare the most relevant precompiler directives in the stdout */

  if (!ThisTask)
    {
      printf("[%s] This is pinocchio V5.0, running on %d MPI tasks\n\n",fdate(),NTasks);
#ifdef _OPENMP
      printf( "Using %d OpenMP threads\n", internal.nthreads_omp );
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
#ifdef NORADIATION
      printf("Radiation is not included in the Friedmann equations\n");
#else
      printf("Radiation is included in the Friedmann equations\n");
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
#error WHITENOISE is not implemented yet
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
