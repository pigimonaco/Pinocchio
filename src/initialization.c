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

int init_cosmology(void);
int set_smoothing(void);
int generate_densities(void);
int set_subboxes(void);
int set_plc(void);
unsigned int gcd(unsigned int, unsigned int);
int          set_fft_decomposition(void);


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
  /* reads the P(k) tables from CAMB */
  if (read_power_table_from_CAMB())
    return 1;
#endif

  /* call to cosmo.c: initialization of cosmological functions */
  if (initialize_cosmology())
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

#ifdef SCALE_DEPENDENT_GROWTH
  /* initialize scale-dependent growing modes and fomegas */
  if (initialize_ScaleDependentGrowth())
    return 1;
#endif

  /* computes the number of sub-boxes for fragmentation */
  SET_WTIME;  
  if (set_subboxes())
    {
      if (!ThisTask)
	printf("Pinocchio done!\n");
      MPI_Finalize();
      exit (0);
    }
  ASSIGN_WTIME(partial, set_subboxes);  

  /* initializes quantities needed for the on-the-fly reconstruction of PLC */
  SET_WTIME;  
  if (set_plc())
    return 1;
  ASSIGN_WTIME(partial, set_plc);  

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
  if (initialize_fft())
    return 1;
  ASSIGN_WTIME(partial, fft_initialization);
  
  /* generation of initial density field */
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

  params.ParticleMass = 2.775499745e11 * params.Hubble100 *
    params.Hubble100 * params.Omega0 * pow(params.InterPartDist,3.);
  
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
  double var_min, var_max, rmin;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=-1;   /* here we want to use the standard, scale-independent growing mode */
#endif

  var_min    = 1.686 / NSIGMA / GrowingMode( outputs.zlast );
  var_min   *= var_min;  
  /* rmax       = Radius(var_min); */
  rmin       = params.InterPartDist / 6.0;
  var_max    = MassVariance(rmin);
  Smoothing.Nsmooth = (log10(var_max) - log10(var_min)) / STEP_VAR + 2;

  if ( Smoothing.Nsmooth <=0 )
    {
      if (!ThisTask)
	dprintf( VERR, 0, "I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
      Smoothing.Nsmooth=1;
    }

  if ( !ThisTask )
    dprintf( VMSG, 0, "Min variance: %f12.6, max variance: %f12.6, number of smoothing radii: %d\n",
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

  int ismooth;
  for ( ismooth = 0; ismooth < Smoothing.Nsmooth-1; ismooth++ )
    {
      Smoothing.Variance[ismooth] = pow(10., log10(var_min)+STEP_VAR*ismooth);
      Smoothing.Radius[ismooth]   = Radius(Smoothing.Variance[ismooth]);
    }
  Smoothing.Radius[ismooth]   = 0.0;
  Smoothing.Variance[ismooth] = var_max;

  if (!ThisTask)
    for ( ismooth = 0; ismooth < Smoothing.Nsmooth; ismooth++ )
      dprintf(VMSG, 0, "           %2d)  Radius=%10f, Variance=%10f\n",
	      ismooth+1,Smoothing.Radius[ismooth],Smoothing.Variance[ismooth]);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


int generate_densities()
{

  cputime.dens = MPI_Wtime();

  if ( !ThisTask )
    dprintf( VMSG, 0, "[%s] Generating density in Fourier space\n",fdate());

#ifdef WHITENOISE

  if ( Ngrids > 1 )
    {
      if ( !ThisTask )
	dprintf( VMSG, 0, "Sorry, this works only with a single grid\n");
      return 1;
    }

  if ( read_white_noise() )
    return 1;

#else

  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    if (GenIC_large(igrid))
      return 1;

#endif

  cputime.dens = MPI_Wtime() - cputime.dens;
    if ( !ThisTask )
      dprintf(VMSG, 0, "[%s] Done generating density in Fourier space, cputime = %f s\n",fdate(), cputime.dens);

  return 0;
}


int set_grids()
{
  /* initialization of fftw quantities on grids (one for the moment) */
  Ngrids = 1;

  MyGrids=(grid_data*)malloc(Ngrids * sizeof(grid_data));

  MyGrids[0].GSglobal[_x_] = params.GridSize[0];
  MyGrids[0].GSglobal[_y_] = params.GridSize[1];
  MyGrids[0].GSglobal[_z_] = params.GridSize[2];
  
  MyGrids[0].Ntotal = (unsigned long long)MyGrids[0].GSglobal[_x_] * 
    (unsigned long long)MyGrids[0].GSglobal[_y_] * 
    (unsigned long long)MyGrids[0].GSglobal[_z_];

  MyGrids[0].BoxSize        = params.BoxSize_htrue;
  MyGrids[0].lower_k_cutoff = 0.0;
  MyGrids[0].upper_k_cutoff = NYQUIST * PI;

  /* allocates pointers */
  cvector_fft        = (fftw_complex**) malloc(Ngrids * sizeof(fftw_complex*));
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
  // seedtable=(unsigned int**)malloc(Ngrids * sizeof(unsigned int*));
 
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
  int    NAll, ic, intersection, axis;
  double Largest_r, Smallest_r, smallestr, largestr, x[3], l[3], d, mod;
  FILE  *fout;
  char   filename[BLENGTH];

  /* ordering of coordinates to accomodate for rotation caused by fft ordering */
#ifdef ROTATE_BOX
  static int rot[3] = {1,2,0};
#else
  static int rot[3] = {0,1,2};
#endif

  if ( params.StartingzForPLC < 0.0 )
    {
      plc.Nreplications    = 0;
      plc.Fstart=plc.Fstop = -1.0;
      plc.Nmax             = 0;
      if ( !ThisTask )
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
      plc.center[0] = gsl_rng_uniform(random_generator) * MyGrids[0].GSglobal[_x_];
      plc.center[1] = gsl_rng_uniform(random_generator) * MyGrids[0].GSglobal[_y_];
      plc.center[2] = gsl_rng_uniform(random_generator) * MyGrids[0].GSglobal[_z_];
      plc.zvers[0]  = 1.0;
      plc.zvers[1]  = 1.0;
      plc.zvers[2]  = 1.0;
    }
  
  /* normalization of the cone axis direction */
  mod = sqrt(plc.zvers[0]*plc.zvers[0] + plc.zvers[1]*plc.zvers[1] + plc.zvers[2]*plc.zvers[2]);
  for ( int ir = 0; ir < 3; ir++ )
    plc.zvers[ir] /= mod;

  /* here we define a system where zvers is the z-axis */
  if ( plc.zvers[2] == 1.0 )
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
      mod = sqrt( plc.zvers[0]*plc.zvers[0] + plc.zvers[1]*plc.zvers[1] );
      plc.xvers[0] =  plc.zvers[1] / mod;
      plc.xvers[1] = -plc.zvers[0] / mod;
      plc.xvers[2] =  0.0;
      /* y axis will be the cross product of z and x */
      plc.yvers[0] = plc.zvers[1]*plc.xvers[2] - plc.zvers[2]*plc.xvers[1];
      plc.yvers[1] = plc.zvers[2]*plc.xvers[0] - plc.zvers[0]*plc.xvers[2];
      plc.yvers[2] = plc.zvers[0]*plc.xvers[1] - plc.zvers[1]*plc.xvers[0];
    }

  /* initialization to compute the number of realizations */
  NAll = (int)(ProperDistance(params.StartingzForPLC) / MyGrids[0].BoxSize) + 2;
  plc.Fstart = 1.0 + params.StartingzForPLC;
  plc.Fstop  = 1.0 + params.LastzForPLC;
  /* if F is the inverse collapse time: */
  /* plc.Fstart = 1./GrowingMode(params.StartingzForPLC); */
  /* plc.Fstop  = 1./GrowingMode(params.LastzForPLC); */

  Largest_r  = ProperDistance(params.StartingzForPLC) / params.InterPartDist;
  Smallest_r = ProperDistance(params.LastzForPLC) / params.InterPartDist;

  l[0] = (double)(MyGrids[0].GSglobal[_x_]);
  l[1] = (double)(MyGrids[0].GSglobal[_y_]);
  l[2] = (double)(MyGrids[0].GSglobal[_z_]);

  /* first, it counts the number of replications needed */
  plc.Nreplications = 0;
  for ( int ir = -NAll; ir <= NAll; ir++ )
    for ( int jr = -NAll; jr <= NAll; jr++ )
      for ( int kr = -NAll; kr <= NAll; kr++ )
	{
	  x[0] = ir*l[0];
	  x[1] = jr*l[1];
	  x[2] = kr*l[2];

	  intersection = cone_and_cube_intersect(x, l, plc.center, plc.zvers, params.PLCAperture, &smallestr, &largestr, &axis);

	  if (intersection && !(smallestr>Largest_r || largestr<Smallest_r))
	    plc.Nreplications++;
	}

  /* second, it allocates the needed memory for the replications */
  plc.repls = malloc(plc.Nreplications * sizeof(replication_data));

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
  {
    int this = 0;
    for ( int ir = -NAll; ir <= NAll; ir++ )
      for ( int jr = -NAll; jr <= NAll; jr++ )
	for ( int kr = -NAll; kr <= NAll; kr++ )
	  {
	    x[0] = ir*l[0];
	    x[1] = jr*l[1];
	    x[2] = kr*l[2];
	    
	    intersection = cone_and_cube_intersect(x, l, plc.center, plc.zvers, params.PLCAperture, &smallestr, &largestr, &axis);
	    
	    if ( intersection && !(smallestr>Largest_r || largestr<Smallest_r) )
	      {
		if (!ThisTask)
		  fprintf(fout," %3d  %3d %3d %3d   %10.6f %10.6f   %d  %d\n",this,ir,jr,kr,smallestr,largestr,intersection,axis);
		plc.repls[this].i = ir;
		plc.repls[this].j = jr;
		plc.repls[this].k = kr;
		plc.repls[this].F1   = -largestr;
		plc.repls[this++].F2 = -smallestr;
	      }	  
	  }
    if ( !ThisTask )
      fclose(fout);
  }

  /* fourth it transforms distances to redshifts */
  for ( double z = 100.; z >= 0.0; z -= 0.01 )
    {
      d = ProperDistance(z) / params.InterPartDist;
      for ( int this = 0; this < plc.Nreplications; this++ )
	{
	  if ( plc.repls[this].F1 <= 0.0 && d<-plc.repls[this].F1 )
	    plc.repls[this].F1 =z + 0.01 + 1.0;
	  if ( plc.repls[this].F2 <= 0.0 && d<-plc.repls[this].F2 )
	    plc.repls[this].F2 = z - 0.01 + 1.0;
	}
    }
  for ( int this = 0; this < plc.Nreplications; this++ )
    {
      if ( plc.repls[this].F1 <= 0.0 )
	plc.repls[this].F1 = 1.0;
      if ( plc.repls[this].F2 <= 0.0 )
	plc.repls[this].F2 = 1.0;
    }

  plc.Nmax = subbox.Npart / 10;

  if ( !ThisTask )
    {
      dprintf( VMSG, 0, "The Past Light Cone will be reconstruct from z=%f to z=%f\n",
	     params.StartingzForPLC,params.LastzForPLC);
      if (params.PLCProvideConeData)
	dprintf( VMSG, 0, "Cone data have been provided in the parameter file\n");
      else
	dprintf( VMSG, 0, "Cone data have been decided by the code\n");
      dprintf( VMSG, 0, "Past Light Cone will be centred on point [%f,%f,%f] (true Mpc)\n",
	     plc.center[0]*params.InterPartDist,plc.center[1]*params.InterPartDist,plc.center[2]*params.InterPartDist);
      dprintf( VMSG, 0, "The cone vertex will be pointed toward [%f,%f,%f]\n",plc.zvers[0],plc.zvers[1],plc.zvers[2]);
      dprintf( VMSG, 0, "It will have an aperture of %f degrees\n",params.PLCAperture);
#ifdef ROTATE_BOX
      if (params.PLCProvideConeData)
	dprintf( VMSG, 0, "(NB: rotation has been applied to the provided coordinates)\n");
#endif
      dprintf( VMSG, 0, "The proper distance at the starting redshift, z=%f, is: %f Mpc\n",
	     params.StartingzForPLC, Largest_r*params.InterPartDist);
      dprintf( VMSG, 0, "The proper distance at the stopping redshift, z=%f, is: %f Mpc\n",
	     params.LastzForPLC, Smallest_r*params.InterPartDist);
      dprintf( VMSG, 0, "The reconstruction will be done for %f < z < %f\n",params.LastzForPLC,params.StartingzForPLC);
      dprintf( VMSG, 0, "The corresponding F values are: Fstart=%f, Fstop=%f\n",plc.Fstart,plc.Fstop);
      dprintf( VMSG, 0, "The box will be replicated %d times to construct the PLC\n",plc.Nreplications);
      for (ic=0; ic<plc.Nreplications; ic++)
	dprintf( VMSG, 0, "   Replication %2d: shift (%2d,%2d,%2d), from F=%f to F=%f\n",
	       ic,plc.repls[ic].i,plc.repls[ic].j,plc.repls[ic].k,
	       plc.repls[ic].F1,plc.repls[ic].F2);
      dprintf( VMSG, 0, "Task 0 will use plc.Nmax=%d\n",plc.Nmax);
      dprintf( VMSG, 0, "\n");
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

  double dP = sqrt( (P[0]-V[0]) * (P[0]-V[0]) +
		    (P[1]-V[1]) * (P[1]-V[1]) +
		    (P[2]-V[2]) * (P[2]-V[2]));
  if (dP == 0.0)
    return 1.0;
  
  double cosDU = ( D[0]*U[0] + D[1]*U[1] + D[2]*U[2]);
  double cosDP = ( D[0]*(P[0]-V[0]) + D[1]*(P[1]-V[1]) + D[2]*(P[2]-V[2]) ) / dP;
  double cosUP = ( U[0]*(P[0]-V[0]) + U[1]*(P[1]-V[1]) + U[2]*(P[2]-V[2]) ) / dP;

  if ( cosDP-cosDU*cosUP ==0.0 )
    return 0.0;
  double tmax = (cosDU-cosDP*cosUP) / (cosDP-cosDU*cosUP);
  if ( tmax < 0 )
    tmax = 0.0;
  else if ( tmax > *L / dP )
    tmax = *L / dP;

  return ( cosDP + tmax*cosDU) / sqrt(1.0 + tmax*tmax + 2.0*tmax*cosUP);
}

int cone_and_cube_intersect(double *Oc, double *L, double *V, double *D, double theta, double *rmin, double *rmax, int *axis)
{
  int    ivec[3];
  double F, Fmax, costh, U[3], P[3];

  /* This routine returns 1 if the cone with vertex V, axis direction D and semi-aperture theta (deg)
     intersects the cube starting from point Oc and with edges of lenght L aligned with the axes.
     It returns 0 if the two solids do not intersect.
     It also computes the smallest and largest distances of the cube from the cone vertex V.
   */


  /* Computation of rmin and rmax */
  {
    double rmax2 = 0.0;
    /* max distance from vertices */
    for (int i = 0; i < 2; i++ )
      for (int j = 0; j < 2; j++ )
	for ( int k = 0; k < 2; k++ )
	  {
	    double r2 = (Oc[0]+i*L[0]-V[0])*(Oc[0]+i*L[0]-V[0]) +
	      (Oc[1]+j*L[1]-V[1]) * (Oc[1]+j*L[1]-V[1]) +
	      (Oc[2]+k*L[2]-V[2]) * (Oc[2]+k*L[2]-V[2]);
	    if (r2 > rmax2)
	      rmax2 = r2;
	  }
    
    *rmax = sqrt(rmax2);
  }

  
  /* min distance from cube faces */
  /* and intersection of axis with cube faces */
  {
    double _axis_ = 0;
    double rmin2 = 1e32;
    for (int dim = 0; dim < 3; dim++)  /* three dimension (normals to cube faces) */
      for (int i = 0; i < 2; i++)      /* two faces per dimension */
	{
	  double proj = Oc[dim] - V[dim] + i*L[dim];
	  int dim1 = (dim+1) % 3;
	  int dim2 = (dim+2) % 3;
	  
	  /* minimum distance */
	  double r2 = proj*proj;                /* the normal component always contributes */
	  if (V[dim1] < Oc[dim1])               /* only of their projection is outside the face */
	    r2 += (V[dim1]-Oc[dim1])*(V[dim1]-Oc[dim1]);
	  else if ( V[dim1] > Oc[dim1]+L[dim1] )
	    r2 += (V[dim1]-Oc[dim1]-L[dim1]) * (V[dim1]-Oc[dim1]-L[dim1]);
	  if ( V[dim2] < Oc[dim2] )
	    r2 += (V[dim2]-Oc[dim2])*(V[dim2]-Oc[dim2]);
	  else if ( V[dim2] > Oc[dim2]+L[dim2] )
	    r2 += (V[dim2]-Oc[dim2]-L[dim2])*(V[dim2]-Oc[dim2]-L[dim2]);

	  if (r2 < rmin2)
	    rmin2 = r2;
	  
	  /* axis intersection */
	  double x;
	  if ( (x = proj / D[dim]) > 0.0 &&
	       V[dim1] + x*D[dim1] >= Oc[dim1] &&
	       V[dim1] + x*D[dim1] <= Oc[dim1] + L[dim1] &&
	       V[dim2] + x*D[dim2] >= Oc[dim2] &&
	       V[dim2] + x*D[dim2] <= Oc[dim2] + L[dim2] )
	    _axis_ += 1<<(dim+i*3);
	}

    *axis = _axis_;
    *rmin = sqrt(rmin2);
  }

  /* step 1: if the whole sky is required, only rmin and rmax are needed */
  if ( theta >= 180. )
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
  Fmax  = -10.0;
  costh = cos( theta / 180. * PI );
  for ( int i = 0; i < 2; i++ )
    for ( int j = 0; j < 2; j++ )
      for ( int k = 0; k < 2; k++ )
	{
	  ivec[0] = i;
	  ivec[1] = j;
	  ivec[2] = k;
	  
	  for ( int dim = 0; dim < 3; dim++ )
	    if ( !ivec[dim] )
	      {
		U[dim]       = 1.0;
		U[(dim+1)%3] = 0.0;
		U[(dim+2)%3] = 0.0;
		P[0]         = Oc[0] + ivec[0]*L[0];
		P[1]         = Oc[1] + ivec[1]*L[1];
		P[2]         = Oc[2] + ivec[2]*L[2];
		F            = maxF(P, V, U, D, L+dim)-costh;
		if (F > Fmax)
		  Fmax = F;
	      }
	}

  /* if the nearest vertex is inside the cone then exit */
  if ( Fmax > 0 )
    return 4;

  /* at this point the cone and the cube do not intersect */
  return 0;

}


#else

int set_plc()
{
  if (!ThisTask)
    dprintf( VERR, 0, "PLC flag at compilation was not set, no Past Light Cone output will be given\n\n");

  return 0;
}

#endif


/* division in sub-boxes */
int set_subboxes()
{

  int          i1,j1,k1, surface, tt;
  int          NN, NN1, N1, N2, N3;
  double       size,sizeG,cc;
  unsigned int TotalNP_pertask;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=-1;
#endif

  /* mass of the largest halo expected in the box */
  params.Largest = 1.e18;
  cc             = 1.0/pow(params.BoxSize_htrue,3.0);
  double aa      = AnalyticMassFunction(params.Largest,outputs.zlast);
  while (aa*params.Largest < cc)
    {
      params.Largest *= 0.99;
      aa              = AnalyticMassFunction(params.Largest,outputs.zlast);
    }
  size  = SizeForMass(params.Largest);
  sizeG = size / params.InterPartDist;

  /*  
      The number of loadable subbox particles is equal 
      to the number allowed by the specified MaxMemPerParticle.
      The boundary layer is set to its maximum value.
  */

  /* finds the optimal number of sub-boxes to use for the fragmentation */
  TotalNP_pertask = (unsigned int)(MyGrids[0].Ntotal/(unsigned long long)NTasks);
  surface         = TotalNP_pertask;
  NSlices         = 1;   // Qui ho spento il calcolo del numero di slice

  double fact = 2*sizeG * 2*sizeG;
  for ( int k = 1; k <= NTasks; k++ )
    for ( int j = 1; j <= NTasks/k; j++ )
      for ( int i = 1; i <= NTasks/k/j; i++ )
	/* the three indices must be exact divisors of the three grid lengths */
	if ( i*j*k == NTasks)
	  {
	    /* number of particles in the sub-box */
	    N1 = find_length(MyGrids[0].GSglobal[_x_],i,0);
	    N2 = find_length(MyGrids[0].GSglobal[_y_],j,0);
	    N3 = find_length(MyGrids[0].GSglobal[_z_],k*NSlices,0);

	    int this = (i>1? 2*(N2*N3) : 0) + 
	      (j>1? 2*(N1*N3) : 0) +
	      (k>1? 2*(N1*N2) : 0);
	    tt = this;
	    if (N1/2 < sizeG)
	      this += (int)((double)tt * fact / (N1*N1));
	    if (N2/2 < sizeG)
	      this += (int)((double)tt * fact / (N2*N2));
	    if (N3/2 < sizeG)
	      this += (int)((double)tt * fact / (N3*N3));

	    if (this < surface)
	      {
		surface = this;
		i1      = i; 
		j1      = j; 
		k1      = k; 
	      }
	  }

  subbox.nbox_x = i1;
  subbox.nbox_y = j1;
  subbox.nbox_z_thisslice = k1;
  subbox.nbox_z_allslices = k1*NSlices;

  /* this will be mybox for the first slice */
			    
  NN             = subbox.nbox_y*subbox.nbox_z_thisslice;
  subbox.mybox_x = ThisTask/NN;
  NN1            = ThisTask-subbox.mybox_x*NN;
  subbox.mybox_y = NN1/subbox.nbox_z_thisslice;
  subbox.mybox_z = NN1-subbox.mybox_y*subbox.nbox_z_thisslice;

  subbox.Lgrid_x = find_length(MyGrids[0].GSglobal[_x_],subbox.nbox_x,subbox.mybox_x);
  subbox.Lgrid_y = find_length(MyGrids[0].GSglobal[_y_],subbox.nbox_y,subbox.mybox_y);
  subbox.Lgrid_z = find_length(MyGrids[0].GSglobal[_z_],subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.pbc_x = (subbox.nbox_x==1);
  subbox.pbc_y = (subbox.nbox_y==1);
  subbox.pbc_z = (subbox.nbox_z_allslices==1);

  subbox.safe_x = (subbox.pbc_x ? 0 : (find_length(MyGrids[0].GSglobal_x,subbox.nbox_x,0)-1)/2);
  subbox.safe_y = (subbox.pbc_y ? 0 : (find_length(MyGrids[0].GSglobal_y,subbox.nbox_y,0)-1)/2);
  subbox.safe_z = (subbox.pbc_z ? 0 : (find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,0)-1)/2);

  subbox.Lgwbl_x = subbox.Lgrid_x + 2*subbox.safe_x;
  subbox.Lgwbl_y = subbox.Lgrid_y + 2*subbox.safe_y;
  subbox.Lgwbl_z = subbox.Lgrid_z + 2*subbox.safe_z;

  subbox.start_x = find_start(MyGrids[0].GSglobal[_x_],subbox.nbox_x,subbox.mybox_x);
  subbox.start_y = find_start(MyGrids[0].GSglobal[_y_],subbox.nbox_y,subbox.mybox_y);
  subbox.start_z = find_start(MyGrids[0].GSglobal[_z_],subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.stabl_x = subbox.start_x - subbox.safe_x;
  subbox.stabl_y = subbox.start_y - subbox.safe_y;
  subbox.stabl_z = subbox.start_z - subbox.safe_z;

  /* 
     Npart: total number of particles in the whole sub-volume 
     Ngood: total number of particles in the well reconstructed region
     Npredpeaks: a guess of the maximum number of peaks in the subbox
     Nalloc: number of particles for which memory has been allocated (set in allocate_main_memory)
     Nstored: number of actually stored particles
  */

  subbox.Npart = subbox.Lgwbl_x * subbox.Lgwbl_y * subbox.Lgwbl_z;
  subbox.Ngood = subbox.Lgrid_x * subbox.Lgrid_y * subbox.Lgrid_z;
  subbox.PredNpeaks = subbox.Ngood/5;
  subbox.Nstored=0;
  /* this is the size of frag_map*/
  subbox.maplength = subbox.Npart/UINTLEN + (subbox.Npart%UINTLEN!=0);
  subbox.Nalloc = set_main_memory(TotalNP_pertask);

  /* messagges */
  if (!ThisTask)
    {
      dprintf( VMSG, 0, "\n");
      dprintf( VMSG, 0, "FRAGMENTATION:\n");
      dprintf( VMSG, 0, "Number of sub-boxes per dimension:     %d %d %d\n",subbox.nbox_x,subbox.nbox_y,subbox.nbox_z_allslices);
      dprintf( VMSG, 0, "Periodic boundary conditions:          %d %d %d\n",subbox.pbc_x,subbox.pbc_y,subbox.pbc_z);
      dprintf( VMSG, 0, "Core 0 will work on a grid:            %d %d %d\n",subbox.Lgwbl_x,subbox.Lgwbl_y,subbox.Lgwbl_z);
      dprintf( VMSG, 0, "The resolved box will be:              %d %d %d\n",subbox.Lgrid_x,subbox.Lgrid_y,subbox.Lgrid_z);
      dprintf( VMSG, 0, "Boundary layer:                        %d %d %d\n",subbox.safe_x,subbox.safe_y,subbox.safe_z);
      dprintf( VMSG, 0, "Number of total particles for core 0:  %d\n",subbox.Npart);
      dprintf( VMSG, 0, "Number of good particles for core 0:   %d\n",subbox.Ngood);
      dprintf( VMSG, 0, "Particles that core 0 will allocate:   %d\n",subbox.Nalloc);
      dprintf( VMSG, 0, "Largest allowed overhead:              %f\n",(float)subbox.Nalloc/(float)subbox.Ngood);
      dprintf( VMSG, 0, "Best number of surface particles: %d %f\n",surface,(float)surface/(float)TotalNP_pertask);  // POI SI LEVA
      dprintf( VMSG, 0, "Largest halo expected in this box at z=%f: %e Msun\n",
  	     outputs.zlast, params.Largest);
      dprintf( VMSG, 0, "   its Lagrangian size: %f Mpc (%6.2f grid points)\n",size,sizeG);
      if ((!subbox.pbc_x && sizeG>subbox.safe_x) || 
	  (!subbox.pbc_y && sizeG>subbox.safe_y) || 
	  (!subbox.pbc_z && sizeG>subbox.safe_z))
	{
	  dprintf( VMSG, 0, "WARNING: the boundary layer on some dimension smaller than the predicted size of the largest halos,\n");
	  dprintf( VMSG, 0, "         the most massive halos may be inaccurate\n");
	}
//      dprintf( VMSG, 0, "Allowed overhead: %f\n",subbox.overhead);
    }


  if (subbox.Nalloc < 0)
    {
      dprintf( VXERR, 0, "ERROR on Task %d: negative Nalloc, please increase MaxMemPerParticle\n",ThisTask);
      return 1;
    }

  // AGGIUNGERE CHECK SU OVERHEAD, USCIRE SE E` <1

  /* NSlices>1 is incompatible with WriteSnapshot */
  if (NSlices > 1 && params.WriteSnapshot)
    {
      params.WriteSnapshot=0;
      if (!ThisTask)
	dprintf( VERR, 0, "Sorry, but snapshots cannot be written if fragmentation is done in slices\n");
    }

  /* initialization of quantities required by compute_mf */
   if (params.OutputInH100)
    mf.hfactor = params.Hubble100;
  else
    mf.hfactor = 1.0;
  mf.hfactor4  = pow(mf.hfactor,4.);
  mf.vol       = (double)MyGrids[0].GSglobal[_x_]*(double)MyGrids[0].GSglobal[_y_] * 
    (double)MyGrids[0].GSglobal[_z_]*pow(params.InterPartDist,3.0);
  mf.mmin      = log10(params.MinHaloMass*params.ParticleMass)-0.001*DELTAM;
  mf.mmax      = log10(params.Largest)+3.0*DELTAM;
  mf.NBIN      = (int)((mf.mmax-mf.mmin)/DELTAM) +1;
  mf.ninbin    = (int*)calloc(mf.NBIN,sizeof(int));
  mf.ninbin_local    = (int*)calloc(mf.NBIN,sizeof(int));
  mf.massinbin = (double*)calloc(mf.NBIN,sizeof(double));
  mf.massinbin_local = (double*)calloc(mf.NBIN,sizeof(double));

  if (!ThisTask)
    {
      dprintf( VMSG, 0, "\nThe mass function will be computed from Log M=%f to Log M=%f (%d bins)\n",
      	     mf.mmin, mf.mmax, mf.NBIN);
      dprintf( VMSG, 0, "HO ABBASSATO DELTAM A 0.01, RIPORTARLO A 0.05!!!\n");
      dprintf( VMSG, 0, "\n");
      fflush(stdout);      
    }

  return 0;
}


int find_start(int L,int n,int ibox)
{
  int LL,MM;

  if (n == 1)
    return 0;
  else
    {
      LL = L/n;
      MM = L%n;
      if ( ibox == 0 )
        return 0;
      else if ( ibox <= MM)
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

  if (n == 1) 
    return L;
  else
    {
      LL = L/n;
      MM = L%n;
      if ( ibox < MM )
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


