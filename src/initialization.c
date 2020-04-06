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

int init_cosmology(void);
int set_smoothing(void);
int generate_densities(void);
int set_subboxes(void);
int set_plc(void);
#ifdef SCALE_DEPENDENT
int set_scaledep_GM(void);
#endif

/* // LEVARE */
/* void test_inverseGM() */
/* { */
/*   double D; */
/*   FILE *fd=fopen("testinv.txt","w"); */
/*   for (int ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++) */
/*     /\* for (int i=0; i<NBINS; i++) *\/ */
/*     /\*   fprintf(fd," %d %d  %g %g\n", *\/ */
/*     /\* 	      ismooth,i,SPLINE_INVGROW[ismooth]->x[i],SPLINE_INVGROW[ismooth]->y[i]); *\/ */
/*     for (double z=100.; z>= 0.0; z-=1.) */
/*       { */
/*     	D=GrowingMode(z,params.k_for_GM); */
/*     	fprintf(fd," %f %d     %g  %g\n", */
/*     		z,ismooth,D,InverseGrowingMode(D,ismooth)); */
/*       } */
/*   fclose(fd); */
/* } */


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
  
  /* computes the number of sub-boxes for fragmentation */
  if (set_subboxes())
    {
      if (!ThisTask)
	printf("Pinocchio done!\n");
      MPI_Finalize();
      exit (0);
    }

#ifdef SCALE_DEPENDENT
  /* computes the growth rates for displacements */
  if (set_scaledep_GM())
    return 1;
#endif

  /* initializes quantities needed for the on-the-fly reconstruction of PLC */
  if (set_plc())
    return 1;

  /* this barrier is set to have correct stdout in case the run cannot start */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* allocations of memory for fmax and memory tests */
  if (allocate_main_memory())
    return 1;

  /* initialization of fft plans */
  if (initialize_fft())
    return 1;

  /* generation of initial density field */
  if (generate_densities())
    return 1;

  cputime.init=MPI_Wtime()-cputime.init;

  if (!ThisTask)
    printf("[%s] initialization done, cpu time = %14.6f\n", fdate(), cputime.init);

  return 0;
}


int set_parameters()
{
  int i;

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
  params.InterPartDist = params.BoxSize_htrue/params.GridSize[0];

  params.ParticleMass = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0 
    * pow(params.InterPartDist,3.);
  strcpy(params.DataDir,"Data/");

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
      printf("Flag for this run: %s\n\n",params.RunFlag);
      printf("PARAMETER VALUES from file %s:\n",params.ParameterFile);
      printf("Omega0                      %f\n",params.Omega0);
      printf("OmegaLambda                 %f\n",params.OmegaLambda);    
      printf("OmegaBaryon                 %f\n",params.OmegaBaryon);
      if (strcmp(params.TabulatedEoSfile,"no"))
	{
	  printf("Dark Energy EoS will be read from file %s\n",params.TabulatedEoSfile);
	}
      else
	{
	  printf("DE EoS parameters           %f %f\n",params.DEw0,params.DEwa);
	}

      printf("Hubble100                   %f\n",params.Hubble100);
      printf("Sigma8                      %f\n",params.Sigma8);
      printf("PrimordialIndex             %f\n",params.PrimordialIndex);
      printf("RandomSeed                  %d\n",params.RandomSeed);
      printf("OutputList                  %s\n",params.OutputList);
      printf("Number of outputs           %d\n",outputs.n);
      printf("Output redshifts           ");
      for (i=0; i<outputs.n; i++)
	printf(" %f ",outputs.z[i]);
      printf("\n");
      printf("GridSize                    %d %d %d\n",params.GridSize[0],params.GridSize[1],params.GridSize[2]);
      printf("BoxSize (true Mpc)          %f\n",params.BoxSize_htrue);
      printf("BoxSize (Mpc/h)             %f\n",params.BoxSize_h100);
      printf("Particle Mass (true Msun)   %g\n",params.ParticleMass);
      printf("Particle Mass (Msun/h)      %g\n",params.ParticleMass*params.Hubble100);
      printf("Inter-part dist (true Mpc)  %f\n",params.InterPartDist);
      printf("Inter-part dist (Mpc/h)     %f\n",params.InterPartDist*params.Hubble100);
      printf("MinHaloMass (particles)     %d\n",params.MinHaloMass);
      printf("BoundaryLayerFactor         %f\n",params.BoundaryLayerFactor);
      printf("MaxMem per task (Mb)        %d\n",params.MaxMem);
      printf("MaxMem per particle (b)     %f\n",params.MaxMemPerParticle);
      printf("CatalogInAscii              %d\n",params.CatalogInAscii);
      printf("NumFiles                    %d\n",params.NumFiles);
      printf("DoNotWriteCatalogs          %d\n",params.DoNotWriteCatalogs);
      printf("DoNotWriteHistories         %d\n",params.DoNotWriteHistories);
      printf("WriteSnapshot               %d\n",params.WriteSnapshot);
      printf("WriteTimelessSnapshot       %d\n",params.WriteTimelessSnapshot);
      printf("OutputInH100                %d\n",params.OutputInH100);
      printf("WriteFmax                   %d\n",params.WriteFmax);
      printf("WriteVmax                   %d\n",params.WriteVmax);
      printf("WriteRmax                   %d\n",params.WriteRmax);
      switch(params.AnalyticMassFunction)
	{
	case 0:
	  printf("Using Press & Schechter (1974) for the analytic mass function\n");
	  break;
	case 1:
	  printf("Using Sheth & Tormen (2001) for the analytic mass function\n");
	  break;
	case 2:
	  printf("Using Jenkins et al. (2001) for the analytic mass function\n");
	  break;
	case 3:
	  printf("Using Warren et al. (2006) for the analytic mass function\n");
	  break;
	case 4:
	  printf("Using Reed et al. (2007) for the analytic mass function\n");
	  break;
	case 5:
	  printf("Using Crocce et al. (2010) for the analytic mass function\n");
	  break;
	case 6:
	  printf("Using Tinker et al. (2008) for the analytic mass function\n");
	  break;
	case 7:
	  printf("Using Courtin et al. (2010) for the analytic mass function\n");
	  break;
	case 8:
	  printf("Using Angulo et al. (2012) for the analytic mass function\n");
	  break;
	case 9:
	  printf("Using Watson et al. (2013) for the analytic mass function\n");
	  break;
	case 10:
	  printf("Using Crocce et al. (2010) with forced universality for the analytic mass function\n");
	  break;
	default:
	  printf("Unknown value for AnalyticMassFunction, Using Watson et al. (2013)\n");
	  params.AnalyticMassFunction=9;
	  break;
	}
      printf("\n");

      printf("\n");
      printf("GENIC parameters:\n");
      printf("InputSpectrum_UnitLength_in_cm %f\n",params.InputSpectrum_UnitLength_in_cm);
      printf("FileWithInputSpectrum          %s\n",params.FileWithInputSpectrum);
      printf("WDM_PartMass_in_kev            %f\n",params.WDM_PartMass_in_kev);
      printf("\n");
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

  var_min    = pow(1.686/NSIGMA / GrowingMode(outputs.zlast,params.k_for_GM),2.0);
  rmin       = params.InterPartDist/6.;
  var_max    = MassVariance(rmin);
  Smoothing.Nsmooth = (log10(var_max)-log10(var_min))/STEP_VAR+2;

  if (Smoothing.Nsmooth<=0)
    {
      if (!ThisTask)
	printf("I am afraid that nothing is predicted to collapse in this configuration.\nI will work with no smoothing\n");
      Smoothing.Nsmooth=1;
    }

  if (!ThisTask)
    printf("Min variance: %f12.6, max variance: %f12.6, number of smoothing radii: %d\n",
	   var_min,var_max,Smoothing.Nsmooth);

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
    printf("[%s] Generating density in Fourier space\n",fdate());

#ifdef WHITENOISE

  if (Ngrids>1)
    {
      if (!ThisTask)
	printf("Sorry, this works only with a single grid\n");
      return 1;
    }

  if (read_white_noise())
    return 1;

#else

  int igrid;
  for (igrid=0; igrid<Ngrids; igrid++)
    if (GenIC(igrid))
      return 1;

#endif

  cputime.dens=MPI_Wtime()-cputime.dens;
  if (!ThisTask)
    printf("[%s] Done generating density in Fourier space, cputime = %f s\n",fdate(), cputime.dens);

  /* // LEVARE!!! */
  /* write_in_cvector(0, kdensity[0]); */
  /* (void)reverse_transform(0); */
  /* write_from_rvector(0, density[0]); */
  /* double Variance=0.0; */
  /* for (int i=0; i<MyGrids[0].total_local_size; i++) */
  /*   Variance+=density[0][i]*density[0][i]; */
  /* Variance/=(double)MyGrids[0].total_local_size; */
  /* printf("VARIANCE on the grid: %g\n",Variance); */

  /* double VVVV=VarianceOnGrid(0,1.0,0.0); */
  /* printf("CALCOLO: %g  %f  %f  %f   %10f    \n",VVVV, Variance/VVVV, Radius(Variance), params.InterPartDist, Radius(Variance)/params.InterPartDist); */

  /* WindowFunctionType=1; */
  /* if (initialize_MassVariance()) */
  /*   return 1; */

  /* double SSSS=MassVariance(params.InterPartDist/PI); */
  /* printf("SKS: %g  %g\n",SSSS, Variance/SSSS); */
  /* return 1; */

  return 0;
}


int set_grids()
{
  /* initialization of fftw quantities on grids (one for the moment) */

  int igrid;

  Ngrids=1;

  MyGrids=(grid_data*)malloc(Ngrids * sizeof(grid_data));

  MyGrids[0].GSglobal_x = params.GridSize[0];
  MyGrids[0].GSglobal_y = params.GridSize[1];
  MyGrids[0].GSglobal_z = params.GridSize[2];

  MyGrids[0].BoxSize = params.BoxSize_htrue;
  MyGrids[0].lower_k_cutoff=0.;
  MyGrids[0].upper_k_cutoff=NYQUIST * PI;

  /* allocates pointers */
  cvector_fft=(fftw_complex**)malloc(Ngrids * sizeof(fftw_complex*));
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
  seedtable=(unsigned int**)malloc(Ngrids * sizeof(unsigned int*));
 
  for (igrid=0; igrid<Ngrids; igrid++)
    if (set_one_grid(igrid))
      return 1;

  return 0;
}



#ifdef PLC
#define NSAFE 2.0

double myz;
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
      plc.center[0]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal_x;
      plc.center[1]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal_y;
      plc.center[2]=gsl_rng_uniform(random_generator)*MyGrids[0].GSglobal_z;
      plc.zvers[0]=1.0;
      plc.zvers[1]=1.0;
      plc.zvers[2]=1.0;
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

  l[0]=(double)(MyGrids[0].GSglobal_x);
  l[1]=(double)(MyGrids[0].GSglobal_y);
  l[2]=(double)(MyGrids[0].GSglobal_z);

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

  plc.Nmax = subbox.Npart / 10;

  if (!ThisTask)
    {
      printf("The Past Light Cone will be reconstruct from z=%f to z=%f\n",
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

  int i,j,k, i1,j1,k1, i2,j2,k2,NS2, NN,NN1,N1,N2,N3,ssafe;
  double tdis,size,this,BytesPerParticle,FmaxBPP,FragG_BPP,FragP_BPP,
    TotalNP,TotalNP_pertask,ratio,smallest,cc,MemPerTask;

  /* typical displacement at zlast */
  tdis = GrowingMode(outputs.zlast,params.k_for_GM) * sqrt( DisplVariance(params.InterPartDist) );

  /* mass of the largest halo expected in the box */
  params.Largest=1.e18;
  cc=1./pow(params.BoxSize_htrue,3.0);

  double aa=AnalyticMassFunction(params.Largest,outputs.zlast);
  while (aa*params.Largest<cc)
    {
      params.Largest*=0.99; 
      aa=AnalyticMassFunction(params.Largest,outputs.zlast);
    }

  /* boundary layer */
  size=SizeForMass(params.Largest);
  subbox.SafetyBorder = params.BoundaryLayerFactor * size;
  subbox.safe = (int)(subbox.SafetyBorder/params.InterPartDist)+1;

  if (!ThisTask)
    {
      printf("\n");
      printf("Determination of the boundary layer\n");
      printf("   growing mode at z=%f: %f\n",outputs.zlast, GrowingMode( outputs.zlast, params.k_for_GM));
      printf("   largest halo expected in this box at z=%f: %e Msun\n",
	     outputs.zlast, params.Largest);
      printf("   its Lagrangian size: %f Mpc\n",size);
      printf("   typical displacement: %f \n",tdis);
      printf("   the boundary layer will be %f, a factor of %f with respect to the typical displacement\n",
	     subbox.SafetyBorder, subbox.SafetyBorder/tdis);
    }

  /* finds the optimal number of sub-boxes to use for the fragmentation */
  ssafe=2.*subbox.safe;
  FmaxBPP = (double)sizeof(product_data) + 10.0*(double)sizeof(double) + 
    (double)sizeof(int) * (double)NTasks / (double)MyGrids[0].GSglobal_z;
  TotalNP = (double)MyGrids[0].GSglobal_x * (double)MyGrids[0].GSglobal_y * (double)MyGrids[0].GSglobal_z;
  TotalNP_pertask = TotalNP/(double)NTasks;

  FragP_BPP=(double)sizeof(product_data);
  FragG_BPP=3.0*(double)sizeof(int)+(double)(sizeof(group_data) +sizeof(histories_data))/10.0;
#ifdef PLC
  FragG_BPP+=(double)sizeof(plcgroup_data)/10.;
#endif

  smallest=1.e10;
  NSlices=0;

  do
    {
      ++NSlices;

      BytesPerParticle=1.e10;
      for (k=1; k<=NTasks; k++)
	for (j=1; j<=NTasks/k; j++)
	  for (i=1; i<=NTasks/k/j; i++)
	    /* the three indices must be exact divisors of the three grid lengths */
	    if (i*j*k==NTasks)
	      {
		/* number of particles in the sub-box */
		N1 = find_length(MyGrids[0].GSglobal_x,i,0);
		N2 = find_length(MyGrids[0].GSglobal_y,j,0);
		N3 = find_length(MyGrids[0].GSglobal_z,k*NSlices,0);
		if (N1<ssafe || N2<ssafe || N3<ssafe)
		  continue;
		NN = (N1 + (i==1? 0 : ssafe))
		  *  (N2 + (j==1? 0 : ssafe)) 
		  *  (N3 + (k*NSlices==1? 0 : ssafe));

		ratio = (double)NN/TotalNP_pertask;
		if (NSlices>1)
		  this=(double)sizeof(product_data) + ratio * (FragP_BPP + FragG_BPP);
		else		    
		  this=( (double)sizeof(product_data) > ratio * FragG_BPP ?
			 (double)sizeof(product_data) : ratio * FragG_BPP) +
		    ratio * FragP_BPP;
		if (this<FmaxBPP)
		  this=FmaxBPP;
		
		if (this<smallest)
		  {
		    smallest=this;
		    i2=i; j2=j; k2=k; NS2=NSlices;
		  }

		if (this < BytesPerParticle)
		  {
		    BytesPerParticle=this;
		    NN1=NN;
		    i1=i;
		    j1=j;
		    k1=k;
		  }
	      }
      if (BytesPerParticle>1000.)
	break;
    }
  while (BytesPerParticle>params.MaxMemPerParticle);


  if (BytesPerParticle>1000.)
    {
      if (!ThisTask)
	{
	  printf("ERROR: no possible division of sub-boxes found up to Nslices=%d\n", 
		 NSlices);
	  printf("lowest possible value of memory per particle is %f ",smallest);
	  printf("found on a subdivision %d-%d-%d on %d slices\n",i2,j2,k2,NS2);
	  printf("please decrease BoundaryLayerFactor or increase MaxMemPerParticle\n");
	  fflush(stdout);
	}
      return 1;
    }

  subbox.nbox_x=i1;
  subbox.nbox_y=j1;
  subbox.nbox_z_thisslice=k1;
  subbox.nbox_z_allslices=k1*NSlices;

  subbox.safe_x = (subbox.nbox_x>1 ? subbox.safe : 0);
  subbox.safe_y = (subbox.nbox_y>1 ? subbox.safe : 0);
  subbox.safe_z = (subbox.nbox_z_allslices>1 ? subbox.safe : 0);

  subbox.pbc_x = (subbox.nbox_x==1);
  subbox.pbc_y = (subbox.nbox_y==1);
  subbox.pbc_z = (subbox.nbox_z_allslices==1);

  /* this will be mybox for the first slice */
  NN=subbox.nbox_y*subbox.nbox_z_thisslice;
  subbox.mybox_x=ThisTask/NN;
  NN1=ThisTask-subbox.mybox_x*NN;
  subbox.mybox_y=NN1/subbox.nbox_z_thisslice;
  subbox.mybox_z=NN1-subbox.mybox_y*subbox.nbox_z_thisslice;

  subbox.Lgrid_x = find_length(MyGrids[0].GSglobal_x,subbox.nbox_x,subbox.mybox_x);
  subbox.Lgrid_y = find_length(MyGrids[0].GSglobal_y,subbox.nbox_y,subbox.mybox_y);
  subbox.Lgrid_z = find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.Lgwbl_x = subbox.Lgrid_x + 2*subbox.safe_x; 
  subbox.Lgwbl_y = subbox.Lgrid_y + 2*subbox.safe_y;
  subbox.Lgwbl_z = subbox.Lgrid_z + 2*subbox.safe_z;

  subbox.Npart = subbox.Lgwbl_x * subbox.Lgwbl_y * subbox.Lgwbl_z;

  subbox.start_x = find_start(MyGrids[0].GSglobal_x,subbox.nbox_x,subbox.mybox_x);
  subbox.start_y = find_start(MyGrids[0].GSglobal_y,subbox.nbox_y,subbox.mybox_y);
  subbox.start_z = find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);

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
#ifdef TIMELESS_SNAPSHOT
  if (NSlices>1 && params.WriteTimelessSnapshot)
    {
      params.WriteTimelessSnapshot=0;
      if (!ThisTask)
	printf("Sorry, but timeless snapshots cannot be written if fragmentation is done in slices\n");
    }
#endif

  /* initialization of quantities required by compute_mf */
   if (params.OutputInH100)
    mf.hfactor=params.Hubble100;
  else
    mf.hfactor=1.0;
  mf.hfactor4=pow(mf.hfactor,4.);
  mf.vol=(double)MyGrids[0].GSglobal_x*(double)MyGrids[0].GSglobal_y
    *(double)MyGrids[0].GSglobal_z*pow(params.InterPartDist,3.0);
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
      printf("\n");
      printf("FRAGMENTATION:\n");
      if (NSlices>1)
	printf("The box will be fragmented in %d slices\n",NSlices);
      else
	printf("The box will be fragmented in one go\n");
      printf("Number of sub-boxes per dimension: %d %d %d\n",subbox.nbox_x,subbox.nbox_y,subbox.nbox_z_allslices);
      printf("Boundary layer (true Mpc):         %f\n",subbox.SafetyBorder);
      printf("Boundary layer (gridpoints):       %d\n",subbox.safe);
      printf("Core 0 will work on a grid:        %d %d %d\n",subbox.Lgwbl_x,subbox.Lgwbl_y,subbox.Lgwbl_z);
      printf("Number of particles for core 0:    %d\n",subbox.Npart);
      printf("The resolved box will be:          %d %d %d\n",subbox.Lgrid_x,subbox.Lgrid_y,subbox.Lgrid_z);
      printf("Periodic boundary conditions:      %d %d %d\n",subbox.pbc_x,subbox.pbc_y,subbox.pbc_z);
      printf("Required bytes per fft particle:   %f\n",BytesPerParticle);
      printf("The overhead for fragmentation is: %f\n",subbox.overhead);
      printf("Required memory per task:          %4.0fMb - Maxmem=%dMb\n", MemPerTask*1024.,params.MaxMem);
      printf("\nThe mass function will be computed from Log M=%f to Log M=%f (%d bins)\n",
      	     mf.mmin, mf.mmax, mf.NBIN);
      printf("\n");
    }

  if (MemPerTask > params.MaxMem/1024.0)
    {
      if (!ThisTask)
      printf("ERROR: your requirements overshoot the available memory per MPI task\n");
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

  // CI del collasso ellissoidale


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
      printf("density, smoothing %d, iter=%d\n",ismooth,iter); // LEVARE
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
      printf("velocities, smoothing %d, iter=%d\n",ismooth,iter); // LEVARE
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
