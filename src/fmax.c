/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan, 
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025
 
 github: https://github.com/pigimonaco/Pinocchio
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

int compute_second_derivatives(double, int);
int Fmax_PDF(void);


/* Computation of collapse times and displacements */
int compute_fmax(void)
{
  int ismooth,ThisGrid;
  double cputmp,cpusm;

  cputime.fmax=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] First part: computation of collapse times\n",fdate());

  ThisGrid=0;

  /*
    TO EXTEND TO MULTIPLE GRIDS here we should

    1) initialize fft for the large-scale grid (the same as the small-scale one?)
    2) compute first derivatives of the large grid
    3) compute second derivatives of the large grid
    4) distribute to all processors and store the few used grid points
    5) re-initialize fft for the small-scale grid
  */

  ScaleDep.order=0; /* collapse times must be computed as for LambdaCDM */
  ScaleDep.redshift=0.0;


  /****************************
   * CYCLE ON SMOOTHING RADII *
   ****************************/

  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {
      if (!ThisTask)
	printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n", 
	       fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
	       sqrt(Smoothing.Variance[ismooth]) );

      cpusm=MPI_Wtime();

      /* 
	 Part 1:
	 Compute second derivatives of the potential
      */

      if (!ThisTask)
	printf("[%s] Computing second derivatives\n",fdate());

      cputmp=MPI_Wtime();

      if (compute_second_derivatives(Smoothing.Radius[ismooth] ,ThisGrid))
	return 1;

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done second derivatives, cpu time = %f s\n",fdate(),cputmp);
      cputime.deriv += cputmp;

      /* 
	 Part 2:
	 Compute collapse times
      */

      if (!ThisTask)
	printf("[%s] Computing collapse times\n",fdate());

      cputmp=MPI_Wtime();

#ifdef TABULATED_CT
      /* initialize spline for interpolating collapse times */
      if (initialize_collapse_times(ismooth,0))
	return 1;

      if (!ThisTask)
	{
	  if (strcmp(params.CTtableFile,"none"))
	    printf("[%s] Collapse times read from file %s\n",fdate(),params.CTtableFile);
	  else
	    printf("[%s] Collapse times computed for interpolation, cpu time =%f s\n",fdate(),cputmp);
	}
#endif // TABULATED_CT

      if (
#if defined(GPU_OMP) || defined(GPU_OMP_FULL)
	  compute_collapse_times_gpu(ismooth)
#else
	  compute_collapse_times(ismooth)
#endif // GPU_OMP || GPU_OMP_FULL

	  )
	return 1;

#ifdef TABULATED_CT
      /* this is needed only for debug options, to be removed in the official code */
      if (reset_collapse_times(ismooth))
	return 1;
#endif // TABULATED_CT

      // cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing collapse times, cpu time = %f s\n",fdate(),cputmp);
      // cputime.coll+=cputmp;

      /*
	End of cycle on smoothing radii
      */

      cpusm = MPI_Wtime()-cpusm;

      if (!ThisTask)
	printf("[%s] Completed, R=%6.3f, expected sigma: %7.4f, computed sigma: %7.4f, cpu time = %f s\n",fdate(),
	       Smoothing.Radius[ismooth], sqrt(Smoothing.Variance[ismooth]), 
	       sqrt(Smoothing.TrueVariance[ismooth]), 
	       cpusm );

      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }



#if defined(CUSTOM_INTERPOLATION) || defined(GPU_OMP) 
  custom_cubic_spline_free(host_spline);

#endif // defined(CUSTOM_INTERPOLATION) || defined(GPU_OMP)

  
#if defined(GPU_OMP)
  /*---------------- Free GPU/CPU memory ----------------------*/

  /*----- Free GPU spline ----*/
  omp_target_free(internal.device.gpu_main_memory, devID);
  free(host_products.Rmax);
  free(host_products.Fmax);
  
#elif defined(GPU_OMP_FULL) // GPU_OMP 
  
  /*-------------------- Free GPU products and second derivatives -------------------------- */
  #pragma omp target exit data map(delete: gpu_products.Rmax[0:MyGrids[0].total_local_size], \
                                           gpu_products.Fmax[0:MyGrids[0].total_local_size], \
                                           gpu_products) device(devID)

   for (int igrid = 0; igrid < Ngrids; igrid++) 
   {
    for (int i = 0; i < 6; i++) 
    {
      #pragma omp target exit data map(delete: second_derivatives[igrid][i][0: MyGrids[igrid].total_local_size]) device(devID)
    }
   }
 
  #pragma omp target exit data map(delete: second_derivatives[0:Ngrids][0:6]) device(devID)
  
  /*-------------------- Free GPU Splines and Smoothing -------------------------- */
  #pragma omp target exit data map(delete: Smoothing.Radius[0:Smoothing.Nsmooth],  \
                                           Smoothing.Variance[0:Smoothing.Nsmooth],\
                                           Smoothing.TrueVariance[0:Smoothing.Nsmooth]) device(devID)

   /* Free gpu_spline */
  #pragma omp target exit data map(delete: gpu_spline[0:1],		               \
                                           gpu_spline->d2y_data[0:NBINS],    \
                                           gpu_spline->coeff_a[0:NBINS - 1], \
                                           gpu_spline->coeff_b[0:NBINS - 1], \
                                           gpu_spline->coeff_c[0:NBINS - 1], \
                                           gpu_spline->coeff_d[0:NBINS - 1]) device(devID)

   
  #if defined(TABULATED_CT)
  #define CT_NBINS_XY (50) 
  #define CT_NBINS_D (100)
  
  /* Free GPU memory */
  for (int i = 0; i < CT_NBINS_XY; ++i) {
    for (int j = 0; j < CT_NBINS_XY; ++j) {
        int index = i * CT_NBINS_XY + j;

        #pragma omp target exit data map(delete: CT_Spline[index]->d2y_data[0:CT_NBINS_D],    \
                                                 CT_Spline[index]->coeff_a[0:CT_NBINS_D - 1], \
                                                 CT_Spline[index]->coeff_b[0:CT_NBINS_D - 1], \
                                                 CT_Spline[index]->coeff_c[0:CT_NBINS_D - 1], \
	                                               CT_Spline[index]->coeff_d[0:CT_NBINS_D - 1], \
	                                               CT_Spline[index][0:1]) device(devID)
      }
  }

  #pragma omp target exit data map(delete: CT_Spline[0:CT_NBINS_XY * CT_NBINS_XY]) device(devID)

  // Ensure you free memory properly when no longer needed on the CPU
   for (int i = 0; i < CT_NBINS_XY; ++i)
   {
      for (int j = 0; j < CT_NBINS_XY; ++j)
      {
        const int index = i * CT_NBINS_XY + j;

	if (CT_Spline[index]->x != NULL)
	  free(CT_Spline[index]->x);
	if (CT_Spline[index]->y != NULL)
	  free(CT_Spline[index]->y);

	if (CT_Spline[index]->d2y_data != NULL)
	  free(CT_Spline[index]->d2y_data);

	if (CT_Spline[index]->coeff_a != NULL)
	  free(CT_Spline[index]->coeff_a);

	if (CT_Spline[index]->coeff_b != NULL)
	  free(CT_Spline[index]->coeff_b);

	if (CT_Spline[index]->coeff_c != NULL)
	  free(CT_Spline[index]->coeff_c);

	if (CT_Spline[index]->coeff_d != NULL)
	  free(CT_Spline[index]->coeff_d);
      }
   }

   free(CT_Spline);
#endif // CT_TABLE
#endif // FULL_GPU_OMP

  /********************************
   * COMPUTATION OF DISPLACEMENTS *
   ********************************/

  /* computes the displacements for all particles at zero smoothing, 
     assuming second derivatives are already in place 
     displacements are computed to the redshift of the first (or only) segment
  */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("\n[%s] Computing displacements  for redshift %f\n",fdate(),ScaleDep.z[0]);

  if (compute_displacements(1,0,ScaleDep.z[0]))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done computing displacements, cpu time = %f s\n",fdate(),cputmp);

  /***************/
  /* END OF FMAX */
  /***************/

  if (finalize_fft())
    return 1;

  if (Fmax_PDF())
    return 1;

  cputime.fmax = MPI_Wtime() - cputime.fmax;
  if (!ThisTask)
    printf("[%s] Finishing fmax, total fmax cpu time = %14.6f\n"
	   "\t\t IO       : %14.6f (%14.6f total time without I/O)\n"
	   "\t\t FFT      : %14.6f\n"
	   "\t\t COLLAPSE : %14.6f\n",
	   fdate(), cputime.fmax, cputime.io, cputime.fmax-cputime.io, cputime.fft, cputime.coll);

  return 0;
}


int compute_first_derivatives(double R, int ThisGrid, int order, double* vector)
{
  /* computes second derivatives of the potential */

  double timetmp = 0;
  
  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for (int ia=1; ia<=3; ia++)
    {

      if (!ThisTask)
	printf("[%s] Computing 1st derivative: %d\n",fdate(),ia);

      double tmp = MPI_Wtime();
      write_in_cvector(ThisGrid, vector);
      timetmp += MPI_Wtime() - tmp;

      if (compute_derivative(ThisGrid,ia,0))
	return 1;

      tmp = MPI_Wtime();
      write_from_rvector_to_products(ThisGrid, ia-1, order);
      timetmp += MPI_Wtime() - tmp;
    }

  cputime.mem_transf += timetmp;
  return 0;
}


int compute_second_derivatives(double R, int ThisGrid)
{

  /* computes second derivatives of the potential */

  double timetmp = 0;
  
  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for ( int ia = 1; ia <= 3; ia++ )
    for ( int ib = ia; ib <= 3; ib++ )
      {

	int ider=( ia == ib ? ia : ia+ib+1 );

	if (!ThisTask)
	  printf("[%s] Computing 2nd derivative: %d\n",fdate(),ider);

	double tmp = MPI_Wtime();
	write_in_cvector(ThisGrid, kdensity[ThisGrid]);
	timetmp += MPI_Wtime() - tmp;

	if (compute_derivative(ThisGrid,ia,ib))
	  return 1;

	tmp = MPI_Wtime();
	write_from_rvector(ThisGrid, second_derivatives[ThisGrid][ider-1]);
	timetmp += MPI_Wtime() - tmp;
      }

  cputime.mem_transf += timetmp;
  return 0;
}


char *fdate()
{
  /* returns a 24-char string with the full date and time
     (the year is moved at the middle of the string) */

  time_t current_time;
  char *string;
  int n;

  /* format from ctime:
    0123456789
              0123456789
	                0123
    Www Mmm dd hh:mm:ss yyyy
  */

  current_time=time(NULL);
  string=ctime(&current_time);

  for (n=0; n<10; n++)
    *(date_string+n)=*(string+n);
  for (n=10; n<15; n++)
    *(date_string+n)=*(string+n+9);
  for (n=10; n<19; n++)
    *(date_string+n+5)=*(string+n);
  *(date_string+24)='\0';

  return date_string;
}


int compute_displacements(int compute_sources, int recompute_sd, double redshift)
{
  /* This routine computes the displacement fields at a specified redshift
     if recompute_sd is true, it recomputes the second derivatives */

  double cputmp;

#ifdef TWO_LPT

  if (recompute_sd)
    {
      /* second derivatives are needed for LPT displacements, if not already in place they must be recomputed */
      cputmp=MPI_Wtime();
      if (!ThisTask)
	printf("\n[%s] Computing second derivatives\n",fdate());

      ScaleDep.order=0;  /* sources are computed as for LambdaCDM at z=0 */
      ScaleDep.redshift=0.0;

      if (compute_second_derivatives(0.0, 0))
	return 1;

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done second derivatives, cpu time = %f s\n",fdate(),cputmp);
      cputime.deriv+=cputmp;
    }

  /* computes the 2LPT and 3LPT displacement fields and stores them in the products */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("\n[%s] Computing LPT displacements\n",fdate());

  if (compute_LPT_displacements(compute_sources, redshift))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done LPT displacements, cpu time = %f s\n",fdate(),cputmp);
  cputime.lpt+=cputmp;

#endif

  /* computes Zeldovich displacements */
  cputmp=MPI_Wtime();
  if (!ThisTask)
    printf("[%s] Computing first derivatives\n",fdate());

  /* for scale-dependent growth we compute the displacements for the first redshift segment 
     in case there is only one segment, this implies to compute them at the final redshift */
  ScaleDep.redshift=redshift;
  ScaleDep.order=1;   /* here we need the first-order growth */

  if (compute_first_derivatives(0.0, 0, 1, kdensity[0]))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done first derivatives, cpu time = %f s\n",fdate(),cputmp);
  cputime.deriv+=cputmp;

  /* /\* Store Zeldovich displacements in the products *\/ */
  /* cputmp=MPI_Wtime(); */
  /* if (!ThisTask) */
  /*   printf("[%s] Storing velocities\n",fdate()); */

  /* if (store_velocities()) */
  /*   return 1; */

  /* cputmp=MPI_Wtime()-cputmp; */
  /* if (!ThisTask) */
  /*   printf("[%s] Done storing velocities, cpu time = %f s\n",fdate(),cputmp); */
  /* cputime.vel+=cputmp; */

  return 0;
}


#include <sys/stat.h>

int dump_products()
{
  /* dumps fmax products on files to reuse them */

  struct stat dr;
  FILE *file;
  char fname[LBLENGTH];

  /* Task 0 checks that the DumpProducts directory exists, or creates it 
     and writes a summary file to check that one does not start from a different run */
  if (!ThisTask)
    {
      if(stat(params.DumpDir,&dr))
	  {
	    printf("Creating directory %s\n",params.DumpDir);
	    if (mkdir(params.DumpDir,0755))
	      {
		printf("ERROR IN CREATING DIRECTORY %s (task 0)\n",params.DumpDir);
		return 1;
	      }
	  }

      sprintf(fname,"%ssummary",params.DumpDir);
      file=fopen(fname,"w");
      fprintf(file,"%d   # NTasks\n",NTasks);
      fprintf(file,"%d   # random seed\n",params.RandomSeed);
      fprintf(file,"%d   # grid size\n",params.GridSize[0]);
      fprintf(file,"%d   # length of product_data\n",(int)sizeof(product_data));
      fclose(file);

      /* dumps the true variances */
      sprintf(fname,"%sTrueVariance",params.DumpDir);
      file=fopen(fname,"wb");
      fwrite(Smoothing.TrueVariance, sizeof(double), Smoothing.Nsmooth, file);
      fclose(file);
    }

  /* all tasks must wait for Task 0 to open the directory (if needed) */
  MPI_Barrier(MPI_COMM_WORLD);

  /* each task writes a separate dump file */
  sprintf(fname,"%sTask.%d",params.DumpDir,ThisTask);
  file=fopen(fname,"wb");
  if (file==0x0)
    {
      printf("ERROR on Task %d: could not open file %s\n",ThisTask, fname);
      return 1;
    }

  // NOTE: with LEAPFROG one needs to swap also kdensity and kvectors
  fwrite(products, sizeof(product_data), MyGrids[0].total_local_size, file);
  fclose(file);

  return 0;
}


int read_dumps()
{
  /* reads in fmax products from files */
  FILE *file;
  char fname[LBLENGTH],buf[SBLENGTH];
  int myNTasks, myRandomSeed, myGridSize, myPDlength;

  /* Task 0 checks that the DumpProducts/summary exists and is compatible with the present run */
  if (!ThisTask)
    {
      sprintf(fname,"%ssummary",params.DumpDir);
      file=fopen(fname,"r");
      (void)fgets(buf, SBLENGTH, file);
      sscanf(buf,"%d",&myNTasks);
      (void)fgets(buf, SBLENGTH, file);
      sscanf(buf,"%d",&myRandomSeed);
      (void)fgets(buf, SBLENGTH, file);
      sscanf(buf,"%d",&myGridSize);
      (void)fgets(buf, SBLENGTH, file);
      sscanf(buf,"%d",&myPDlength);
      fclose(file);
      
      int error=0;
      if (NTasks != myNTasks)
	{
	  printf("ERROR: the number of tasks in %s does not match - %d vs %d\n",
		 fname, myNTasks, NTasks);
	  error++;
	}
      if (params.RandomSeed != myRandomSeed)
	{
	  printf("ERROR: the random seed in %s does not match - %d vs %d\n",
		 fname, myRandomSeed, params.RandomSeed);
	  error++;
	}
      if (params.GridSize[0] != myGridSize)
	{
	  printf("ERROR: the grid size in %s does not match - %d vs %d\n",
		 fname, myGridSize, params.GridSize[0]);
	  error++;
	}
      if (myPDlength != (int)sizeof(product_data))
	{
	  printf("ERROR: the length of product_data in %s does not match - %d vs %d\n",
		 fname, myPDlength, (int)sizeof(product_data));
	  error++;
	}
      if (error)
	return 1;

      /* reads and broadcasts true variances */
      sprintf(fname,"%sTrueVariance",params.DumpDir);
      file=fopen(fname,"rb");
      if (file==0x0)
	{
	  printf("ERROR on Task 0: could not open file %s\n",fname);
	  return 1;
	}
      fread(Smoothing.TrueVariance, sizeof(double), Smoothing.Nsmooth, file);
      fclose(file);
    }
  MPI_Bcast(Smoothing.TrueVariance, Smoothing.Nsmooth, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  /* each task reads a separate dump file */
  sprintf(fname,"%sTask.%d",params.DumpDir,ThisTask);
  file=fopen(fname,"rb");
  if (file==0x0)
    {
      printf("ERROR on Task %d: could not open file %s\n",ThisTask, fname);
      return 1;
    }

  // NOTE: with LEAPFROG one needs to read in also kdensity and kvectors
  fread(products, sizeof(product_data), MyGrids[0].total_local_size, file);
  fclose(file);

  return 0;
}


int Fmax_PDF(void)
{

  unsigned long long my_counter[NBINS], counter[NBINS], coll;

  for (int i=0; i<NBINS; i++)
    my_counter[i]=0;

  for (int i=0; i<MyGrids[0].total_local_size; i++)
    {
      int xF = (int)(products[i].Fmax*10.);
      if (xF<0)
	xF=0;
      if (xF>=NBINS)
	xF=NBINS-1;
      my_counter[xF]++;
    }

  MPI_Reduce(my_counter, counter, NBINS, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      coll=0;
      for (int i=10; i<NBINS; i++)
	coll+=counter[i];
      printf("[%s] Number of collapsed particles to z=0: %Lu\n",fdate(),coll);

      char filename[LBLENGTH];
      sprintf(filename,"pinocchio.%s.FmaxPDF.out",params.RunFlag);
      FILE *file=fopen(filename,"w");
      fprintf(file, "# Fmax PDF over %Lu particles\n",MyGrids[0].Ntotal);
      fprintf(file, "# 1-2: F interval\n");
      fprintf(file, "# 3: number of particles in that interval\n");
      fprintf(file, "#\n");
      for (int i=0; i<NBINS; i++)
	fprintf(file, " %6.1f   %6.1f  %Lu\n",(double)i/10., 
		(i==NBINS-1? 999.0 : (double)(i+1)/10.), counter[i]);
      fclose(file);
    }

  return 0;
}
