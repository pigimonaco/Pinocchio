/*****************************************************************
 *                        PINOCCHI0  V4.0                        *
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

int compute_second_derivatives(double, int);
int compute_first_derivatives(double, int);


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

    1) initialize fftw for the large-scale grid (the same as the small-scale one?)
    2) compute first derivatives of the large grid
    3) compute second derivatives of the large grid
    4) distribute to all processors and store the few used grid points
    5) re-initialize fftw for the small-scale grid
  */


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

      if (compute_collapse_times(ismooth))
	return 1;

      cputmp = MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing collapse times, cpu time = %f s\n",fdate(),cputmp);
      cputime.coll += cputmp;

#ifdef TWO_LPT
#ifndef SMOOTH_VELOCITIES
      if (ismooth==Smoothing.Nsmooth-1)
#endif
	{
	  /* computes sources and displacements for 2LPT (and 3LPT in case...) */
	  if (!ThisTask)
	    printf("[%s] Computing LPT displacements\n",fdate());

	  cputmp=MPI_Wtime();

	  if (compute_LPT_displacements(ismooth))
	    return 1;

	  cputmp=MPI_Wtime()-cputmp;
	  if (!ThisTask)
	    printf("[%s] Done LPT displacements, cpu time = %f s\n",fdate(),cputmp);
	  cputime.lpt+=cputmp;
	}
#endif

      /* 
	 Part 3:
	 Compute first derivatives of the potential
      */

#ifndef SMOOTH_VELOCITIES
      if (ismooth==Smoothing.Nsmooth-1)
#endif
	{
	  if (!ThisTask)
	    printf("[%s] Computing first derivatives\n",fdate());

	  cputmp=MPI_Wtime();

	  if (compute_first_derivatives(Smoothing.Radius[ismooth] ,ThisGrid))
	    return 1;

	  cputmp=MPI_Wtime()-cputmp;
	  if (!ThisTask)
	    printf("[%s] Done first derivatives, cpu time = %f s\n",fdate(),cputmp);
	  cputime.deriv += cputmp;
	  
      /* 
	 Part 4:
	 Store velocities of collapsed particles
      */

	  if (!ThisTask)
	    printf("[%s] Computing velocities\n",fdate());

	  cputmp=MPI_Wtime();

	  if (compute_velocities(ismooth))
	    return 1;

	  cputmp=MPI_Wtime()-cputmp;
	  if (!ThisTask)
	    printf("[%s] Done computing velocities, cpu time = %f s\n",fdate(),cputmp);
	  cputime.vel+=cputmp;
	}

      /*
	End of cycle on smoothing radii
      */

      cpusm = MPI_Wtime()-cpusm;

      if (!ThisTask)
	printf("Smoothing radius completed: %d, R=%6.3f, expected rms: %7.4f, computed rms: %7.4f, cpu time = %f s\n",
	       ismooth, Smoothing.Radius[ismooth], sqrt(Smoothing.Variance[ismooth]), 
	       sqrt(Smoothing.TrueVariance[ismooth]), cpusm );

      fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
    }

  /****************************************************/
  /* FINAL PART                                       */
  /* Output density, Fmax and Vel* fields if required */
  /****************************************************/

  /* Writes the results in files in the Data/ directory */
  cputmp=MPI_Wtime();

#ifndef TEST_ONLY
  if (write_fields())
    return 1;

  cputime.io += MPI_Wtime()-cputmp;
#endif

  if (finalize_fft())
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


int compute_first_derivatives(double R, int ThisGrid)
{
  /* computes second derivatives of the potential */

  double timetmp = 0;
  
  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for (int ia = 1; ia <= 3; ia++)
    {

      if (!ThisTask)
	printf("[%s] Computing 1st derivative: %d\n",fdate(),ia);

      double tmp = MPI_Wtime();
      write_in_cvector(ThisGrid, kdensity[ThisGrid]);
      timetmp += MPI_Wtime() - tmp;
      
      if (compute_derivative(ThisGrid,ia,0))
	return 1;

      tmp = MPI_Wtime();
      write_from_rvector(ThisGrid, first_derivatives[ThisGrid][ia-1]);
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


int compute_displacements(void)
{
#ifdef TWO_LPT
  if (compute_second_derivatives(0.0, 0))
    return 1;
  if (compute_LPT_displacements(0))
    return 1;

  for ( int i = 0; i < 3; i++ )
    for (int j = 0; j < MyGrids[0].total_local_size; j++ )
      second_derivatives[0][i+3][j] = second_derivatives[0][i][j];

  for ( int i = 0; i < 3; i++ )
    VEL2_for_displ[i] = second_derivatives[0][i+3];

#endif

  if (compute_first_derivatives(0.0, 0))
    return 1;

  for ( int i = 0; i < 3; i++ )
    VEL_for_displ[i] = first_derivatives[0][i];

  return 0;
}
