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

#ifdef SCALE_DEPENDENT
  ScaleDep.order=0; /* collapse times must be computed as for LambdaCDM */
#endif


  /****************************
   * CYCLE ON SMOOTHING RADII *
   ****************************/
  for (ismooth=0; ismooth<Smoothing.Nsmooth; ismooth++)
    {

      cpusm=MPI_Wtime();

      if (!ThisTask)
	printf("\n[%s] Starting smoothing radius %d of %d (R=%9.5f, sigma=%9.5f)\n", 
	       fdate(), ismooth+1, Smoothing.Nsmooth, Smoothing.Radius[ismooth],
	       sqrt(Smoothing.Variance[ismooth]) );

      /* 
	 Part 1:
	 Compute second derivatives of the potential
      */

      cputmp=MPI_Wtime();

      if (!ThisTask)
	printf("[%s] Computing second derivatives\n",fdate());

      if (compute_second_derivatives(Smoothing.Radius[ismooth] ,ThisGrid))
	return 1;

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done second derivatives, cpu time = %f s\n",fdate(),cputmp);

      /* 
	 Part 2:
	 Compute collapse times
      */

      cputmp=MPI_Wtime();

      if (!ThisTask)
	printf("[%s] Computing collapse times\n",fdate());

#ifdef TABULATED_CT
      /* initialize spline for interpolating collapse times */
      if (initialize_collapse_times(ismooth))
	return 1;

      cputmp=MPI_Wtime()-cputmp;

      if (!ThisTask)
	printf("[%s] Collapse times computed for interpolation, cpu time =%f s\n",fdate(),cputmp);

      cputime.coll+=cputmp;
      cputmp=MPI_Wtime();
#endif

      if (compute_collapse_times(ismooth))
	return 1;

#ifdef TABULATED_CT
      /* this is needed only for debug options, to be removed in the official code */
      if (reset_collapse_times(ismooth))
	return 1;
#endif

      cputmp=MPI_Wtime()-cputmp;
      if (!ThisTask)
	printf("[%s] Done computing collapse times, cpu time = %f s\n",fdate(),cputmp);
      cputime.coll+=cputmp;

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

  /********************************
   * COMPUTATION OF DISPLACEMENTS *
   ********************************/

#ifdef TWO_LPT
  ismooth=Smoothing.Nsmooth-1;

  /* computes sources and displacements for 2LPT and 3LPT */
  cputmp=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] Computing LPT displacements\n",fdate());

  if (compute_LPT_displacements(ismooth))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done LPT displacements, cpu time = %f s\n",fdate(),cputmp);
  cputime.lpt+=cputmp;

#endif

#ifndef RECOMPUTE_DISPLACEMENTS
  /* displacements are computed during fragmentation in this case */

  cputmp=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] Computing first derivatives\n",fdate());

  if (compute_first_derivatives(Smoothing.Radius[ismooth] ,ThisGrid))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done first derivatives, cpu time = %f s\n",fdate(),cputmp);

  /* Store velocities of collapsed particles */
  cputmp=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] Storing velocities\n",fdate());

  if (compute_velocities(ismooth))
    return 1;

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done computing velocities, cpu time = %f s\n",fdate(),cputmp);
  cputime.vel+=cputmp;
#endif

  /****************************************************/
  /* FINAL PART                                       */
  /* Output density, Fmax and Vel* fields if required */
  /****************************************************/

  /* Writes the results in files in the Data/ directory if required */
  cputmp=MPI_Wtime();

  if (write_fields())
    return 1;

  cputime.io+=MPI_Wtime()-cputmp;

  if (finalize_fftw())
    return 1;

  cputime.fmax = MPI_Wtime() - cputime.fmax;
  if (!ThisTask)
    printf("[%s] Finishing fmax, total cpu time = %14.6f\n", fdate(), cputime.fmax);

  return 0;
}


int compute_first_derivatives(double R, int ThisGrid)
{
  /* computes second derivatives of the potential */

  int ia;

  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for (ia=1; ia<=3; ia++)
    {

      if (!ThisTask)
	printf("[%s] Computing 1st derivative: %d\n",fdate(),ia);

      write_in_cvector(ThisGrid, kdensity[ThisGrid]);

      if (compute_derivative(ThisGrid,ia,0))
	return 1;

      write_from_rvector(ThisGrid, first_derivatives[ThisGrid][ia-1]);

    }

  return 0;
}


int compute_second_derivatives(double R, int ThisGrid)
{

  /* computes second derivatives of the potential */

  int ia,ib,ider;

  /* smoothing radius in grid units */
  Rsmooth = R / MyGrids[ThisGrid].CellSize;

  for (ia=1; ia<=3; ia++)
    for (ib=ia; ib<=3; ib++)
      {

	ider=( ia==ib? ia : ia+ib+1 );

	if (!ThisTask)
	  printf("[%s] Computing 2nd derivative: %d\n",fdate(),ider);

	write_in_cvector(ThisGrid, kdensity[ThisGrid]);

	if (compute_derivative(ThisGrid,ia,ib))
	  return 1;

	write_from_rvector(ThisGrid, second_derivatives[ThisGrid][ider-1]);

      }

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
  int i;

#ifdef TWO_LPT
  int j;

  if (compute_second_derivatives(0.0, 0))
    return 1;
  if (compute_LPT_displacements(0))
    return 1;

  for (i=0; i<3; i++)
    for (j=0; j<MyGrids[0].total_local_size; j++)
      second_derivatives[0][i+3][j]=second_derivatives[0][i][j];

  for (i=0; i<3; i++)
    VEL2_for_displ[i]=second_derivatives[0][i+3];

#endif

  if (compute_first_derivatives(0.0, 0))
    return 1;

  for (i=0; i<3; i++)
    VEL_for_displ[i]=first_derivatives[0][i];

  return 0;
}
