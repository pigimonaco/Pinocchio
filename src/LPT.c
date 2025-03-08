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
#ifdef TWO_LPT

int compute_LPT_displacements(int compute_sources, double redshift)
{
  double time;

  /*  
      second derivative components are sorted this way
       0 -> (1,1)
       1 -> (2,2)
       2 -> (3,3)
       3 -> (1,2)
       4 -> (1,3)
       5 -> (2,3)
  */

  if (compute_sources)
    {

      /**************************/
      /* Computation of sources */
      /**************************/

      ScaleDep.order=0; /* sources are computed as for LambdaCDM */
      ScaleDep.redshift=0.0;

  /* Computes the source for 2 LPT, 3LPT_1 and 3LPT_2 term */
  /* loop on all local particles */
#ifdef _OPENMP
#pragma omp parallel
      {

#pragma omp for nowait
#endif
	for (int index=0; index<MyGrids[0].total_local_size; index++)
	  {
	      
	    /* HERE IT WILL BE NECESSARY TO SUM DIFFERENT SOURCES IF MULTIPLE GRIDS ARE USED */
	      
	    /* source term for 2LPT */
	    source_2LPT[index] = 
	      second_derivatives[0][0][index] * second_derivatives[0][1][index] +
	      second_derivatives[0][0][index] * second_derivatives[0][2][index] +
	      second_derivatives[0][1][index] * second_derivatives[0][2][index] -
	      second_derivatives[0][3][index] * second_derivatives[0][3][index] -
	      second_derivatives[0][4][index] * second_derivatives[0][4][index] -
	      second_derivatives[0][5][index] * second_derivatives[0][5][index];
		
#ifdef THREE_LPT
	    source_3LPT_1[index] = 3.0 *(second_derivatives[0][0][index]*
					 (second_derivatives[0][1][index]*second_derivatives[0][2][index]
					  -second_derivatives[0][5][index]*second_derivatives[0][5][index])
					 - second_derivatives[0][3][index]*
					 (second_derivatives[0][3][index]*second_derivatives[0][2][index]
					  -second_derivatives[0][4][index]*second_derivatives[0][5][index])
					 + second_derivatives[0][4][index]*
					 (second_derivatives[0][3][index]*second_derivatives[0][5][index]
					  -second_derivatives[0][4][index]*second_derivatives[0][1][index]));
		
	    source_3LPT_2[index] = 2.0 *     /* this factor is needed because nabla2phi is half the theoretical one */
	      (second_derivatives[0][0][index] + second_derivatives[0][1][index] + second_derivatives[0][2][index]) * 
	      source_2LPT[index];
#endif
	  }
#ifdef _OPENMP
      }
#endif

      /* forward FFT for 2LPT source */
      write_in_rvector(0, source_2LPT);

      if (!ThisTask)
	printf("[%s] source for 2LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] source for 2LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft += time;

      write_from_cvector(0, kvector_2LPT);

#ifdef THREE_LPT

      /* second derivatives of second-order potential */
      for ( int ia = 1; ia <= 3; ia++ )
	for ( int ib = ia; ib <= 3; ib++ )
	  {
	    int ider = ( ia==ib? ia : ia+ib+1 );

	    if (!ThisTask)
	      printf("[%s] Computing 2nd derivative of 2LPT source: %d\n",fdate(),ia);

	    write_in_cvector(0, kvector_2LPT);

	    if (compute_derivative(0,ia,ib))
	      return 1;

#ifdef _OPENMP
#pragma omp parallel
	    {
#pragma omp for nowait
#endif
	      for (int index=0; index<MyGrids[0].total_local_size; index++)
		/* the first 2 factor is needed because nabla2phi is half the theoretical one */
		source_3LPT_2[index] -= 2.0 * (ider<=3? 1.0 : 2.0) *
		  rvector_fft[0][index] * second_derivatives[0][ider-1][index];
#ifdef _OPENMP
	    }
#endif
	  }

      /* forward FFT for 3LPT_1 source */
      write_in_rvector(0, source_3LPT_1);

      if (!ThisTask)
	printf("[%s] source for 3LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] sources for 3LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft+=time;

      write_from_cvector(0, kvector_3LPT_1);


      /* forward FFT for 3LPT_2 source */
      write_in_rvector(0, source_3LPT_2);

      if (!ThisTask)
	printf("[%s] source for 3LPT: starting fft\n",fdate());

      time = forward_transform(0);

      if (!ThisTask)
	printf("[%s] source for 3LPT: done fft, cputime = %f\n",fdate(),time);

      cputime.fft+=time;

      write_from_cvector(0, kvector_3LPT_2);

#endif
    }

  /* displacements for the three terms */
  if (!ThisTask)
    printf("[%s] Computing displacements\n",fdate());

  ScaleDep.order=2; /* here we need the 2LPT  */
  ScaleDep.redshift=redshift;

  compute_first_derivatives(0., 0, 2, kvector_2LPT);

/*   for ( int ia = 1; ia <= 3; ia++ ) */
/*     { */
/*       if (!ThisTask) */
/* 	printf("[%s] Computing 1st derivative of 2LPT source: %d\n",fdate(),ia); */

/*       write_in_cvector(0, kvector_2LPT); */

/*       if (compute_derivative(0,ia,0)) */
/* 	return 1; */

/*       write_from_rvector(0, first_derivatives[0][ia-1]); */
/*     } */

/* #ifdef _OPENMP */
/* #pragma omp parallel */
/*   { */
/*     /\* assigns displacement to particles *\/ */
    
/* #pragma omp for nowait */
/* #endif */
/*     for (int index=0; index<MyGrids[0].total_local_size; index++) */
/*       for ( int ia = 0; ia < 3; ia++ ) */
/* 	products[index].Vel_2LPT[ia] = first_derivatives[0][ia][index]; */

/* #ifdef _OPENMP */
/*   } */
/* #endif */

#ifdef THREE_LPT

  if (!ThisTask)
    printf("[%s] Computing 3LPT_1 displacements\n",fdate());

  ScaleDep.order=3; /* here we need the 3LPT_1  */

  compute_first_derivatives(0., 0, 3, kvector_3LPT_1);

  if (!ThisTask)
    printf("[%s] Computing 3LPT_2 displacements\n",fdate());

  ScaleDep.order=4; /* here we need the 3LPT_2  */

  compute_first_derivatives(0., 0, 4, kvector_3LPT_2);

#endif


  /* bye! */
  return 0;
}

#endif
