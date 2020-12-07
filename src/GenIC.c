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

/*

  The code contained in this file has been adapted from the 2LPTic code
  by Roman Scoccimarro, downloadable from:

  http://cosmo.nyu.edu/roman/2LPT/

  The 2LPTic code is a 2LPT extension of the N-GenIC code by V. Springel,
  downloadable from:

  http://www.gadgetcode.org/right.html#ICcode

  and distributed under GNU/GPL license.  The part of the 2LPTic code
  we use here is entirely contained in N-GenIC.

*/


#include "pinocchio.h"

//#define DEBUG

int GenIC(int ThisGrid)
{
  int i, j, k, ii, jj;
  int Nmesh, Nsample, Local_nx, Local_x_start, VectorLength;
  double Box, fac;
  double kvec[3], kmag, kmag2, p_of_k;
  double delta, phase, ampl;

#ifdef DEBUG
  FILE *numbers;
  char filename[300];
#endif


  Local_nx      = MyGrids[ThisGrid].GSlocal_k_y;
  Local_x_start = MyGrids[ThisGrid].GSstart_k_y;
  VectorLength  = MyGrids[ThisGrid].total_local_size_fft;
  Nmesh         = MyGrids[ThisGrid].GSglobal_x;
  Nsample       = MyGrids[ThisGrid].GSglobal_x;
  Box           = MyGrids[ThisGrid].BoxSize;

  fac = pow(1./Box,1.5);

  gsl_rng_set(random_generator, params.RandomSeed);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	seedtable[ThisGrid][i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[ThisGrid][j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[ThisGrid][(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[ThisGrid][(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[ThisGrid][i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[ThisGrid][j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	seedtable[ThisGrid][(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	seedtable[ThisGrid][(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }


  /* the vector is initialized to zero */
  for (i = 0; i < VectorLength; i++)
    kdensity[ThisGrid][i]=0.;

#ifdef DEBUG
  sprintf(filename,"random_numbers_pin.%d",ThisTask);
  numbers=fopen(filename,"w");
#endif

  /* this is the main loop on grid modes */
  for(i = 0; i < Nmesh; i++)
    {
      ii = Nmesh - i;
      if(ii == Nmesh)
	ii = 0;
      if((i >= Local_x_start && i < (Local_x_start + Local_nx)) ||
	 (ii >= Local_x_start && ii < (Local_x_start + Local_nx)))
	{
	  for(j = 0; j < Nmesh; j++)
	    {
	      if(i < Nmesh / 2)
		kvec[0] = i * 2 * PI / Box;
	      else
		kvec[0] = -(Nmesh - i) * 2 * PI / Box;

	      if(j < Nmesh / 2)
		kvec[1] = j * 2 * PI / Box;
	      else
		kvec[1] = -(Nmesh - j) * 2 * PI / Box;

	      gsl_rng_set(random_generator, seedtable[ThisGrid][i * Nmesh + j]);

	      for(k = 0; k < Nmesh / 2; k++)
		{
		  if(k < Nmesh / 2)
		    kvec[2] = k * 2 * PI / Box;
		  else
		    kvec[2] = -(Nmesh - k) * 2 * PI / Box;

		  phase = gsl_rng_uniform(random_generator) * 2 * PI;
		  do
		    ampl = gsl_rng_uniform(random_generator);
		  while(ampl == 0);

		  if(i == Nmesh / 2 || j == Nmesh / 2 || k == Nmesh / 2)
		    continue;
		  if(i == 0 && j == 0 && k == 0)
		    continue;

		  kmag2 = kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2];
		  kmag = sqrt(kmag2);

		  if(kmag * Box / (2 * PI) > NYQUIST * Nsample / 2)
		    continue;

		  p_of_k = PowerSpectrum(kmag);
#ifndef NO_RANDOM_MODULES
		  p_of_k *= -log(ampl);
#endif
		  delta = fac * sqrt(p_of_k);

		  if(k > 0)
		    {
		      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
			{

			  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ] =  delta * cos(phase);
			  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1] =  delta * sin(phase);

#ifdef DEBUG
			  fprintf(numbers,"%6d %6d %6d %6d %12f %12f %20g %20g %12f %12f\n",i,ii,j,k,phase/2/PI,ampl,
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ],
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1], kmag, PowerSpectrum(kmag));
#endif
			}

		      /*  ORIGINAL
			for(axes = 0; axes < 3; axes++)
			{
			cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
			-kvec[axes] / kmag2 * delta * sin(phase);
			cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
			kvec[axes] / kmag2 * delta * cos(phase);
			}
		      */
		    }
		  else	/* k=0 plane needs special treatment */
		    {
		      if(i == 0)
			{
			  if(j >= Nmesh / 2)
			    continue;
			  else
			    {
			      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				{
				  jj = Nmesh - j;	/* note: j!=0 surely holds at this point */
				  
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ] = delta * cos(phase);
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1] = delta * sin(phase);

				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)  ] = delta * cos(phase);
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)+1] = -delta * sin(phase);

#ifdef DEBUG
				  fprintf(numbers,"%6d %6d %6d %6d %12f %12f %20g %20g %12f %12f\n",i,ii,j,k,phase/2/PI,ampl,
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ],
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1], kmag, PowerSpectrum(kmag));

				  fprintf(numbers,"%6d %6d %6d %6d %12f %12f %20g %20g %12f %12f\n",i,ii,jj,k,phase/2/PI,ampl,
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)  ],
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)+1], kmag, PowerSpectrum(kmag));
#endif

				  /*  ORIGINAL
				    for(axes = 0; axes < 3; axes++)
				    {
				    cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				    -kvec[axes] / kmag2 * delta * sin(phase);
				    cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				    kvec[axes] / kmag2 * delta * cos(phase);
					  
				    cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].re =
				    -kvec[axes] / kmag2 * delta * sin(phase);
				    cdisp[axes][((i - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k].im =
				    -kvec[axes] / kmag2 * delta * cos(phase);
				    }
				  */
				}
			    }
			}
		      else	/* here comes i!=0 : conjugate can be on other processor! */
			{
			  if(i >= Nmesh / 2)
			    continue;
			  else
			    {
			      /* ii has already been assigned 
			      ii = Nmesh - i;
			      if(ii == Nmesh)
				ii = 0;
			      */
			      jj = Nmesh - j;
			      if(jj == Nmesh)
				jj = 0;

			      if(i >= Local_x_start && i < (Local_x_start + Local_nx))
				{

				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ] = delta * cos(phase);
				  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1] = delta * sin(phase);

#ifdef DEBUG				  
				  fprintf(numbers,"%6d %6d %6d %6d %12f %12f %20g %20g %12f %12f\n",i,ii,j,k,phase/2/PI,ampl,
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)  ],
					  kdensity[ThisGrid][2*(((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k)+1], kmag, PowerSpectrum(kmag));
#endif
				}

			      /*  ORIGINAL
				for(axes = 0; axes < 3; axes++)
				{
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].re =
				-kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((i - Local_x_start) * Nmesh + j) * (Nmesh / 2 + 1) + k].im =
				kvec[axes] / kmag2 * delta * cos(phase);
				}
			      */

			      if(ii >= Local_x_start && ii < (Local_x_start + Local_nx))
				{
				  
				  kdensity[ThisGrid][2*(((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)  ] = delta * cos(phase);
				  kdensity[ThisGrid][2*(((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)+1] = -delta * sin(phase);

#ifdef DEBUG
				  fprintf(numbers,"%6d %6d %6d %6d %12f %12f %20g %20g %12f %12f\n",ii,i,jj,k,phase/2/PI,ampl,
					  kdensity[ThisGrid][2*(((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)  ],
					  kdensity[ThisGrid][2*(((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) + k)+1], kmag, PowerSpectrum(kmag));
#endif
				}

			      /*  ORIGINAL
				for(axes = 0; axes < 3; axes++)
				{
				cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
				k].re = -kvec[axes] / kmag2 * delta * sin(phase);
				cdisp[axes][((ii - Local_x_start) * Nmesh + jj) * (Nmesh / 2 + 1) +
				k].im = -kvec[axes] / kmag2 * delta * cos(phase);
				}
			      */
			    }
			}
		    }
		}
	    }
	}
    }

#ifdef DEBUG
  fclose(numbers);
#endif

  /* This is needed to achieve proper normalization for pinocchio */

  fac=pow((double)Nmesh, 3.0);
  for (i=0; i<VectorLength; i++)
    kdensity[ThisGrid][i] *= fac;

  return 0;
}


double VarianceOnGrid(int ThisGrid, double Time, double ThisRadius)
{

  int i,j,k,ii,jj,kk;
  double fundamental = 2. * PI / MyGrids[ThisGrid].BoxSize;
  double nyquist = NYQUIST * MyGrids[ThisGrid].GSglobal_x/2 * fundamental;
  double d3k=pow(fundamental,3.);
  double w, D, kmag;
  double Variance=0.0;

  for (i=0; i<=MyGrids[ThisGrid].GSglobal_x; i++)
    {
      ii = i-MyGrids[ThisGrid].GSglobal_x/2;
      for (j=0; j<=MyGrids[ThisGrid].GSglobal_y; j++)
	{
	  jj = j-MyGrids[ThisGrid].GSglobal_y/2;
	  for (k=0; k<=MyGrids[ThisGrid].GSglobal_z; k++)
	    {
	      kk = k-MyGrids[ThisGrid].GSglobal_z/2;
	     
	      
	      kmag=sqrt((double)(ii*ii+jj*jj+kk*kk))*fundamental;
	      if (kmag>0.0 && kmag <= nyquist)
		{
		  w = WindowFunction(kmag * ThisRadius);
		  D = GrowingMode(1./Time-1.,kmag);
		  Variance += PowerSpectrum(kmag) * D*D * w*w * d3k;
		}
	    }
	}
    }

  Variance /= pow(2.*PI,3.);

  return Variance;
}
