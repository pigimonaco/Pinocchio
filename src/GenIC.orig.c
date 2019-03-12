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

//#define DPRINT

#ifdef DPRINT
#define Dprintf(...) printf(__VA_ARGS__);
#else
#define Dprintf(...)
#endif


#define GRID MyGrids[ThisGrid]

int GenIC(int ThisGrid)
{
  int           i, j, k;
  int           Nmesh, Nsample, VectorLength;
  int           local_n[3], local_start[3];
  double        Box, fac;
  unsigned int *SEEDTABLE;

  SEEDTABLE = (unsigned int*) malloc( sizeof(unsigned int) * GRID.GSglobal[_x_] * GRID.GSglobal[_y_]);

  
  /* local_n[0]     = GRID.GSlocal_k[_y_]; */
  /* local_start[0] = GRID.GSstart_k[_y_]; */
  
  /* local_n[1]     = GRID.GSlocal_k[_x_]; */
  /* local_start[1] = GRID.GSstart_k[_x_]; */

  local_n[0]     = GRID.GSlocal_k[_x_];
  local_start[0] = GRID.GSstart_k[_x_];
  
  local_n[1]     = GRID.GSlocal_k[_y_];
  local_start[1] = GRID.GSstart_k[_y_];

  local_n[2]     = GRID.GSlocal_k[_z_];
  local_start[2] = GRID.GSstart_k[_z_];
  
  VectorLength   = GRID.total_local_size_fft;
  Nmesh          = GRID.GSglobal[_x_];
  Nsample        = GRID.GSglobal[_x_];
  Box            = GRID.BoxSize;

  fac = pow(1./Box,1.5);
  
  Dprintf(" Task %d: %d -> %d | %d -> %d | %d -> %d | %d %g | %g %g %g %g\n", ThisTask,
	  local_start[0], local_start[0]+local_n[0],
	  local_start[1], local_start[1]+local_n[1],
	  local_start[2], local_start[2]+local_n[2],
	  VectorLength, fac,
	  PowerSpectrum(0.0001),
	  PowerSpectrum(0.001),
	  PowerSpectrum(0.01),
	  PowerSpectrum(0.1));
  
  gsl_rng_set(random_generator, params.RandomSeed);

  for(i = 0; i < Nmesh / 2; i++)
    {
      for(j = 0; j < i; j++)
	SEEDTABLE[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	SEEDTABLE[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	SEEDTABLE[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	SEEDTABLE[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	SEEDTABLE[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	SEEDTABLE[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i; j++)
	SEEDTABLE[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);

      for(j = 0; j < i + 1; j++)
	SEEDTABLE[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
    }
  
  gsl_rng_set(random_generator, params.RandomSeed);

  // the vector is initialized to zero
  for (i = 0; i < VectorLength; i++)
    kdensity[ThisGrid][i]=0.;


  int    ii, jj, kk, RR;
  int    iii, jjj;
  int    addr, addr_j, Nmesh_2, Nmesh_odd;
  double sign;
  double kvec[3];
  double delta, phase, ampl;

  double p_of_k, kmag, kmag2_i, kmag2_ij, kmag2_local;
  
  gsl_rng *k0_generator;
  
  Nmesh_2   = Nmesh / 2;
  Nmesh_odd = (Nmesh % 2);  

  // a second random generator, used on the plane k = 0
  
  k0_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  
  // this is the main loop on grid modes
  for(i = 0; i < local_n[_x_]; i++)
    {

      // ii is the x grid coordinate
      ii = local_start[_x_] + i;
      if(ii == Nmesh_2)
	continue;
      
      if(ii < Nmesh_2)
	kvec[0] = ii * 2 * PI / Box;
      else
	kvec[0] = -(Nmesh - ii) * 2 * PI / Box;

      kmag2_i = kvec[0] * kvec[0];
      iii = ii;
      
      
      for(j = 0; j < local_n[_y_]; j++)	    
	{

	  // jj is the y grid coordinate
	  jj = local_start[_y_] + j;
	  if(jj == Nmesh_2)
	    continue;
	  jjj = jj;
	  
	  if(jj < Nmesh_2)
	    kvec[1] = jj * 2 * PI / Box;
	  else
	    kvec[1] = -(Nmesh - jj) * 2 * PI / Box;

	  kmag2_ij = kmag2_i + kvec[1] * kvec[1];


	  // initialize the random generator with the seed found in the
	  // seedtable at the (i,j) coordinate in grid coordinates
	  gsl_rng_set(random_generator, SEEDTABLE[ii * Nmesh + jj]);


	  // since the chain of random numbers is subsequent, to
	  // reproduce all the time the same sequence, this
	  // process must generate all the random numbers on the
	  // k column that are below its starting k coordinate
	  if(local_start[_z_] > 0)
	    for(RR = 0; RR < local_start[_z_]; RR++)
	      {
		phase = gsl_rng_uniform(random_generator);
		do
		  phase = gsl_rng_uniform(random_generator);
		while(phase == 0);
	      }

	  // inner loop
	  for(k = 0; (k < local_n[_z_]) && (k+local_start[_z_] < Nmesh_2); k++)
	    {
	      // kk is the z grid coordinate
	      kk = local_start[_z_] + k;

	      // generate phase and amplitude
	      phase = gsl_rng_uniform(random_generator) * 2 * PI;
	      do
		ampl = gsl_rng_uniform(random_generator);
	      while(ampl == 0);	      	      

	      // blind points on the cube
	      if(ii == 0 && jj == 0 && kk == 0)
		continue;
	      if(kk == Nmesh_2)
		continue;
	      
	      if(kk < Nmesh_2)
		kvec[2] = kk * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - kk) * 2 * PI / Box;

	      kmag2_local = kmag2_ij + kvec[2] * kvec[2];
	      kmag = sqrt(kmag2_local);

	      if(kmag * Box / (2 * PI) > NYQUIST * Nsample / 2)
		continue;		  	      

	      p_of_k = PowerSpectrum(kmag);	      

	      sign = 1.0;	      
	      addr_j = j;	     

	      // special simmetries on the plane k = 0
	      if( kk == 0 )
		{

		  // some points are left empty
		  if( (ii == 0) && (jj == Nmesh_2 || jj == Nmesh_2 + Nmesh_odd) )
		    continue;
		  if( (ii == Nmesh_2 + Nmesh_odd) )
		    continue;

		  // the helf-plane for ii > Nmesh_2 takes the seeds
		  // from points in the other half-plane
		  if( (ii > Nmesh_2) ||
		      ( ii == 0 && jj > Nmesh/2) )
		    {
		      // simmetries on j
		      jjj = Nmesh -jj;
		      if(jjj == Nmesh)
			jjj = 0;

		      if(Nmesh_odd && jj == Nmesh_2+1)
			{
			  jjj = Nmesh_2+1;
			  addr_j = Nmesh_2;
			}

		      // simmetries on i
		      if(ii > Nmesh_2)
			iii = Nmesh - ii;
		      
		      sign = -1.0;

		      // re-initialize the spare random generator to the right seed
		      gsl_rng_set(k0_generator, SEEDTABLE[iii * Nmesh + jjj]);
		      phase = gsl_rng_uniform(k0_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(k0_generator);
		      while(ampl == 0);
		    }
		  Dprintf("\t (+) @ %d :: %d %d %d :: %u :: kmag: %g ph: %g A: %g D: %g P: %g PS: %g\n",
		  	 ThisTask, ii, jj, kk, SEEDTABLE[iii*Nmesh + jjj],
		  	 kmag, phase, ampl, delta, p_of_k, PowerSpectrum(kmag)); fflush(stdout);
		}
	      else
		{
		  Dprintf("\t (+) @ %d :: %d %d %d :: %u :: kmag: %g ph: %g A: %g D: %g P: %g PS: %g\n",
			  ThisTask, ii, jj, kk, SEEDTABLE[ii * Nmesh + jj],
			  kmag, phase, ampl, delta, p_of_k, PowerSpectrum(kmag)); fflush(stdout);
		}

#ifndef NO_RANDOM_MODULES
	      p_of_k *= -log(ampl);
#endif
	      delta = fac * sqrt(p_of_k);

	      // calculate the storing address in local coordinates
	      addr = 2*((i * local_n[_y_] + addr_j) * local_n[_z_] + k);
	      
	      kdensity[ThisGrid][addr ] =  delta * cos(phase);
	      kdensity[ThisGrid][addr + 1] =  sign * delta * sin(phase);

	    }

	  // since the chain of random numbers is subsequent, [..see above..]
	  // generate all the random numbers on the
	  // k column that are above its final k coordinate
	  if(local_start[_z_]+local_n[_z_] < Nmesh_2)
	    for(RR = local_start[_z_]+local_n[_z_]; RR <= Nmesh_2; RR++)
	      {
		phase = gsl_rng_uniform(random_generator);
		do
		  phase = gsl_rng_uniform(random_generator);
		while(phase == 0);
	      }	  
	}
    }

  gsl_rng_free(k0_generator);
  free(SEEDTABLE);
  
#ifdef DEBUG
  fclose(numbers);
#endif


  int addr_;

  FILE *file;
  char filename[100];

  sprintf(filename, "partial_%d", ThisTask);

  for(ii = 0; ii < NTasks; ii++)
    {
      if(ThisTask == ii)
	{
	  file = fopen(filename, "w");
	  for(i = 0; i < local_n[_x_]; i++)
	    for(j = 0; j < local_n[_y_]; j++)
	      for(k = 0; k < (local_n[_z_]-1); k++)
		{
		  addr = 2*( (i*local_n[_y_] + j)*local_n[_z_] + k);
		  addr_ = ((i+local_start[_x_])*Nmesh + j+local_start[_y_])*(Nmesh_2+1) + k+local_start[_z_];
		  
		  /* fprintf(file, "%d %d %d %d @ %d - %g %g\n", addr_, */
		  /* 	  i+local_start[_x_], */
		  /* 	  j+local_start[_y_], */
		  /* 	  k+local_start[_z_], */
		  /* 	  ThisTask, kdensity[ThisGrid][addr], kdensity[ThisGrid][addr+1]); */
		  fprintf(file, "%d - %g %g\n", addr_, kdensity[ThisGrid][addr], kdensity[ThisGrid][addr+1]);
		}	      
	  fflush(stdout);
	  fclose(file);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  if(ThisTask == 0)
    {
      system("cat partial_* | sort -k1 -n > my.dkdensity");
      system("rm -f partial_*");	  
    }

  
  /* This is needed to achieve proper normalization for pinocchio */

  fac=pow((double)Nmesh, 3.0);
  for (i=0; i<VectorLength; i++)
    kdensity[ThisGrid][i] *= fac;

  return 0;
}


