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

//#define DEBUG

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif


double greens_function(double *, double, int, int);

int set_one_grid(int ThisGrid)
{
  ptrdiff_t L, M, N;
  ptrdiff_t alloc_local, local_n0, local_0_start;
  ptrdiff_t local_n1, local_1_start;

  fftw_mpi_init();

  MyGrids[ThisGrid].norm = 1.0
    /(double)MyGrids[ThisGrid].GSglobal_x
    /(double)MyGrids[ThisGrid].GSglobal_y
    /(double)MyGrids[ThisGrid].GSglobal_z;
  MyGrids[ThisGrid].CellSize = MyGrids[ThisGrid].BoxSize/(double)MyGrids[ThisGrid].GSglobal_x;

  L=MyGrids[ThisGrid].GSglobal_x;
  M=MyGrids[ThisGrid].GSglobal_y;
  N=MyGrids[ThisGrid].GSglobal_z;

  /* get local data size and allocate */
  alloc_local = fftw_mpi_local_size_3d_transposed(L, M, N/2+1, MPI_COMM_WORLD,
						  &local_n0, &local_0_start,
						  &local_n1, &local_1_start);

  MyGrids[ThisGrid].GSlocal_x = MyGrids[ThisGrid].GSglobal_x;
  MyGrids[ThisGrid].GSlocal_y = MyGrids[ThisGrid].GSglobal_y;
  MyGrids[ThisGrid].GSlocal_z = local_n0;

  MyGrids[ThisGrid].GSstart_x = 0;
  MyGrids[ThisGrid].GSstart_y = 0;
  MyGrids[ThisGrid].GSstart_z = local_0_start;

  MyGrids[ThisGrid].GSlocal_k_x = MyGrids[ThisGrid].GSglobal_x;
  MyGrids[ThisGrid].GSlocal_k_y = local_n1;
  MyGrids[ThisGrid].GSlocal_k_z = MyGrids[ThisGrid].GSglobal_z;

  MyGrids[ThisGrid].GSstart_k_x = 0;
  MyGrids[ThisGrid].GSstart_k_y = local_1_start;
  MyGrids[ThisGrid].GSstart_k_z = 0;

  MyGrids[ThisGrid].total_local_size_fft = 2*alloc_local;
  MyGrids[ThisGrid].total_local_size = MyGrids[ThisGrid].GSlocal_x * MyGrids[ThisGrid].GSlocal_y * MyGrids[ThisGrid].GSlocal_z;

  /* this gives is the number of unused x-grid points in rvector_fft */
  if (MyGrids[0].GSglobal_x%2)
    MyGrids[0].off=1;
  else
    MyGrids[0].off=2;

  return 0;
}

int initialize_fft()
{
  ptrdiff_t L, M, N;
  int igrid;

  for (igrid=0; igrid<Ngrids; igrid++)
    {

      L=MyGrids[igrid].GSglobal_x;
      M=MyGrids[igrid].GSglobal_y;
      N=MyGrids[igrid].GSglobal_z;

      /* create plan for out-of-place r2c DFT */
      MyGrids[igrid].forward_plan = fftw_mpi_plan_dft_r2c_3d(L, M, N, rvector_fft[igrid], cvector_fft[igrid], 
							     MPI_COMM_WORLD, FFTW_MEASURE + FFTW_MPI_TRANSPOSED_OUT);

      MyGrids[igrid].reverse_plan = fftw_mpi_plan_dft_c2r_3d(L, M, N, cvector_fft[igrid], rvector_fft[igrid], 
							     MPI_COMM_WORLD, FFTW_MEASURE + FFTW_MPI_TRANSPOSED_IN);
    }

  return 0;
}


double forward_transform(int ThisGrid)
{
  double time;
  
  time=MPI_Wtime();

  fftw_execute(MyGrids[ThisGrid].forward_plan);

  return MPI_Wtime()-time;
}


double reverse_transform(int ThisGrid)
{

  int i;
  double time;
  // double ave, var; //levare

  time=MPI_Wtime();

  fftw_execute(MyGrids[ThisGrid].reverse_plan);

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
      rvector_fft[ThisGrid][i]*=MyGrids[ThisGrid].norm;

#ifdef DEBUG
  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    fprintf(results, " %d  %g \n", i, rvector_fft[ThisGrid][i]);
  fclose(results);
#endif

  return MPI_Wtime()-time;
}


int finalize_fftw()
{
#ifndef RECOMPUTE_DISPLACEMENTS
  int igrid;

  for (igrid=Ngrids-1; igrid>=0; igrid--)
    {
      fftw_destroy_plan(MyGrids[igrid].forward_plan);
      fftw_destroy_plan(MyGrids[igrid].reverse_plan);
    }

  for (igrid=Ngrids-1; igrid>=0; igrid--)
    if (deallocate_fft_vectors(igrid))
      return 1;

  fftw_cleanup();
#endif
  return 0;
}



int compute_derivative(int ThisGrid, int first_derivative, int second_derivative)
{
  int swap, local_x,local_y,local_z,ixx,iyy,izz,index,nxhalf,nyhalf,nzhalf;
  double kx,ky,kz,kxnorm,kynorm,kznorm,green,k_squared,smoothing,tmp,time;
  double diff_comp[4];
#ifdef SCALE_DEPENDENT
  double k_module;
#endif

#ifdef DEBUG
  sprintf(filename,"results.%d-%d.%d",first_derivative,second_derivative,ThisTask);
  results=fopen(filename,"w");
#endif

  /* k vectors */
  kxnorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_x;
  kynorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_y;
  kznorm   = 2.*PI/(double)MyGrids[ThisGrid].GSglobal_z;

  /* Nyquist frequencies */
  nxhalf = MyGrids[ThisGrid].GSglobal_x/2;
  nyhalf = MyGrids[ThisGrid].GSglobal_y/2;
  nzhalf = MyGrids[ThisGrid].GSglobal_z/2;

  /* for first derivatives the real and imaginary parts must be swapped */
  swap=((first_derivative==0 && second_derivative> 0) ||
	(first_derivative> 0 && second_derivative==0));

/*
  NB: ix, iy and iz are global coordinates of the box (from 0 to N-1)
      ix: [0,N/2]
      iy, iz: [0,N-1]
      iyy and izz are unfolded, they run from -N/2+1 to N/2 (ixx = ix)
      local_x, local_y and local_z are local coordinates of the slice
      (for the x and z coordinates they are the same as ix and iz)
*/

  double growth_rate=1.0;

/* loop over k-space indices */
/* This loop is correct for transposed order in FFTW */

  for (local_z = 0; local_z < MyGrids[ThisGrid].GSlocal_k_z; local_z++)
    {
      izz = local_z + MyGrids[ThisGrid].GSstart_k_z;
      if (local_z > nzhalf)
	izz -= MyGrids[ThisGrid].GSglobal_z;
      kz  = kznorm*izz;

      for (local_y = 0; local_y < MyGrids[ThisGrid].GSlocal_k_y; local_y++)
	{
	  iyy = local_y + MyGrids[ThisGrid].GSstart_k_y;
	  if (iyy > nyhalf)
	    iyy -= MyGrids[ThisGrid].GSglobal_y;
	  ky  = kynorm*iyy;

	  for (local_x = 0; local_x <= nxhalf; local_x++)
	    {
	      ixx = local_x;
	      kx  = kxnorm*ixx;

              k_squared  = kx*kx + ky*ky + kz*kz;

#ifdef SCALE_DEPENDENT
	      k_module = sqrt(k_squared);

	      /* In the scale-dependent case the delta(k) must be multiplied 
		 by the relevant growth rate */
	      switch (ScaleDep.order)
		{
		case 0:
		  growth_rate = 1.0;
		  break;
		case 1:
		  growth_rate = GrowingMode(ScaleDep.redshift,k_module);
		  break;
		case 2:
		  growth_rate = GrowingMode_2LPT(ScaleDep.redshift,k_module);
		  break;
		case 3:
		  growth_rate = GrowingMode_3LPT_1(ScaleDep.redshift,k_module);
		  break;
		case 4:
		  growth_rate = GrowingMode_3LPT_2(ScaleDep.redshift,k_module);
		  break;
		default:
		  growth_rate = 1.0;
		  break;
		}
#endif

	      /* corresponding index of real part in vector (imaginary in index + 1) */
	      index = 1 + 2*local_x + (MyGrids[ThisGrid].GSglobal_x+MyGrids[ThisGrid].off)
		*(local_z + local_y* MyGrids[ThisGrid].GSglobal_z);

              if (k_squared != 0.)
		{

		  /* Gaussian smoothing window */
		  smoothing = exp(-0.5 * k_squared * Rsmooth * Rsmooth);

		  /* the components are stored in the vectors that control the differentiation */
		  diff_comp[0] = 1.0;
		  diff_comp[1] = kx;
		  diff_comp[2] = ky;
		  diff_comp[3] = kz;

		  green = greens_function(diff_comp, k_squared, first_derivative, second_derivative);

		  (cvector_fft[ThisGrid][index/2])[0] *=  green * smoothing * growth_rate;
		  (cvector_fft[ThisGrid][index/2])[1] *=  green * smoothing * growth_rate;
		}

              if (swap)
		{
		  tmp                                 = (cvector_fft[ThisGrid][index/2])[1];
		  (cvector_fft[ThisGrid][index/2])[1] = (cvector_fft[ThisGrid][index/2])[0];
		  (cvector_fft[ThisGrid][index/2])[0] = -tmp;
		}
	    }
	}
    }

  if (!ThisTask)
    printf("[%s] compute_derivative: starting fft\n",fdate());

  time=reverse_transform(ThisGrid);

  if (!ThisTask)
    printf("[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);

  cputime.fft+=time;

  return 0;
}


double greens_function(double *diff_comp, double k_squared, int first_derivative, int second_derivative)
{

  /* in this case the greens_function is simply 1 */
  if (first_derivative == -1 && second_derivative == -1)
    return 1.0;

  if (first_derivative==0 && second_derivative==0)
    return -diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;
  else
    return diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;

}


void write_in_cvector(int ThisGrid, double *vector)
{
  int i;

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    *((double*)cvector_fft[ThisGrid]+i)=*(vector+i);

}


void write_from_cvector(int ThisGrid, double *vector)
{
  int i;

  for (i=0; i<MyGrids[ThisGrid].total_local_size_fft; i++)
    *(vector+i)=*((double*)cvector_fft[ThisGrid]+i);

}


void write_in_rvector(int ThisGrid, double *vector)
{
  int lx,ly,lz;
  
  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) =
	    *(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y));

  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=MyGrids[ThisGrid].GSlocal_x; lx<MyGrids[ThisGrid].GSlocal_x+MyGrids[ThisGrid].off; lx++)
	*(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y)) = 0.0;

}


void write_from_rvector(int ThisGrid, double *vector)
{
  int lx,ly,lz;

  for (lz=0; lz<MyGrids[ThisGrid].GSlocal_z; lz++)
    for (ly=0; ly<MyGrids[ThisGrid].GSlocal_y; ly++)
      for (lx=0; lx<MyGrids[ThisGrid].GSlocal_x; lx++)
	*(vector + lx + MyGrids[ThisGrid].GSlocal_x *(ly + lz* MyGrids[ThisGrid].GSlocal_y)) =
          *(rvector_fft[ThisGrid] + lx + (MyGrids[ThisGrid].GSlocal_x + MyGrids[ThisGrid].off) * (ly + lz * MyGrids[ThisGrid].GSlocal_y));

}
