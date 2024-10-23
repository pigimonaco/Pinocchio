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


// -------------------------
// defines and variables
// -------------------------

//#define DEBUG

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif

#define GRID MyGrids[ThisGrid]

typedef ptrdiff_t point[3];
point *starts;

// -------------------------
// functions' prototypes
// -------------------------


double greens_function                  (double *, double, int, int);
int    cubes_order                      (const void *, const void *);

// -------------------------
// code segment
// -------------------------


int cubes_order(const void *A, const void *B)
{
  float a = *(starts[*(int*)A]);
  float b = *(starts[*(int*)B]);

  return (a - b);

  /*
  if( a > b)
    return 1;
  if(a < b)
    return -1;
  return 0;
  */
}




int set_one_grid(int ThisGrid)
{
  ptrdiff_t    alloc_local;
  unsigned int pfft_flags;

  GRID.norm = (double)1.0 /
    ((double)GRID.Ntotal);
  
  GRID.CellSize = (double)GRID.BoxSize / GRID.GSglobal[_x_];

  //pfft_flags  = PFFT_MEASURE;
  pfft_flags  = 0;
  if(params.use_transposed_fft)
    pfft_flags  |= PFFT_TRANSPOSED_OUT;

  alloc_local = pfft_local_size_dft_r2c_3d(GRID.GSglobal, FFT_Comm, pfft_flags,
					   GRID.GSlocal, GRID.GSstart,
					   GRID.GSlocal_k, GRID.GSstart_k);

  
  dprintf(VDBG, ThisTask, "[set grid %02d] task %d %ld "
	  "i: %ld %ld %ld - i start: %ld %ld %ld - "
	  "o: %ld %ld %ld - o start: %ld %ld %ld\n",
	  ThisGrid, ThisTask, alloc_local,
	  GRID.GSlocal[_x_], GRID.GSlocal[_y_], GRID.GSlocal[_z_],
	  GRID.GSstart[_x_], GRID.GSstart[_y_], GRID.GSstart[_z_],
	  GRID.GSlocal_k[_x_], GRID.GSlocal_k[_y_], GRID.GSlocal_k[_z_],
	  GRID.GSstart_k[_x_], GRID.GSstart_k[_y_], GRID.GSstart_k[_z_]);
  
  
  GRID.total_local_size_fft = 2 * alloc_local;
  GRID.total_local_size     = GRID.GSlocal[_x_] * GRID.GSlocal[_y_] * GRID.GSlocal[_z_];  
  
  MyGrids[0].off = 0;
  // order the sub-blocks by row-major order (i, j, k), k first, then j, then i

  /* int           i; */
  /* unsigned int index; */

  /* starts         = (point*)malloc(sizeof(point) * NTasks); */
  
  /* MPI_Allgather(GRID.GSlocal, sizeof(point), MPI_BYTE, starts, sizeof(point), MPI_BYTE, MPI_COMM_WORLD); */
  /* for(i = 0; i < NTasks; i++) */
  /*   { */
  /*     index            = (starts[i][_x_]*GRID.GSglobal[_y_] + starts[i][_y_])*GRID.GSglobal[_z_] + starts[i][_z_]; */
  /*     starts[i][_x_]   = index; */
  /*     cubes_ordering[i] = i; */
  /*   } */

  /* qsort(cubes_ordering, NTasks, sizeof(int), cubes_order); */

  /* free(starts); */
  
  return 0;
}




int compute_fft_plans()
{
  ptrdiff_t DIM[3];
  int       pfft_flags;
  int       ThisGrid;
  
  dprintf(VMSG, 0, "[%s] Computing fft plans\n",fdate());

  for (ThisGrid = 0; ThisGrid < Ngrids; ThisGrid++)
    {

      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];

      /* create plan for out-of-place DFT */

      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE ;
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_OUT;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
      pfft_plan_with_nthreads(internal.nthreads_fft);
#endif
      GRID.forward_plan = pfft_plan_dft_r2c_3d(DIM, rvector_fft[ThisGrid], cvector_fft[ThisGrid],
					       FFT_Comm, PFFT_FORWARD, pfft_flags);


      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];

      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE;
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_IN;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
      pfft_plan_with_nthreads(internal.nthreads_fft);
#endif
      GRID.reverse_plan = pfft_plan_dft_c2r_3d(DIM, cvector_fft[ThisGrid], rvector_fft[ThisGrid], 
					       FFT_Comm, PFFT_BACKWARD, pfft_flags);
    }


  dprintf(VMSG, 0, "[%s] fft plans done\n",fdate());

  return 0;
}


double forward_transform(int ThisGrid)
{
  double time;
  
  time=MPI_Wtime();

  pfft_execute(GRID.forward_plan);

  return MPI_Wtime()-time;
}


double reverse_transform(int ThisGrid)
{

  int i;
  double time;

  time=MPI_Wtime();

  pfft_execute(GRID.reverse_plan);

  dvec         NORM   = {GRID.norm, GRID.norm, GRID.norm, GRID.norm};
  unsigned int mysize = GRID.total_local_size_fft / 4;

// #pragma GCC ivdep  
//   for (i = 0; i < mysize; i++)
//     rvector_fft[ThisGrid][i] *= NORM;
    
  for (i = GRID.total_local_size_fft - GRID.total_local_size_fft%4 ; i < GRID.total_local_size_fft; i++)
    rvector_fft[ThisGrid][i] *= GRID.norm;

  // non-vector code 
  for (i = 0 ; i < GRID.total_local_size_fft; i++)
      rvector_fft[ThisGrid][i] *= GRID.norm;

  return MPI_Wtime() - time;
}


int finalize_fft()
{
#ifndef RECOMPUTE_DISPLACEMENTS
  int ThisGrid;

  for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--)
    {
      pfft_destroy_plan(GRID.forward_plan);
    }

  /* for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--) */
  /*   { */
  /*     pfft_free(GRID.forward_plan); */
  /*     pfft_free(GRID.reverse_plan); */
  /*   } */

  pfft_cleanup();
  MPI_Comm_free(&FFT_Comm);
#endif

  return 0;
}


int compute_derivative(int ThisGrid, int first_derivative, int second_derivative)
{
  int    swap, local[3], start[3], C[3], N[3], Nhalf[3];
  double knorm[3];
#ifdef SCALE_DEPENDENT
  double k_module;
#endif
  
#ifdef DEBUG
  sprintf(filename,"results.%d-%d.%d",first_derivative,second_derivative,ThisTask);
  results=fopen(filename,"w");
#endif
  
  // note: for PFFT the transposed order is different when
  // the subdivision is done either in 1 or higher dimensions.
  // In 1D the memory order is [ y, x, z ], while in 2D and 3D
  // it is [ y, z, x].

  
  if(!params.use_transposed_fft)
    // this is the non-transposed order x, y, z    
    for(int i = 0; i < 3; i++)
      C[i] = i;
  else
    {
      if(internal.tasks_subdivision_dim > 1)
	C[0] = _y_, C[1] = _z_, C[2] = _x_;
      else
	C[0] = _y_, C[1] = _x_, C[2] = _z_;
    }

  for(int i = 0; i < 3; i++)
    {
      // grid number
      N[i]     = GRID.GSglobal[C[i]];
      // Nyquist frequencies
      Nhalf[i] = N[i] / 2;
      // k vectors
      knorm[i] = 2.*PI / (double)GRID.GSglobal[ C[i] ];
      
      local[i] = GRID.GSlocal_k[C[i]];
      start[i] = GRID.GSstart_k[C[i]];

    }

  // for first derivatives the real and imaginary parts must be swapped
  swap=((first_derivative==0 && second_derivative> 0) ||
	(first_derivative> 0 && second_derivative==0));

  double growth_rate=1.0;

  // loop over k-space indices

  int idx;
  for (idx = 0; idx < local[_x_]; idx++)
    {
      int    ii[3];
      
      ii[_x_] = idx + start[_x_]; // this is the kx index
      
      if (ii[_x_] > Nhalf[_x_])   // kx-N if kx > Nx/2
	ii[_x_] -= N[_x_];
      
      double k_x  = knorm[_x_] * ii[_x_]; 
      double k2_0 = k_x * k_x;

      int idy;
      for (idy = 0; idy < local[_y_]; idy++)
	{
	  ii[_y_] = idy + start[_y_]; // this is the ky index
	  
	  if (ii[_y_] > Nhalf[_y_])
	    ii[_y_] -= N[_y_];        // ky-N if ky > Ny/2
	  
	  double k_y  = knorm[_y_] * ii[_y_];
	  double k2_1 = k2_0 + k_y * k_y;

	  int idz;
	  for (idz = 0; idz < local[_z_]; idz++)
	    {
	      ii[_z_] = idz + start[_z_]; // this is the kz index
	      
	      if (ii[_z_] > Nhalf[_z_])
		ii[_z_] -= N[_z_];        // kz-N if ky > Nz/2
	      
	      double k_z  = knorm[_z_] * ii[_z_];

              double k_squared  = k2_1 + k_z * k_z;

#ifdef SCALE_DEPENDENT
	      double k_module = sqrt(k_squared);

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

	      int index = 2*(( idx * local[_y_] + idy ) * local[_z_] + idz);
	      
              if (k_squared != 0.)
		{

		  // Gaussian smoothing window
		  double smoothing = exp(-0.5 * k_squared * Rsmooth * Rsmooth);
		  
		  double diff_comp[4];
		  
		  // the components are stored in the vectors that control the differentiation
      
		  diff_comp[0] = 1.0;
		  diff_comp[1] = k_x;
		  diff_comp[2] = k_y;
		  diff_comp[3] = k_z;

		  double green = greens_function(diff_comp, k_squared, first_derivative, second_derivative);

		  (cvector_fft[ThisGrid][index/2])[0] *=  green * smoothing * growth_rate;
		  (cvector_fft[ThisGrid][index/2])[1] *=  green * smoothing * growth_rate;
		}

              if (swap)
		{
		  double tmp                          = (cvector_fft[ThisGrid][index/2])[1];
		  (cvector_fft[ThisGrid][index/2])[1] = (cvector_fft[ThisGrid][index/2])[0];
		  (cvector_fft[ThisGrid][index/2])[0] = -tmp;
		}
	    }
	}
    }

  /* int idy, idz, index; */
  /* for (idx = 0; idx < local[_x_]; idx++) */
  /*   for (idy = 0; idy < local[_y_]; idy++) */
  /*     for (idz = 0; idz < local[_z_]; idz++) */
  /* 	{ */
  /* 	  index = 2*(( idx * local[_y_] + idy ) * local[_z_] + idz); */
  /* 	  printf(" %3d %3d %3d  %6d  %8f %8f\n", */
  /* 		 idx,idy,idz,index,cvector_fft[ThisGrid][index/2][0],cvector_fft[ThisGrid][index/2][1]); */
  /* 	}   */

  /* char name[60]; */
  /* char coord[4] = {' ', 'x', 'y', 'z'}; */
  /* sprintf(name, "c_before_derivative.1.%c", coord[first_derivative]); */
  /* dump_cvector((double*)cvector_fft[ThisGrid], params.use_transposed_fft*internal.tasks_subdivision_dim, */
  /* 	       MyGrids[ThisGrid].GSglobal[_x_], */
  /* 	       MyGrids[ThisGrid].GSlocal_k, */
  /* 	       MyGrids[ThisGrid].GSstart_k, name, 0); */
  

  if (!ThisTask)
    dprintf(VMSG, 0, "[%s] compute_derivative: starting fft\n",fdate());

  {
    double time = reverse_transform(ThisGrid);
    
    if (!ThisTask)
      dprintf(VMSG, 0, "[%s] compute_derivative: done fft, cpu time = %f\n",fdate(),time);
    
    cputime.fft+=time;
  }

  /* for (int local_x = 0; local_x < MyGrids[ThisGrid].GSlocal[_x_]; local_x++) */
  /*   for (int local_y = 0; local_y < MyGrids[ThisGrid].GSlocal[_y_]; local_y++) */
  /*     for (int local_z = 0; local_z < MyGrids[ThisGrid].GSlocal[_z_]; local_z++) */
  /* 	{ */
  /* 	  index=local_x + MyGrids[ThisGrid].GSlocal[_x_] * (local_y + local_z * MyGrids[ThisGrid].GSlocal[_y_]); */
  /* 	  printf(" %3d %3d %3d  %6d  %8f\n",local_x,local_y,local_z,index,rvector_fft[ThisGrid][index]); */
  /* 	} */

  fflush(stdout);

  return 0;
}


double greens_function(double *diff_comp, double k_squared, int first_derivative, int second_derivative)
{

  if (first_derivative == -1 && second_derivative == -1)
    /* in this case the greens_function is simply 1 */
    return 1.0;

  if (first_derivative==0 && second_derivative==0)
    return -diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;
  else
    return diff_comp[first_derivative]*diff_comp[second_derivative]/k_squared;

}


void write_in_cvector(int ThisGrid, double * restrict vector)
{
  dvec * restrict target = (dvec*)cvector_fft[ThisGrid];
  dvec * restrict source = (dvec*)vector;
  int mysize = GRID.total_local_size_fft/DVEC_SIZE;

#if !defined(_OPENMP)  
#pragma GCC ivdep
#endif
#ifdef _OPENMP
#pragma omp for simd schedule(static)
#endif
  for ( int i = 0; i < mysize; i++ )
    *(target + i) = *(source + i);

  for (int i = mysize*DVEC_SIZE; i < GRID.total_local_size_fft; i++ )
    *((double*)cvector_fft[ThisGrid] + i) = *(vector+i);

  // non-vector code   
  /* for ( i = 0; i < GRID.total_local_size_fft; i++ ) */
  /*   *((double*)cvector_fft[ThisGrid] + i) = *(vector + i); */

}




void write_from_cvector(int ThisGrid, double * restrict vector)
{
  dvec * restrict source = (dvec*)cvector_fft[ThisGrid];
  dvec * restrict target = (dvec*)vector;

  int mysize = GRID.total_local_size_fft/DVEC_SIZE;

#if !defined(_OPENMP)
#pragma GCC ivdep
#endif
#ifdef _OPENMP
#pragma omp for simd schedule(static)  
#endif
  for ( int i = 0; i < mysize; i++ )
    *(target + i) = *(source + i);

  for (int i = mysize*DVEC_SIZE; i < GRID.total_local_size_fft; i++ )
    *(vector + i) = *((double*)cvector_fft[ThisGrid] + i);

  // non-vector code   
  /* for ( i = 0; i < GRID.total_local_size_fft; i++ ) */
  /*   *(vector + i) = *((double*)cvector_fft[ThisGrid] + i); */

}



void write_in_rvector(int ThisGrid, double * restrict vector)
{
  dvec * restrict target = (dvec*)rvector_fft[ThisGrid];
  dvec * restrict source = (dvec*)vector;
  int mysize = GRID.total_local_size / DVEC_SIZE;

#if !defined(_OPENMP)
#pragma GCC ivdep
#endif  
#ifdef _OPENMP
#pragma omp for simd schedule(static)  
#endif
  for( int i = 0; i < mysize; i++ )
    *(target + i) = *(source + i);
  
  for( int i = mysize*DVEC_SIZE ; i < GRID.total_local_size; i++ )
    *((double*)rvector_fft[ThisGrid] + i) = *(vector + i);

  // non-vector code   
  /* for(i = 0; i < GRID.total_local_size; i++) */
  /*   *(rvector_fft[ThisGrid] + i) = *(vector + i); */

}


void write_from_rvector(int ThisGrid, double * restrict vector)
{
  dvec * restrict source = (dvec*)rvector_fft[ThisGrid];
  dvec * restrict target = (dvec*)vector;
  int mysize = GRID.total_local_size / DVEC_SIZE;

#if !defined(_OPENMP)
#pragma GCC ivdep
#endif
#ifdef _OPENMP
#pragma omp for simd schedule(static)  
#endif
  for( int i = 0; i < mysize; i++ )
    *(target + i) = *(source + i);

  for( int i = mysize*DVEC_SIZE; i < GRID.total_local_size; i++ )
    *(vector + i) = *((double*)rvector_fft[ThisGrid] + i);

  // non-vector code 
  /* for( i = 0 ; i < GRID.total_local_size; i++) */
  /*   *(vector + i) = *(rvector_fft[ThisGrid] + i); */

}


int store_velocities()
{
  /* LUCA: vettorializziamo e ompizziamo? */

  /* loop on all particles */
  for (int index = 0; index < MyGrids[0].total_local_size; index++)
    for (int i = 0; i < 3; i++)
      products[index].Vel[i]=first_derivatives[0][i][index];

  return 0;
}


void dump_cvector(double * restrict vector, int T, int Nmesh, ptrdiff_t * restrict local_n, ptrdiff_t * restrict local_start,  char *filename, int ThisGrid)
// this routine dump in a unique files the 3d complex-space cube
// the quantity of each "particle" (i.e. cell) is written with its unique ordinal number in x,y,z
//   order
// - T          : flags whether the memory order is transposed or not, as follows:
//                   0: non transposed
//                   1: transposed as in 1D pfft order
//                 2-3: transposed as in 2d/3D pfft order
// - NMesh_2    : Ngrid/2
// - local_n    : the extension of this tasks domain in the 3 coordinates
// - local_start: the starting point of this task domain in the 3 coordinates
//
{

  int ii, i, j, k, addr, addr_;
  int Nmesh_2 = Nmesh / 2;
  FILE *file;
  char myfilename[300];
  
  sprintf(myfilename, "partial_%d", ThisTask);

  for(ii = 0; ii < NTasks; ii++)
    {
      if( (ThisTask == ii) &&
	  (local_start[_z_] < Nmesh_2-1) )	
	{
	  file = fopen(myfilename, "w");

	  if(T == 1)
	    // transposed order, 1D decomposition
	    {
	      for(j = 0; j < local_n[_y_]; j++)
	  	for(i = 0; i < local_n[_x_]; i++)
	  	  for(k = 0; k < (local_n[_z_]) && (k + local_start[_z_] <= Nmesh_2); k++)
	  	    {
		      if(internal.dump_kdensity == 1)
			// write in x,y order of the non-transposed case
			addr_   = ((i+local_start[_x_])*Nmesh + j+local_start[_y_])*(Nmesh_2+1) + k+local_start[_z_];
		      else
			// write in natural order
			addr_   = ((j+local_start[_y_])*Nmesh + i+local_start[_x_])*(Nmesh_2+1) + k+local_start[_z_];
	  	      addr    = 2*((j * local_n[_x_] + i) * local_n[_z_] + k);
	  	      fprintf(file, "%d - %g %g\n", addr_,
	  		      vector[addr], vector[addr+1]);
	  	    }
	      
	    }
	  else if(T > 1)
	    // transposed order, 2D or 3D decomposition
	    {
	      for(j = 0; j < local_n[_y_]; j++)
	  	for(k = 0; k < (local_n[_z_]) && (k + local_start[_z_] <= Nmesh_2); k++)
	  	  for(i = 0; i < local_n[_x_]; i++)
	  	    {
		      if(internal.dump_kdensity == 1)
			// write in x,y order of the non-transposed case
			addr_   = ((i+local_start[_x_])*Nmesh + j+local_start[_y_])*(Nmesh_2+1) + k+local_start[_z_];
		      else
			// write in natural order
			addr_   = ((j+local_start[_y_])*Nmesh + k+local_start[_z_])*(Nmesh_2+1) + i+local_start[_x_];
	  	      addr    = 2*((j * local_n[_z_] + k) * local_n[_x_] + i);
	  	      fprintf(file, "%d - %g %g\n", addr_,
	  		      vector[addr], vector[addr+1]);
	  	    }
	      
	    }
	  else
	    // not transposed
	    {
	      for(i = 0; i < local_n[_x_]; i++)
	  	for(j = 0; j < local_n[_y_]; j++)
	  	  for(k = 0; k < (local_n[_z_]) && (k + local_start[_z_] <= Nmesh_2); k++)
	  	    {
		      if(internal.dump_kdensity == 2)
			// write in x,y order of the transposed case (1-D decomposition)
			addr_   = ((j+local_start[_y_])*Nmesh + i+local_start[_x_])*(Nmesh_2+1) + k+local_start[_z_];
		      else
			// write in natural order
			addr_   = ((i+local_start[_x_])*Nmesh + j+local_start[_y_])*(Nmesh_2+1) + k+local_start[_z_];
	  	      addr    = 2*( (i*local_n[_y_] + j)*local_n[_z_] + k);
	  	      fprintf(file, "%d - %g %g\n", addr_,
	  		      vector[addr], vector[addr+1]);
	  	    }
	      
	    }

	  fflush(file);
	  fclose(file);
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if(ThisTask == 0)
    {
      int err;
      sprintf(myfilename, "cat partial_* | sort -k1 -n > %s", filename);
      err = system(myfilename);
      err = system("rm -f partial_*");
      if(err != 0)
	dprintf(VXX, 0, "something got wrong while writing cvector dump file\n");
    }
     
  
  return;

}


void dump_rvector(double * restrict vector, int Nmesh, ptrdiff_t * restrict local_n, ptrdiff_t * restrict local_start,  char *filename, int ThisGrid)
// this routine dump in a unique files the 3d complex-space cube
// the quantity of each "particle" (i.e. cell) is written with its unique ordinal number in x,y,z
//   order
// - NMesh_2    : Ngrid/2
// - local_n    : the extension of this tasks domain in the 3 coordinates
// - local_start: the starting point of this task domain in the 3 coordinates
//
{

  int ii, i, j, k, addr, addr_;
  //int Nmesh_2 = Nmesh / 2;
  FILE *file;
  char myfilename[300];
  
  sprintf(myfilename, "partial_%d", ThisTask);

  for(ii = 0; ii < NTasks; ii++)
    {
      if( (ThisTask == ii) )
	{
	  file = fopen(myfilename, "w");
	  
	  for(i = 0; i < local_n[_x_]; i++)
	    for(j = 0; j < local_n[_y_]; j++)
	      for(k = 0; k < local_n[_z_]; k++)
		{
		  /* if (params.use_transposed_fft) */
		  /*   { */
		  /*     addr_   = ((j+local_start[_y_])*Nmesh + k+local_start[_z_])*Nmesh + i+local_start[_x_]; */
		  /*     addr    = ( (k*local_n[_y_] + j)*local_n[_x_] + i); */
		  /*   } */
		  /* else */
		    { // OK, AGGIUDICATO
		      addr_   = ((i+local_start[_x_])*Nmesh + j+local_start[_y_])*Nmesh + k+local_start[_z_];
		      addr    = ( (i*local_n[_y_] + j)*local_n[_z_] + k);
		    }

		  fprintf(file, "%d  %d %d %d  %g\n", addr_, 
			  (int)(i + local_start[_x_]), (int)(j + local_start[_y_]), (int)(k + local_start[_z_]),
			  vector[addr]);
		}

	  fflush(file);
	  fclose(file);
	}

      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if(ThisTask == 0)
    {
      int err;
      sprintf(myfilename, "cat partial_* | sort -k1 -n > %s", filename);
      err = system(myfilename);
      err = system("rm -f partial_*");
      if(err != 0)
	dprintf(VXX, 0, "something got wrong while writing cvector dump file\n");
    }
     
  
  return;

}
