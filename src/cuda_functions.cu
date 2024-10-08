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

#if defined(CUFFTMP)

#include "pinocchio.h"
#include <cuda_runtime.h>
#include <complex.h>
#include <cuComplex.h>
#include <cufftMp.h>

// -------------------------
// defines and variables
// -------------------------

//#define DEBUG

cufftDoubleReal        *real_grid;
cufftDoubleComplex *complex_grid;
cudaStream_t stream{};
//cufftHandle forward_plan;
//cufftHandle reverse_plan;
cudaLibXtDesc *forward_desc;
cudaLibXtDesc *reverse_desc;
MPI_Comm comm = MPI_COMM_WORLD;

#ifdef DEBUG
  FILE *results;
  char filename[300];
#endif

#define GRID MyGrids[ThisGrid]
#define int64 long long int

struct Box3D {
    int64 lower[3];
    int64 upper[3];
    int64 strides[3];
};

// -------------------------
// functions' prototypes
// -------------------------


//double greens_function                  (double *, double, int, int);
//int    cubes_order                      (const void *, const void *);

// -------------------------
// code segment
// -------------------------

/*
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
 
}
*/

/*
GIOVANNI
Redefine these functions involving FFT calls in a .cu extension file in which we can use
the NVIDIA cufftMP library. The reason is that cufftMP needs some C++ abstract descriptors
which don't have any counterpart in standard C
*/

/*
Write a function to assign a device to each MPI task
This may become redundant since it will be called by each function,
but for these first steps it's necessary so we know which task 
is using which GPU. This step is planned to be skipped in the future
*/

/*C++ extensions*/
auto make_box = [](int64 lower[3], int64 upper[3], int64 strides[3]) {
		  Box3D box;
		  for(int i = 0; i < 3; i++) {
		    box.lower[i] = lower[i];
		    box.upper[i] = upper[i];
		    box.strides[i] = strides[i];
		  }
		  return box;
		};


void cuda_init()
{
  int ndevices;
  cudaGetDeviceCount(&ndevices);
  cudaSetDevice(ThisTask % ndevices);

  if (ThisTask == 0)
    if (ndevices == 0){
      printf("No accelerators found!");
      return;
    }
}

  
int set_one_grid(int ThisGrid)
{

  cuda_init();
  
  cufftResult_t status; //Check whether cuda functions are returning
  
  //ptrdiff_t    alloc_local;
  //unsigned int pfft_flags;

  //Create the FFT cuda stream 
  //cudaStreamCreate(&stream);
  
  //Create the plans
  //cufftHandle forward_plan;
  //cufftHandle reverse_plan;

  /*
  status = cufftCreate(&GRID.forward_plan);
  if (status != CUFFT_SUCCESS) {printf("!!! forward plan cufftCreate ERROR %d !!!\n", status);}
  
  status = cufftCreate(&GRID.reverse_plan);
  if (status != CUFFT_SUCCESS) {printf("!!! reverse plan cufftCreate ERROR %d !!!\n", status);}
  */
  
  GRID.norm = (double)1.0 /
    ((double)GRID.Ntotal);
  
  GRID.CellSize = (double)GRID.BoxSize / GRID.GSglobal[_x_];

  /*Initialization*/
  int64 nx               = (int64)GRID.GSglobal[_x_];
  int64 ny               = (int64)GRID.GSglobal[_y_];
  int64 nz               = (int64)GRID.GSglobal[_z_];
  int64 nz_real          = nz;
  int64 nz_complex       = (((int64)(GRID.GSglobal[_z_])/2)+1);
  int64 nz_real_padded   = 2*nz_complex;

  Box3D box_real, box_complex;
  
  {
    // Input data are X-slabs. 
    // Strides are packed and in-place (i.e., real is padded)
    int64 lower[3]   = {nx / NTasks * (ThisTask),   0,  0};
    int64 upper[3]   = {nx / NTasks * (ThisTask+1), ny, nz_real};
    int64 strides[3] = {(upper[1]-lower[1])*nz_real_padded, nz_real_padded, 1};
    box_real = make_box(lower, upper, strides);
  }
  
  {
    // Output data are Y-slabs.
    // Strides are packed
    int64 lower[3]   = {0,  ny / NTasks * (ThisTask),   0};
    int64 upper[3]   = {nx, ny / NTasks * (ThisTask+1), nz_complex};
    int64 strides[3] = {(upper[1]-lower[1])*(upper[2]-lower[2]), (upper[2]-lower[2]), 1};
    box_complex = make_box(lower, upper, strides);
  }
  /*End of initialization*/
  
  //MPI_Comm comm = MPI_COMM_WORLD;
  
  //Create the FFT cuda stream and 
  //attach the MPI default communicator to the distributed FFT
  /*
  status = cufftMpAttachComm(GRID.forward_plan, CUFFT_COMM_MPI, &comm);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftMpAttachComm ERROR %d !!!\n", status);}

  status = cufftSetStream(GRID.forward_plan, stream);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftSetStream ERROR %d !!!\n", status);}

  status = cufftMpAttachComm(GRID.reverse_plan, CUFFT_COMM_MPI, &comm);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftMpAttachComm ERROR %d !!!\n", status);}

  status = cufftSetStream(GRID.reverse_plan, stream);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftSetStream ERROR %d !!!\n", status);}
    
  
  printf("[DBG CHECK START] %ld, %ld, %ld, %lld, %lld, %lld\n", GRID.GSstart[_x_], GRID.GSstart[_y_], GRID.GSstart[_z_],
	 real_lower[0], real_lower[1], real_lower[2]);
  printf("[DBG CHECK START_K] %ld, %ld, %ld, %lld, %lld, %lld\n", GRID.GSstart_k[_x_], GRID.GSstart_k[_y_], GRID.GSstart_k[_z_],
	 c_lower[0], c_lower[1], c_lower[2]);
  printf("[DBG CHECK LOCAL] %ld, %ld, %ld, %lld, %lld, %lld\n", GRID.GSlocal[_x_], GRID.GSlocal[_y_], GRID.GSlocal[_z_],
	 real_upper[0], real_upper[1], real_upper[2]);
  printf("[DBG CHECK LOCAL_K] %ld, %ld, %ld, %lld, %lld, %lld\n", GRID.GSlocal_k[_x_], GRID.GSlocal_k[_y_], GRID.GSlocal_k[_z_],
	 c_upper[0], c_upper[1], c_upper[2]);
  
  printf("Numbers: %lld, %lld\n", nz_complex, nz_real_padded);
  */

  /*
  //Distribute r2c plan
  status = cufftXtSetDistribution(GRID.forward_plan, 3, box_real.lower, box_real.upper, box_complex.lower,
				  box_complex.upper, box_real.strides, box_complex.strides);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtSetDistribution r2c ERROR %d !!!\n", status);}

  //Distribute c2r plan
  status = cufftXtSetDistribution(GRID.reverse_plan, 3, box_real.lower, box_real.upper, box_complex.lower,
				  box_complex.upper, box_real.strides, box_complex.strides);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtSetDistribution c2r ERROR %d !!!\n", status);}
  */
  //pfft_flags  = PFFT_MEASURE;
  //pfft_flags  = 0;
  //if(params.use_transposed_fft)
  // pfft_flags  |= PFFT_TRANSPOSED_OUT;

  /*
  alloc_local = pfft_local_size_dft_r2c_3d(GRID.GSglobal, FFT_Comm, pfft_flags,
					   GRID.GSlocal, GRID.GSstart,
					   GRID.GSlocal_k, GRID.GSstart_k);
  */
  
  dprintf(VDBG, ThisTask, "[set grid %02d] task %d %lld "
	  "i: %lld %lld %lld - i start: %lld %lld %lld - "
	  "o: %lld %lld %lld - o start: %lld %lld %lld\n",
	  ThisGrid, ThisTask, nz_complex,
	  box_real.upper[0], box_real.upper[1], box_real.upper[2],
	  box_real.lower[0], box_real.lower[1], box_real.lower[2],
	  box_complex.upper[0], box_complex.upper[1], box_complex.upper[2],
	  box_complex.lower[0], box_complex.lower[1], box_complex.lower[2]);
  
  
  GRID.total_local_size_fft = 2 * nz_complex;
  GRID.total_local_size     = nx * ny * nz_complex;
  
  MyGrids[0].off = 0;
  // order the sub-blocks by row-major order (i, j, k), k first, then j, then i

  /* int           i; */
  /* intint index; */

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
  cuda_init();
  //ptrdiff_t DIM[3];
  long long int DIM[3];
  //int       pfft_flags;
  int       ThisGrid;

  cufftResult_t status; //Check whether cuda functions are returning
    
  dprintf(VMSG, 0, "[%s] Computing fft plans\n",fdate());

  for (ThisGrid = 0; ThisGrid < Ngrids; ThisGrid++)
    {

      //cufftHandle GRID.forward_plan;
      //cufftHandle GRID.reverse_plan;
  
      cudaStreamCreate(&stream);
  
      status = cufftCreate(&GRID.forward_plan);
      if (status != CUFFT_SUCCESS) {printf("!!! forward plan cufftCreate ERROR %d !!!\n", status);}
  
      status = cufftCreate(&GRID.reverse_plan);
      if (status != CUFFT_SUCCESS) {printf("!!! reverse plan cufftCreate ERROR %d !!!\n", status);}
  
      status = cufftMpAttachComm(GRID.forward_plan, CUFFT_COMM_MPI, &comm);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftMpAttachComm ERROR %d !!!\n", status);}

      status = cufftSetStream(GRID.forward_plan, stream);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftSetStream ERROR %d !!!\n", status);}

      status = cufftMpAttachComm(GRID.reverse_plan, CUFFT_COMM_MPI, &comm);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftMpAttachComm ERROR %d !!!\n", status);}

      status = cufftSetStream(GRID.reverse_plan, stream);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftSetStream ERROR %d !!!\n", status);}

      /*Initialization*/
      int64 nx               = (int64)GRID.GSglobal[_x_];
      int64 ny               = (int64)GRID.GSglobal[_y_];
      int64 nz               = (int64)GRID.GSglobal[_z_];
      int64 nz_real          = nz;
      int64 nz_complex       = (((int64)(GRID.GSglobal[_z_])/2)+1);
      int64 nz_real_padded   = 2*nz_complex;

      Box3D box_real, box_complex;
  
      {
	// Input data are X-slabs. 
	// Strides are packed and in-place (i.e., real is padded)
	int64 lower[3]   = {nx / NTasks * (ThisTask),   0,  0};
	int64 upper[3]   = {nx / NTasks * (ThisTask+1), ny, nz_real};
	int64 strides[3] = {(upper[1]-lower[1])*nz_real_padded, nz_real_padded, 1};
	box_real = make_box(lower, upper, strides);
      }
  
      {
	// Output data are Y-slabs.
	// Strides are packed
	int64 lower[3]   = {0,  ny / NTasks * (ThisTask),   0};
	int64 upper[3]   = {nx, ny / NTasks * (ThisTask+1), nz_complex};
	int64 strides[3] = {(upper[1]-lower[1])*(upper[2]-lower[2]), (upper[2]-lower[2]), 1};
	box_complex = make_box(lower, upper, strides);
      }
      /*End of initialization*/



      
      //Distribute r2c plan
      status = cufftXtSetDistribution(GRID.forward_plan, 3, box_real.lower, box_real.upper, box_complex.lower,
				      box_complex.upper, box_real.strides, box_complex.strides);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftXtSetDistribution r2c ERROR %d !!!\n", status);}

      //Distribute c2r plan
      status = cufftXtSetDistribution(GRID.reverse_plan, 3, box_real.lower, box_real.upper, box_complex.lower,
				      box_complex.upper, box_real.strides, box_complex.strides);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftXtSetDistribution c2r ERROR %d !!!\n", status);}

      
      //Allocate the descriptors
      /*
      status = cufftXtMalloc(GRID.forward_plan, &forward_desc, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMalloc forward ERROR %d !!!\n", status);}

      status = cufftXtMalloc(GRID.reverse_plan, &reverse_desc, CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMalloc reverse ERROR %d !!!\n", status);}
      */
      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];

      // create plan for out-of-place DFT 

      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE ;
      /*
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_OUT;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
      pfft_plan_with_nthreads(internal.nthreads_fft);
#endif
      //GRID.forward_plan = pfft_plan_dft_r2c_3d(DIM, rvector_fft[ThisGrid], cvector_fft[ThisGrid],
      //				       FFT_Comm, PFFT_FORWARD, pfft_flags);
      */
      
      //Make the forward FFT plan
      
      size_t workspace_f;
      status = cufftMakePlan3d(GRID.forward_plan, DIM[_x_], DIM[_y_], DIM[_z_], CUFFT_D2Z, &workspace_f);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftMakePlan3d ERROR %d !!!\n", status);}

      
      DIM[_x_] = GRID.GSglobal[_x_];
      DIM[_y_] = GRID.GSglobal[_y_];
      DIM[_z_] = GRID.GSglobal[_z_];
      
      //pfft_flags = PFFT_MEASURE | PFFT_TUNE;
      //pfft_flags = PFFT_MEASURE;
      /*
      pfft_flags = 0;
      if(params.use_transposed_fft)
	pfft_flags |= PFFT_TRANSPOSED_IN;
#ifdef USE_FFT_THREADS
      fftw_plan_with_nthreads(internal.nthreads_fft);
      pfft_plan_with_nthreads(internal.nthreads_fft);
#endif
      GRID.reverse_plan = pfft_plan_dft_c2r_3d(DIM, cvector_fft[ThisGrid], rvector_fft[ThisGrid], 
					       FFT_Comm, PFFT_BACKWARD, pfft_flags);
      */

      //Make the reverse FFT plan
      
      size_t workspace_r;
      status = cufftMakePlan3d(GRID.reverse_plan, DIM[_x_], DIM[_y_], DIM[_z_], CUFFT_Z2D, &workspace_r);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftMakePlan3d ERROR %d !!!\n", status);}
      
    }
      
  dprintf(VMSG, 0, "[%s] fft plans done\n",fdate());

  return 0;
}


double forward_transform(int ThisGrid)
{
  cuda_init();
  double time;
  cufftDoubleComplex *complex_grid;
  cufftDoubleReal *real_grid;
  cudaError_t mmm;
  cufftResult_t status;


  //Define the descriptors
  cudaLibXtDesc *forward_desc;
  cudaLibXtDesc *reverse_desc;
  
  // Alloco fftwgrid su GPU utilizzando cudaMalloc
  /*
  long long int DIM[3];
  DIM[_x_] = GRID.GSglobal[_x_];
  DIM[_y_] = GRID.GSglobal[_y_];
  DIM[_z_] = GRID.GSglobal[_z_];

  size_t workspace_f;
  status = cufftMakePlan3d(GRID.forward_plan, DIM[_x_], DIM[_y_], DIM[_z_], CUFFT_D2Z, &workspace_f);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftMakePlan3d ERROR %d !!!\n", status);}
  */
  long long unsigned size_finta_fft = (long long unsigned)MyGrids[0].total_local_size;

  mmm=cudaMalloc(&real_grid, (size_t)(size_finta_fft*sizeof(cufftDoubleReal)));
  if (mmm != cudaSuccess) {printf("!!! cudaMalloc real ERROR %d !!!\n", mmm);}
  
  mmm=cudaMalloc(&complex_grid, (size_t)(size_finta_fft*sizeof(cufftDoubleComplex)));
  if (mmm != cudaSuccess) {printf("!!! cudaMalloc complex ERROR %d !!!\n", mmm);}

  //Allocate the descriptors
     
  status = cufftXtMalloc(GRID.forward_plan, &forward_desc, CUFFT_XT_FORMAT_DISTRIBUTED_INPUT);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMalloc forward ERROR %d !!!\n", status);}

  status = cufftXtMalloc(GRID.reverse_plan, &reverse_desc, CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMalloc reverse ERROR %d !!!\n", status);}
     
       
  time=MPI_Wtime();

  // Copy data from the CPU to the GPU.
  // The CPU data is distributed according to CUFFT_XT_FORMAT_DISTRIBUTED_INPUT
  status = cufftXtMemcpy(GRID.forward_plan, forward_desc, real_grid, CUFFT_COPY_HOST_TO_DEVICE);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMemcpyH2D ERROR %d !!!\n", status);}
  
  status = cufftXtExecDescriptor(GRID.forward_plan, forward_desc, forward_desc, CUFFT_FORWARD);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftExecDescriptor ERROR %d !!!\n", status);}
  
  // Copy memory back to the CPU. Data is now distributed according to CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT
  status = cufftXtMemcpy(GRID.forward_plan, complex_grid, forward_desc, CUFFT_COPY_DEVICE_TO_HOST);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMemcpyD2H ERROR %d !!!\n", status);}
  
  return MPI_Wtime()-time;
}


double reverse_transform(int ThisGrid)
{
  cuda_init();
  int i;
  double time;
  cufftResult_t status;
  /*
  long long int DIM[3];
  DIM[_x_] = GRID.GSglobal[_x_];
  DIM[_y_] = GRID.GSglobal[_y_];
  DIM[_z_] = GRID.GSglobal[_z_];

  size_t workspace_f;
  status = cufftMakePlan3d(GRID.reverse_plan, DIM[_x_], DIM[_y_], DIM[_z_], CUFFT_D2Z, &workspace_f);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftMakePlan3d ERROR %d !!!\n", status);}
  */
  
  time=MPI_Wtime();

  // Copy data from the CPU to the GPU.
  // The CPU data is distributed according to CUFFT_XT_FORMAT_DISTRIBUTED_INPUT
  status = cufftXtMemcpy(GRID.reverse_plan, reverse_desc, complex_grid, CUFFT_COPY_HOST_TO_DEVICE);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMemcpyH2D ERROR %d !!!\n", status);}
  
  status = cufftXtExecDescriptor(GRID.reverse_plan, reverse_desc, reverse_desc, CUFFT_INVERSE);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftExecDescriptor ERROR %d !!!\n", status);}
  
  // Copy memory back to the CPU. Data is now distributed according to CUFFT_XT_FORMAT_DISTRIBUTED_OUTPUT
  status = cufftXtMemcpy(GRID.reverse_plan, real_grid, reverse_desc, CUFFT_COPY_DEVICE_TO_HOST);
  if (status != CUFFT_SUCCESS) {printf("!!! cufftXtMemcpyD2H ERROR %d !!!\n", status);}
  
  dvec         NORM   = {GRID.norm, GRID.norm, GRID.norm, GRID.norm};
  unsigned int mysize = GRID.total_local_size_fft / 4;

  rvector_fft = (double**)real_grid;
     
#pragma GCC ivdep  
  for (i = 0; i < mysize; i++)
    *((dvec*)rvector_fft[ThisGrid] + i) *= NORM;
    
  for (i = GRID.total_local_size_fft - GRID.total_local_size_fft%4 ; i < GRID.total_local_size_fft; i++)
    rvector_fft[ThisGrid][i] *= GRID.norm;

  // non-vector code 
  /* for (i = 0 ; i < GRID.total_local_size_fft; i++) */
  /*   rvector_fft[ThisGrid][i] *= GRID.norm; */

  return MPI_Wtime() - time;
}


int finalize_fft()
{
  cuda_init();
#ifndef RECOMPUTE_DISPLACEMENTS
  int ThisGrid;

  cufftResult_t status; //Check whether cuda functions are returning
  
  for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--)
    {
      cufftXtFree(forward_desc);
      cufftXtFree(reverse_desc);
      status = cufftDestroy(GRID.forward_plan);
      status = cufftDestroy(GRID.reverse_plan);
      if (status != CUFFT_SUCCESS) {printf("!!! cufftDestroy fftwgrid ERROR %d !!!\n", status);}
    }

  /* for (ThisGrid = Ngrids-1; ThisGrid >= 0; ThisGrid--) */
  /*   { */
  /*     pfft_free(GRID.forward_plan); */
  /*     pfft_free(GRID.reverse_plan); */
  /*   } */

  cudaStreamDestroy(stream);
  //pfft_cleanup();
  //MPI_Comm_free(&FFT_Comm);
#endif

  return 0;
}
#endif //close initial if defined(CUFFTMP)
