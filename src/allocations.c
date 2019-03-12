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

//#define VERBOSE

static size_t cubes_ordering_memory, group_memory, frag_prods_memory, frag_arrays_memory, frag_memory;
static size_t prods_memory, fft_memory, all_memory, first_allocated_memory;
static size_t fmax_memory, fields_memory;

#define ALIGN_MEMORY_BLOCK(M) while( (M) % ALIGN ) M++
//#define ALIGN_MEMORY_BLOCK(M) 


int set_main_memory(unsigned int TotalNP_pertask)
{

  /* here it computes the size of required memory 
     and returns the number of allocatable particles for fragmentation */

  cubes_ordering_memory = NTasks * sizeof(int);
  prods_memory  = MyGrids[0].total_local_size * sizeof(product_data);       // products
  ALIGN_MEMORY_BLOCK( prods_memory );
  
  fields_memory = 0;
  
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    {
      fields_memory += MyGrids[igrid].total_local_size_fft * sizeof(double); // kdensity
      ALIGN_MEMORY_BLOCK( fields_memory );
    }
  
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    {
      fields_memory += 6 * MyGrids[igrid].total_local_size * sizeof(double); // second_derivatives
      ALIGN_MEMORY_BLOCK( fields_memory );
    }
  
#ifdef TWO_LPT
  fields_memory += MyGrids[0].total_local_size_fft * sizeof(double);       // kvector_2LPT
  ALIGN_MEMORY_BLOCK( fields_memory );
#ifdef THREE_LPT
  fields_memory += MyGrids[0].total_local_size_fft * sizeof(double);   // kvector_3LPT_*
  ALIGN_MEMORY_BLOCK( fields_memory );
  fields_memory += MyGrids[0].total_local_size_fft * sizeof(double);   // kvector_3LPT_*
  ALIGN_MEMORY_BLOCK( fields_memory );
#endif
#endif
  
  for ( int igrid =0, fft_memory = 0; igrid < Ngrids; igrid++ )
    {
      fft_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( fft_memory );
      fft_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( fft_memory );
    }
  
  first_allocated_memory = prods_memory + fields_memory;
  fmax_memory            = first_allocated_memory + fft_memory;

  /* this is the memory needed to fragment the collapsed medium, including PLC data;
     we assume that all peaks are at worst 1/10 of the particles in the well resolved region 
     we add to it the memory for the fragmentation map */
  group_memory = subbox.PredNpeaks * (sizeof(group_data) + sizeof(histories_data)) +
    2*subbox.maplength*sizeof(unsigned int);
#ifdef PLC
  group_memory += plc.Nmax * sizeof(plcgroup_data);
#endif
  ALIGN_MEMORY_BLOCK( group_memory );
  
  /* checks that the requested memory is sufficient for fmax */
  /* the +2 is put to have a good margin for rounding of integer divisions etc. */
  if (params.MaxMemPerParticle < (prods_memory + fields_memory + fft_memory)/TotalNP_pertask + 2)
    {
      if (!ThisTask)
	{
	  printf("ERROR: the memory specified in MaxMemPerParticle is insufficient to perform the run,\n");
	  printf("       please raise it to at least %d\n", (int)(prods_memory + fields_memory + fft_memory)/TotalNP_pertask + 2);
	}
      return -1;
    }
  else
    /* -1 is to compensate for remainder in integer division */
    return (TotalNP_pertask*params.MaxMemPerParticle - prods_memory - group_memory)/
      (sizeof(product_data) + 6 * sizeof(int)) -1;

}




int allocate_main_memory(void)
{
  /* 
     This routine allocates all the memory needed to contain information
     on particles and groups.

     ZELDOVICH DISPLACEMENTS:

     name                sizeof                    N            

     products            4 double + 1 int = 36     MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size

     first_derivatives and density will share the same memory as second_derivatives

     these will be allocated by fftw routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     Total for ZEL: 92 + 16 = 108 + overhead for fft matrices

     2LPT DISPLACEMENTS

     name                sizeof                    N            

     products            4 double + 1 int = 36     MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size
     kvector_2LPT        1 double = 8              MyGrids[0].total_local_size_fft

     first_derivatives and density will share the same memory as second_derivatives
     source_2LPT will share the same memory as kvector_2LPT

     this is needed by GenIC:

     ++ NOTE: seedtable allocation/deallocation has been moved in GenIC
     seedtable           1 int = 4                 MyGrids[*].GSglobal_x * MyGrids[*].GSglobal_y

     these will be allocated by fftw routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     Total for 2LPT: 100 + 16 = 116 + overhead for fft matrices

     Products for fragmentation will require the same memory as those for fmax
     plus the overhead for boundary layers.  

     Group catalogs will require less memory than the fmax products. This will be 
     checked with an estimate of the number of halos.

     prods_memory:            memory needed by fmax products
     first_allocated_memory:  memory needed by fmax, exluding fftw vectors
     fft_memory:              memory needed by fftw vectors
     fmax_memory:             memory needed by fmax, including fftw vectors 
     frag_prods_memory:       memory needed by fragment products
     group_memory:            memory needed by group catalogs
     frag_memory:             memory needed by fragment
     all_memory:              total memory needed

   */

  size_t memory;
  struct map{
    int lz;
    double on, am, m1, m2, m3, m4, m5, m6, m7, m8, m9;
  } mymap;
  MPI_Status status;
#ifdef VERBOSE
  size_t *bcast;
  int     bsize, n;
  void   *last;
#endif

  /* Fmax: prods_memory + fields_memory + fft_memory
     fragment: prods_memory + frag_prods + frag_arrays + group_memory */
  
  frag_prods_memory  = subbox.Nalloc * sizeof(product_data);
  frag_arrays_memory = subbox.Nalloc * 6 * sizeof(int);

  frag_memory         = prods_memory + frag_prods_memory +
    frag_arrays_memory + group_memory;

  /* this is the largest amount of memory needed by the code (they are almost equal) */
  all_memory = (fmax_memory > frag_memory? fmax_memory : frag_memory);

  /* map of memory usage for all tasks */

  double overN = (double)(MyGrids[0].total_local_size > subbox.Ngood ? MyGrids[0].total_local_size : subbox.Ngood);
  mymap.lz = MyGrids[0].GSlocal[_z_];
  mymap.am = (double)all_memory/MBYTE;
  mymap.on = overN;
  mymap.m1 = (double)prods_memory/overN;
  mymap.m2 = (double)fields_memory/overN;
  mymap.m3 = (double)fft_memory/overN;
  mymap.m4 = (double)fmax_memory/overN;
  mymap.m5 = (double)frag_prods_memory/overN;
  mymap.m6 = (double)frag_arrays_memory/overN;
  mymap.m7 = (double)group_memory/overN;
  mymap.m8 = (double)frag_memory/overN;
  mymap.m9 = (double)all_memory/overN;

  if (!ThisTask)
    {
      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "Map of memory usage for all MPI tasks\n");
      dprintf(VMSG, 0, "Task N. planes    mem(MB) overhead   products   fields     ffts     fmax  frag pr.  groups fragment  total bytes per particle\n");
    }
  
  for ( int i = 0; i < NTasks; i++ )
    {
      if (!ThisTask)
	{
	  if (i)
	    MPI_Recv((void*)&mymap, sizeof(struct map), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
	  dprintf(VMSG, 0, "%6d  %6d  %8.0f  %6.1f       "
		  "%6.1f   %6.1f   %6.1f   %6.1f   "
		  "%6.1f   %6.1f   %6.1f   %6.1f %6.1f\n",
		  i, mymap.lz, mymap.am, mymap.on,
		  mymap.m1, mymap.m2, mymap.m3,
		  mymap.m4, mymap.m5, mymap.m6,
		  mymap.m7, mymap.m8, mymap.m9 );
	}
      
      else if (ThisTask==i)
	MPI_Send((void*)&mymap, sizeof(struct map), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
  
  if (!ThisTask)
    {
      dprintf(VMSG, 0, "\n");  
      fflush(stdout);
    }

  // it tests that the required memory is lower than the specified limit
  if ((double)all_memory/MBYTE > params.MaxMem)
    {
      dprintf(VXERR, ThisTask, "ERROR on task %d: the run requires more memory than the MaxMem value\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  // it tests than there is enough space to allocate all needed memory
  //main_memory=(char*)aligned_alloc(32, (all_memory + cubes_ordering_memory) * sizeof(char));
  main_memory = (char*)calloc(all_memory + cubes_ordering_memory,  sizeof(char));
  if (main_memory==0x0)
    {
      dprintf(VXERR, ThisTask, "ERROR on task %d: I cannot allocate memory\n",ThisTask);
      fflush(stdout);
      return 1;
    }
  free(main_memory);

  // allocates main memory
  cubes_ordering = (unsigned int*)calloc(NTasks, sizeof(unsigned int));
  main_memory    = (char*)aligned_alloc(32, first_allocated_memory * sizeof(char));
  //main_memory=(char*)calloc(first_allocated_memory, sizeof(char));
  if (main_memory==0x0)
    {
      dprintf(VXERR, ThisTask, "ERROR on taks %d: I cannot allocate memory after first successful attempt\n", ThisTask);
      fflush(stdout);
      return 1;
    }

  // sets pointers to the allocated memory
  products = (product_data *)main_memory;
  
  memory   = prods_memory;                        // this has already been aligned at 32bytes
  
  for (int igrid = 0; igrid < Ngrids; igrid++)
    {
      kdensity[igrid] = (double *)(main_memory + memory);
      memory         += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( memory );
    }
  
  for (int igrid = 0; igrid < Ngrids; igrid++)
    {
      for ( int i = 0; i < 6; i++ )
	{
	  second_derivatives[igrid][i] = (double *)(main_memory + memory);
	  memory                      += MyGrids[igrid].total_local_size * sizeof(double);
	  ALIGN_MEMORY_BLOCK( memory );
	}
      
      for ( int i = 0; i < 3; i++ )
	first_derivatives[igrid][i] = second_derivatives[igrid][i];
      density[igrid] = second_derivatives[igrid][0];
    }
  
#ifdef TWO_LPT
  kvector_2LPT   = (double *)(main_memory + memory);
  source_2LPT    = kvector_2LPT;
  memory        += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( memory );
  
#ifdef THREE_LPT
  kvector_3LPT_1 = (double *)(main_memory + memory);
  memory        += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( memory );
  
  kvector_3LPT_2 = (double *)(main_memory + memory);
  memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( memory );
  
  source_3LPT_1  = kvector_3LPT_1;
  source_3LPT_2  = kvector_3LPT_2;
#endif
#endif

  /* allocates fftw vectors */
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    {
      rvector_fft[igrid] = pfft_alloc_real(MyGrids[igrid].total_local_size_fft);
      cvector_fft[igrid] = pfft_alloc_complex(MyGrids[igrid].total_local_size_fft/2);

      //rvector_fft[igrid] = (double*)      aligned_alloc(32, (MyGrids[igrid].total_local_size_fft      )* sizeof(double));
      //cvector_fft[igrid] = (pfft_complex*)aligned_alloc(32, (MyGrids[igrid].total_local_size_fft/2 + 1)*sizeof(pfft_complex));

      if ( rvector_fft[igrid] == 0x0 || cvector_fft[igrid] == 0x0 )
	{
	  dprintf(VXERR, ThisTask, "ERROR on taks %d: I cannot allocate memory for fft vectors\n", ThisTask);
	  fflush(stdout);
	  return 1;
	}
    }  

  if (!ThisTask)
    dprintf(VXX, ThisTask, "[%s] Memory has been successfully allocated\n",fdate());

#ifdef VERBOSE

  last = main_memory + first_allocated_memory;
  
#ifndef TWO_LPT
  bsize  = 3+3*Ngrids;
#else
  bsize  = 4+3*Ngrids;
#endif
  bsize *= 2;
  bcast  = (size_t*)malloc(bsize*sizeof(size_t));

  n=0;
  bcast[n++] = (size_t)((void*)kdensity[0] - (void*)products);

  for ( int igrid = 1; igrid < Ngrids; igrid++ )
    bcast[n++] = (size_t)((void*)kdensity[igrid] - (void*)kdensity[igrid-1]);

  bcast[n++] = (size_t)((void*)second_derivatives[0][0] - (void*)kdensity[Ngrids-1]);
  for ( int igrid = 1; igrid < Ngrids; igrid++)
    bcast[n++] = (size_t)((void*)second_derivatives[igrid][0] - (void*)second_derivatives[igrid-1][0]);

#ifdef TWO_LPT
  bcast[n++] = (size_t)((void*)kvector_2LPT - (void*)second_derivatives[Ngrids-1][0]);
#endif

  bcast[n++] = (size_t)(last - (void*)products);
  bcast[n++] = fftmemory;

  bcast[n++] = prods_memory;
  for ( int igrid=0; igrid<Ngrids; igrid++ )
    bcast[n++] = MyGrids[igrid].total_local_size_fft * sizeof(double);
  
  for ( int igrid=0; igrid<Ngrids; igrid++)
    bcast[n++] = 6*MyGrids[igrid].total_local_size * sizeof(double);
  
#ifdef TWO_LPT
#ifndef THREE_LPT
  bcast[n++] = MyGrids[0].total_local_size_fft * sizeof(double);
#else
  bcast[n++] = 3*MyGrids[0].total_local_size_fft * sizeof(double);
#endif
#endif
  for ( int igrid = 0; igrid < Ngrids; igrid++ )
    bcast[n++] = MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);
  
  bcast[n++] = memory;
  bcast[n++] = allmemory;


  if (!ThisTask)
    {
      printf("Map of memory allocations for fmax (in byte):\n");
      printf("        TASK ");
      printf("    PRODUCTS");
      for (igrid=0; igrid<Ngrids; igrid++)
	printf("    KDENS G.%d",igrid);
      for (igrid=0; igrid<Ngrids; igrid++)
	printf("    DERIV G.%d",igrid);
#ifdef TWO_LPT
      printf("   KVEC 2LPT ");
#ifdef THREE_LPT
      printf("   KVEC 3LPT ");
#endif
#endif
      for (igrid=0; igrid<Ngrids; igrid++)
	printf("    SEEDS G.%d ",igrid);
      printf(" MAIN MEMORY ");
      printf("- FFT VECTORS / REQ. MEMORY\n");
    }

  for ( int i = 0; i < NTasks; i++)
    {

      if (!ThisTask)
	{
	  if (i)
	    s=MPI_Recv((void*)bcast, bsize*sizeof(size_t), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
	  printf(" FOUND %5d",i);
	  for (n=0; n<bsize/2; n++)
	    printf(" %12zu",bcast[n]);
	  printf("\n");
	  printf(" PRED  %5d",i);
	  for (n=bsize/2; n<bsize; n++)
	    printf(" %12zu",bcast[n]);
	  printf("\n");
	}
      else if (ThisTask==i)
	s=MPI_Send((void*)bcast, bsize*sizeof(size_t), MPI_BYTE, 0, 0, MPI_COMM_WORLD);

    }

  if (!ThisTask)
    {
      printf("\n");  
      fflush(stdout);
    }

  free(bcast);

#endif  // END of VERBOSE

  return 0;
}


int deallocate_fft_vectors(int ThisGrid)
{

  pfft_free(cvector_fft[ThisGrid]);
  pfft_free(rvector_fft[ThisGrid]);

  return 0;
}

int reallocate_memory_for_fragmentation()
{

  size_t memory;
  void *tmp;

  kdensity=0x0;
  density=0x0;
  first_derivatives=0x0;
  second_derivatives=0x0;
#ifdef TWO_LPT
  kvector_2LPT=0x0;
  source_2LPT=0x0;
#ifdef THREE_LPT
  kvector_3LPT_1=0x0;
  source_3LPT_1=0x0;
  kvector_3LPT_2=0x0;
  source_3LPT_2=0x0;
#endif
#endif

  /* products in the sub-box are stored in the memory that has been freed */

  /* enlarge memory allocation if needed */
  if (frag_memory > first_allocated_memory)
    {
      tmp=(void*)realloc(main_memory,frag_memory);
      if (tmp!=0x0)
	{
	  main_memory=tmp;
	  products = (product_data *)main_memory;
	  if (!ThisTask)
	    printf("Task 0 reallocated memory for %f Gb\n",(double)frag_memory/GBYTE);
	}
      else
	{
	  printf("ERROR on task %d: could not reallocate memory from %zu to %zu bytes\n",ThisTask,
		 first_allocated_memory, frag_memory);
	  fflush(stdout);
	  return 1;
	}
    }

  /* the frag structure is located just after memory products */
  frag = (product_data *)(main_memory + prods_memory);
  memory=prods_memory+frag_prods_memory;



  frag_pos = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  indices = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  indicesY = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  sorted_pos = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  group_ID = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  linking_list = (int*)(main_memory+memory);
  memory+=subbox.Nalloc*sizeof(int);
  groups = (group_data *)(main_memory + memory);
  memory += subbox.PredNpeaks*sizeof(group_data);
  wheretoplace_mycat = (void*)(main_memory + memory);
  memory += subbox.PredNpeaks*sizeof(histories_data);
  frag_map = (unsigned int*)(main_memory+memory);
  memory+=subbox.maplength*sizeof(unsigned int);
  frag_map_update = (unsigned int*)(main_memory+memory);
  memory+=subbox.maplength*sizeof(unsigned int);
#ifdef PLC
  plcgroups = (plcgroup_data *)(main_memory + memory);
  memory+=plc.Nmax*sizeof(plcgroup_data);
#endif


  return 0;
}
