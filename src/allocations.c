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

static size_t cubes_ordering_memory, group_memory, frag_prods_memory, frag_memory, prods_memory, fft_memory,
  all_memory, first_allocated_memory, fmax_memory, fields_memory;

#define ALIGN_MEMORY_BLOCK(M) while( (M) % ALIGN ) M++
//#define ALIGN_MEMORY_BLOCK(M) 

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

     ++ NOTE: seedtable allocation/deallocation is being moved in GenIC
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

     IN CASE DO_NOT_REALLOC IS DEFINED:
     first_allocated_memory is large enough to include frag_memory, fftw vectors are an extra overhead

   */

  int igrid,i;
  size_t memory;
  struct map{
    int lz;
    double oh,am,m1,m2,m3,m4,m5,m6,m7,m8;
  } mymap;
  MPI_Status status;
#ifdef VERBOSE
  size_t *bcast;
  int bsize, n;
  void *last;
#endif

  /* Computes total required memory */
  cubes_ordering_memory = NTasks * sizeof(int);
  prods_memory = MyGrids[0].total_local_size * sizeof(product_data);   // products
  ALIGN_MEMORY_BLOCK( prods_memory );
  
  fields_memory=0;
  
  for (igrid = 0; igrid < Ngrids; igrid++)
    {
      fields_memory += MyGrids[igrid].total_local_size_fft * sizeof(double); // kdensity
      ALIGN_MEMORY_BLOCK( fields_memory );
    }
	 
  for (igrid = 0; igrid < Ngrids; igrid++)
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
  fields_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( fields_memory );
#endif
#endif

  first_allocated_memory = prods_memory + fields_memory;

  // adds memory for fft vectors
  for (igrid = 0, fft_memory = 0; igrid < Ngrids; igrid++)
    {
      fft_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( fft_memory );
      fft_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( fft_memory );
    }

  fmax_memory = first_allocated_memory + fft_memory;

  // this is the memory needed to fragment the collapsed medium, including PLC data;
  //   we assume that all peaks are at worst 1/10 of the particles
  frag_prods_memory = subbox.Npart * sizeof(product_data);
  group_memory      = subbox.Npart * 3 * sizeof(int) +
    subbox.Npart/10 * (sizeof(group_data) + sizeof(histories_data));
#ifdef PLC
  group_memory += plc.Nmax * sizeof(plcgroup_data);
#endif

#ifdef TEST_ONLY
  frag_prods_memory = 0;
  group_memory = 0;
#endif


  if (NSlices>1)
    frag_memory = prods_memory + frag_prods_memory + group_memory;
  else
    frag_memory = (group_memory > prods_memory ? group_memory : prods_memory) + frag_prods_memory;

#ifdef DO_NOT_REALLOC
  
  if (frag_memory > first_allocated_memory)
    first_allocated_memory = frag_memory;
  all_memory = first_allocated_memory + fft_memory;
  
#else
  
  // this is the largest amount of memory needed by the code
  all_memory = (fmax_memory > frag_memory? fmax_memory : frag_memory);
  
#endif

  // map of memory usage for all tasks
  mymap.lz = MyGrids[0].GSlocal[_z_];
  mymap.am=(double)all_memory/MBYTE;

  if (MyGrids[0].total_local_size)
    {
      double factor = 1.0 / (double)MyGrids[0].total_local_size;
      
      mymap.oh = (double) subbox.Npart     * factor;
      mymap.m1 = (double) prods_memory     * factor;
      mymap.m2 = (double) fields_memory    * factor;
      mymap.m3 = (double) fft_memory       * factor;
      mymap.m4 = (double) fmax_memory      * factor;
      mymap.m5 = (double) frag_prods_memory* factor;
      mymap.m6 = (double) group_memory     * factor;
      mymap.m7 = (double) frag_memory      * factor;
      mymap.m8 = (double) all_memory       * factor;
    }
  else
    {
      mymap.oh = 0.0;
      mymap.m1 = 0.0;
      mymap.m2 = 0.0;
      mymap.m3 = 0.0;
      mymap.m4 = 0.0;
      mymap.m5 = 0.0;
      mymap.m6 = 0.0;
      mymap.m7 = 0.0;
      mymap.m8 = 0.0;
    }

  if (!ThisTask)
    {
      dprintf(VMSG, 0, "\n");
      dprintf(VMSG, 0, "Map of memory usage for all MPI tasks\n");
#ifdef DO_NOT_REALLOC
      dprintf(VMSG, 0, "  NB: here we avoid reallocations\n");
#endif
      dprintf(VMSG, 0, "Task N. planes    mem(MB) overhead   products   fields     ffts     fmax  frag pr.  groups fragment  total bytes per particle\n");
    }
  
  for (i=0; i<NTasks; i++)
    {
      if (!ThisTask)
	{
	  if (i)
	    MPI_Recv((void*)&mymap, sizeof(struct map), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
	  dprintf(VMSG, 0, "%6d  %6d  %8.0f  %6.1f       "
		  "%6.1f   %6.1f   %6.1f   %6.1f   "
		  "%6.1f   %6.1f   %6.1f   %6.1f\n",
		  i, mymap.lz, mymap.am, mymap.oh,
		  mymap.m1, mymap.m2, mymap.m3, mymap.m4,
		  mymap.m5, mymap.m6, mymap.m7, mymap.m8);
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
  main_memory=(char*)calloc(all_memory + cubes_ordering_memory,  sizeof(char));
  if (main_memory==0x0)
    {
      dprintf(VXERR, ThisTask, "ERROR on task %d: I cannot allocate memory\n",ThisTask);
      fflush(stdout);
      return 1;
    }
  free(main_memory);
  main_memory = 0x0;

  // allocates main memory
  cubes_ordering = (unsigned int*)calloc(NTasks, sizeof(unsigned int));
  main_memory    = (char*)aligned_alloc(ALIGN, first_allocated_memory * sizeof(char));
  //main_memory=(char*)calloc(first_allocated_memory, sizeof(char));
  if (main_memory == 0x0)
    {
      dprintf(VXERR, ThisTask, "ERROR on taks %d: I cannot allocate memory after first successful attempt\n", ThisTask);
      fflush(stdout);
      return 1;
    }
  *main_memory = 0;
  // sets pointers to the allocated memory
  products = (product_data *)main_memory;
  
  memory   = prods_memory;                        // this has already been aligned at 32bytes
  
  for (igrid = 0; igrid < Ngrids; igrid++)
    {
      kdensity[igrid] = (double *)(main_memory + memory);
      memory         += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( memory );
    }
  
  for (igrid = 0; igrid < Ngrids; igrid++)
    {
      for (i = 0; i < 6; i++)
	{
	  second_derivatives[igrid][i] = (double *)(main_memory + memory);
	  memory                      += MyGrids[igrid].total_local_size * sizeof(double);
	  ALIGN_MEMORY_BLOCK( memory );
	}
      
      for (i = 0; i < 3; i++)
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
  for (igrid = 0; igrid < Ngrids; igrid++)
    {
      rvector_fft[igrid] = pfft_alloc_real(MyGrids[igrid].total_local_size_fft);
      cvector_fft[igrid] = pfft_alloc_complex(MyGrids[igrid].total_local_size_fft/2);

      //rvector_fft[igrid] = (double*)      aligned_alloc(32, (MyGrids[igrid].total_local_size_fft      )* sizeof(double));
      //cvector_fft[igrid] = (pfft_complex*)aligned_alloc(32, (MyGrids[igrid].total_local_size_fft/2 + 1)*sizeof(pfft_complex));

      if (rvector_fft[igrid] == 0x0 || cvector_fft[igrid] == 0x0)
	{
	  dprintf(VXERR, ThisTask, "ERROR on taks %d: I cannot allocate memory for fft vectors\n", ThisTask);
	  fflush(stdout);
	  return 1;
	}
    }  

  if (!ThisTask)
    dprintf(VXX, ThisTask, "[%s] Memory has been successfully allocated\n",fdate());

#ifdef VERBOSE

  last = main_memory+first_allocated_memory;
  
#ifndef TWO_LPT
  bsize  = 3+3*Ngrids;
#else
  bsize  = 4+3*Ngrids;
#endif
  bsize *= 2;
  bcast  = (size_t*)malloc(bsize*sizeof(size_t));

  n=0;
  bcast[n++] = (size_t)((void*)kdensity[0]-(void*)products);

  for (igrid=1; igrid<Ngrids; igrid++)
    bcast[n++] = (size_t)((void*)kdensity[igrid]-(void*)kdensity[igrid-1]);

  bcast[n++] = (size_t)((void*)second_derivatives[0][0]-(void*)kdensity[Ngrids-1]);
  for (igrid=1; igrid<Ngrids; igrid++)
    bcast[n++] = (size_t)((void*)second_derivatives[igrid][0]-(void*)second_derivatives[igrid-1][0]);

#ifdef TWO_LPT
  bcast[n++] = (size_t)((void*)kvector_2LPT-(void*)second_derivatives[Ngrids-1][0]);
#endif

  bcast[n++] = (size_t)(last-(void*)products);
  bcast[n++] = fftmemory;

  bcast[n++] = prods_memory;
  for (igrid=0; igrid<Ngrids; igrid++)
    bcast[n++] = MyGrids[igrid].total_local_size_fft * sizeof(double);
  
  for (igrid=0; igrid<Ngrids; igrid++)
    bcast[n++] = 6*MyGrids[igrid].total_local_size * sizeof(double);
  
#ifdef TWO_LPT
#ifndef THREE_LPT
  bcast[n++] = MyGrids[0].total_local_size_fft * sizeof(double);
#else
  bcast[n++] = 3*MyGrids[0].total_local_size_fft * sizeof(double);
#endif
#endif
  for (igrid=0; igrid<Ngrids; igrid++)
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

  for (i=0; i<NTasks; i++)
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

int reallocate_memory_for_fragmentation_1()
{

  int PredNpeaks;
  size_t memory;
  void *tmp;
  
#ifdef VERBOSE
  size_t *bcast;
  int bsize, s, n;
  MPI_Status status;
  void *last;
#endif

  kdensity           =0x0;
  density            =0x0;
  first_derivatives  =0x0;
  second_derivatives =0x0;
  // seedtable allocation/deallocation is being moved in GenIC
  /* seedtable          =0x0; */
  
#ifdef TWO_LPT
  kvector_2LPT       =0x0;
  source_2LPT        =0x0;
#ifdef THREE_LPT
  kvector_3LPT_1     =0x0;
  source_3LPT_1      =0x0;
  kvector_3LPT_2     =0x0;
  source_3LPT_2      =0x0;
#endif
#endif

  // products in the sub-box are stored in the memory that has been freed

#ifndef DO_NOT_REALLOC
  // enlarge memory allocation if needed
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
#endif

  // the frag structure is located just after memory products
  if (NSlices>1)
    {
      frag = (product_data *)(main_memory + prods_memory);
      memory=prods_memory+frag_prods_memory;
    }
  else
    {
      frag = (product_data *)(main_memory + (group_memory > prods_memory ? group_memory : prods_memory));
      memory=0;
    }

  /* QUESTO CREAVA PROBLEMI DI MEMORIA, BISOGNA CAPIRE PERCHE' */
  /* for (i=0; i<frag_prods_memory; i++) */
  /*   *((char*)frag+i)=0; */

  PredNpeaks=subbox.Npart/10;

  indices            = (int*)(main_memory+memory);
  memory            += subbox.Npart*sizeof(int);
  group_ID           = (int*)(main_memory+memory);
  memory            += subbox.Npart*sizeof(int);
  linking_list       = (int*)(main_memory+memory);
  memory            += subbox.Npart*sizeof(int);
  groups             = (group_data *)(main_memory + memory);
  memory            += PredNpeaks*sizeof(group_data);
  wheretoplace_mycat = (void*)(main_memory + memory);
  memory            += PredNpeaks*sizeof(histories_data);
#ifdef PLC
  plcgroups          = (plcgroup_data *)(main_memory + memory);
  memory            += plc.Nmax*sizeof(plcgroup_data);
#endif

  return 0;
}

int reallocate_memory_for_fragmentation_2(int Npeaks)
{
  int i,PredNpeaks;
#ifdef VERBOSE
  size_t *bcast;
  int bsize, s, n;
  MPI_Status status;
#endif


  PredNpeaks=subbox.Npart/10;

  // the number of peaks was supposed to be at most 1/10 of the number of particles
  //   (we need space for Npeaks+2 groups)
  if (Npeaks+2 > PredNpeaks)
    {
      printf("ERROR on task %d: surprisingly, the number of peaks %d exceeds Npart/10 (%d)\n",
	     ThisTask,Npeaks,subbox.Npart/10);
      fflush(stdout);
      return 1;
    }

  if (!NSlices)
    products=0x0;

  // initializes group memory
  for (i = 0; i < subbox.Npart; i++)
    *(indices+i) = i;

  for (i = 0; i < subbox.Npart; i++)
    *(group_ID+i) = 0;

  for (i = 0; i < subbox.Npart; i++)
    *(linking_list + i)=0;

  for (i = 0; i < PredNpeaks*sizeof(group_data); i++)
    *((char*)groups + i)=0;

  for (i = 0; i < PredNpeaks*sizeof(histories_data); i++)
    *((char*)wheretoplace_mycat + i)=0;

#ifdef PLC
  for (i = 0; i < plc.Nmax*sizeof(plcgroup_data); i++)
    *((char*)plcgroups + i)=0;
#endif

  return 0;
}


