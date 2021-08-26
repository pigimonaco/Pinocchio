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

//#define VERBOSE

#define ALIGN_MEMORY_BLOCK(M) while( (M) % ALIGN ) M++

int organize_main_memory()
{

  /* here it computes the size of required memory 
     and returns the number of allocatable particles for fragmentation

     ZELDOVICH DISPLACEMENTS:

     name                sizeof                    N            

     products            4 float + 1 int = 20      MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size

     first_derivatives and density will share the same memory as second_derivatives

     these will be allocated by pfft routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     This is needed by GenIC:
     seedtable           1 int = 4                 MyGrids[*].GSglobal[_x_] * MyGrids[*].GSglobal[_y_]

     Total for ZEL: 76 + 16 = 92 + overhead for fragmentation + smaller things (seedtable etc)

     2LPT/3LPT DISPLACEMENTS

     name                sizeof                    N            

     products (2LPT):    7 float + 1 int = 32      MyGrids[0].total_local_size
     products (3LPT)    13 float + 1 int = 56     MyGrids[0].total_local_size
     kdensity            1 double = 8              MyGrids[*].total_local_size_fft
     second_derivatives  6 double = 48             MyGrids[*].total_local_size
     kvector_2LPT        1 double = 8              MyGrids[0].total_local_size_fft
     kvector_3LPT_1      1 double = 8              MyGrids[0].total_local_size_fft
     kvector_3LPT_2      1 double = 8              MyGrids[0].total_local_size_fft

     first_derivatives and density will share the same memory as second_derivatives
     source_?LPT* will share the same memory as kvector_?LPT*

     This is needed by GenIC:
     seedtable           1 int = 4                 MyGrids[*].GSglobal[_x_] * MyGrids[*].GSglobal[_y_]

     These will be allocated by pfft routines:

     cvector_fft         1 double = 8              MyGrids[*].total_local_size_fft
     rvector_fft         1 double = 8              MyGrids[*].total_local_size_fft

     Total for 2LPT: 32 + 80 + 16 = 128 + overhead for fragmentation + smaller things
     Total for 3LPT: 56 + 80 + 16 = 152 + overhead for fragmentation + smaller things

     Timeless snapshot adds 4 bytes to the count.

     In case or recomputation of displacements: 

     products (ZEL):     7 float + 1 int = 32      MyGrids[0].total_local_size
     products (2LPT):    13 float + 1 int = 56     MyGrids[0].total_local_size
     products (3LPT):    25 float + 1 int = 104    MyGrids[0].total_local_size
     Total for ZEL:   32 + 80 + 16 = 128 + overhead for fragmentation
     Total for 2LPT:  56 + 80 + 16 = 152 + overhead for fragmentation
     Total for 3LPT: 104 + 80 + 16 = 200 + overhead for fragmentation

     Products for fragmentation will require the same memory as those for fmax
     plus the overhead for boundary layers.  

     Group catalogs will require less memory than the fmax products. This will be 
     checked with an estimate of the number of halos.

     memory.prods:            fmax products in fft space
     memory.fields_to_keep    kdensities
     memory.fields            second_derivatives, // seedtable
     memory.first_allocated:  the three above
     memory.fft:              needed by pfft vectors
     memory.fmax_total:       needed by fmax, the sum of fft and first_allocated 
     memory.frag_prods:       fmax products in subbox space
     memory.frag_arrays:      linking list and group ID
     memory.groups:           groups catalog, including histories, map and PLC
     memory.frag_allocated:   needed by fragment, allocated
     memory.frag_total:       needed by fragment, including FFTs
     memory.all_allocated:    total allocated memory (not including FFTs)
     memory.all:              total memory needed (including FFTs)

   */

  int igrid;
  unsigned int myNalloc, Nalloc;

  /***********************************************/
  /* COLLAPSE TIMES PART                         */
  /* this is the memory needed to store products */
  memory.prods = MyGrids[0].total_local_size * sizeof(product_data);       /* products */
  ALIGN_MEMORY_BLOCK( memory.prods );
  memory.fields = memory.fields_to_keep = 0;

  /* these fields are needed after Fmax if displacements are recomputed */
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      memory.fields_to_keep += MyGrids[igrid].total_local_size_fft * sizeof(double); /* kdensity */
      ALIGN_MEMORY_BLOCK( memory.fields_to_keep );
    }
#ifdef TWO_LPT
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double);    /* kvector_2LPT */
  ALIGN_MEMORY_BLOCK( memory.fields_to_keep );
#ifdef THREE_LPT
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double);    /* kvector_3LPT_1 */
  ALIGN_MEMORY_BLOCK( memory.fields_to_keep );
  memory.fields_to_keep += MyGrids[0].total_local_size_fft * sizeof(double);    /* kvector_3LPT_2 */
  ALIGN_MEMORY_BLOCK( memory.fields_to_keep );
#endif
#endif

  /* these fields are needed only by Fmax */
  for (igrid=0; igrid<Ngrids; igrid++)
    for (int dim=0; dim<6; dim++)
      {
	memory.fields += MyGrids[igrid].total_local_size * sizeof(double); /* second_derivatives */
	ALIGN_MEMORY_BLOCK( memory.fields );
      }

  /* for (igrid=0; igrid<Ngrids; igrid++) */
  /*   memory.fields += MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);  /\* seedtable *\/ */

  /* memory allocated by PFFT */
  for (igrid=0, memory.fft=0; igrid<Ngrids; igrid++)
    {
      memory.fft += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( memory.fft );
      memory.fft += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( memory.fft );
    }

  memory.first_allocated = memory.prods + memory.fields_to_keep + memory.fields;
  memory.fmax_total = memory.first_allocated + memory.fft;

  /********************************************************************************/
  /* FRAGMENTATION PART                                                           */
  /* this is the memory needed to fragment the collapsed medium, including PLC data
     and the map for needed particles on the boundary */
  memory.groups = subbox.PredNpeaks * sizeof(group_data);
  ALIGN_MEMORY_BLOCK( memory.groups );
  memory.groups += subbox.PredNpeaks *  + sizeof(histories_data);
  ALIGN_MEMORY_BLOCK( memory.groups );
#ifdef PLC
  memory.groups += plc.Nmax * sizeof(plcgroup_data);
  ALIGN_MEMORY_BLOCK( memory.groups );
#endif
  memory.groups += subbox.maplength*sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK( memory.groups );
  memory.groups += subbox.maplength*sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK( memory.groups );

  /* Here it computes the number of particles it can allocate for the subbox */
  size_t other_mem = memory.prods + memory.groups 
#ifdef RECOMPUTE_DISPLACEMENTS
    + memory.fields_to_keep + memory.fft
#endif
    ;

  if (other_mem > MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle)
    myNalloc=0;
  else
   /* -10 is to compensate for remainder in integer division */
    myNalloc = (MyGrids[0].ParticlesPerTask*params.MaxMemPerParticle - other_mem) 
      / (sizeof(product_data) + FRAGFIELDS * sizeof(int)) -10;

  /* Nalloc will be the smallest among all tasks */
  MPI_Reduce(&myNalloc, &Nalloc, 1, MPI_UNSIGNED, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Nalloc, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

  memory.frag_prods  = Nalloc * sizeof(product_data);
  memory.frag_arrays = Nalloc * FRAGFIELDS * sizeof(int);
  /* this is the memory occupied by fragmentation */
  memory.frag_allocated = memory.prods + memory.frag_prods + memory.groups + memory.frag_arrays
#ifdef RECOMPUTE_DISPLACEMENTS
    + memory.fields_to_keep
#endif
    ;
  memory.frag_total = memory.frag_allocated
#ifdef RECOMPUTE_DISPLACEMENTS
    + memory.fft
#endif
    ;


#ifdef CLASSIC_FRAGMENTATION
  if (myNalloc<subbox.Npart)
    {
      printf("ERROR: Task %d can allocate only %d subbox particles, while %d are needed;\n",
	     ThisTask,myNalloc,subbox.Npart);
      printf("       a large overhead is probably needed, please increase MaxMemPerParticle\n");
      return (size_t)0;
    }
#endif

  /* this is the largest amount of memory needed by the code */
  memory.all_allocated = (memory.first_allocated > memory.frag_allocated ? memory.first_allocated : memory.frag_allocated);
  memory.all = (memory.fmax_total > memory.frag_total ? memory.fmax_total : memory.frag_total);

  double myfraction=(double)Nalloc/(double)subbox.Ngood;
  double minfraction;
  MPI_Reduce(&myfraction, &minfraction, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      if (minfraction<0.5)
	{
	  printf("ERROR: for some tasks the number of allocatable particles is only a factor %8f\n",minfraction);
	  printf("       times the number of good particles in the subbox; this is far too small!\n");
	  printf("       Please increase MaxMemPerParticle\n");
	  fflush(stdout);
	  return (size_t)0;
	}
    }

  return Nalloc;
}



int allocate_main_memory()
{
  /* 
     This routine allocates all the memory needed to contain information
     on particles and groups.
  */

  int igrid,i;
  size_t count_memory;
  struct map{
    int lx,ly,lz;
    double oh,am,m1,m2,m3,m4,m5,m6,m7,m8;
  } mymap;
  MPI_Status status;
#ifdef VERBOSE
  size_t *bcast;
  int bsize, n;
  void *last;
#endif

  /* map of memory usage for all tasks */
  mymap.lx=MyGrids[0].GSlocal[_x_];
  mymap.ly=MyGrids[0].GSlocal[_y_];
  mymap.lz=MyGrids[0].GSlocal[_z_];
  mymap.am=(double)memory.all/MBYTE;
  mymap.oh=(double)subbox.Nalloc/(double)MyGrids[0].ParticlesPerTask;
  mymap.m1=(double)memory.prods/(double)MyGrids[0].ParticlesPerTask;
  mymap.m2=(double)(memory.fields+memory.fields_to_keep)/(double)MyGrids[0].ParticlesPerTask;
  mymap.m3=(double)memory.fft/(double)MyGrids[0].ParticlesPerTask;
  mymap.m4=(double)memory.fmax_total/(double)MyGrids[0].ParticlesPerTask;
  mymap.m5=(double)memory.frag_prods/(double)MyGrids[0].ParticlesPerTask;
  mymap.m6=(double)(memory.groups+memory.frag_arrays)/(double)MyGrids[0].ParticlesPerTask;
  mymap.m7=(double)memory.frag_total/(double)MyGrids[0].ParticlesPerTask;
  mymap.m8=(double)memory.all/(double)MyGrids[0].ParticlesPerTask;

  if (!ThisTask)
    {
      printf("\n");
      printf("Map of memory usage for all MPI tasks\n");
      printf("Task N.   FFT domain      mem(MB) overhead   products   fields     ffts     fmax  frag pr.  groups fragment  total bytes per particle\n");
    }
  for (i=0; i<NTasks; i++)
    {
      if (!ThisTask)
	{
	  if (i)
	    MPI_Recv((void*)&mymap, sizeof(struct map), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
	  printf("%6d  %4d-%4d-%4d  %8.0f  %6.1f       %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f\n",
		 i, mymap.lx,mymap.ly,mymap.lz, mymap.am, mymap.oh, mymap.m1, mymap.m2, mymap.m3, mymap.m4, mymap.m5, mymap.m6, mymap.m7, mymap.m8);
	}
      else if (ThisTask==i)
	MPI_Send((void*)&mymap, sizeof(struct map), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
  if (!ThisTask)
    {
      printf("\n");  
      fflush(stdout);
    }

  // ???
  memory.all+=10;

  /* it tests than there is enough space to allocate all needed memory */
  main_memory=(char*)calloc(memory.all,sizeof(char));
  if (main_memory==0x0)
    {
      printf("ERROR on task %d: I cannot allocate memory\n",ThisTask);
      fflush(stdout);
      return 1;
    }
  free(main_memory);

  /* allocates main memory */
  main_memory=(char*)calloc(memory.first_allocated,sizeof(char));
  if (main_memory==0x0)
    {
      printf("ERROR on taks %d: I cannot allocate memory after first successful attempt\n", ThisTask);
      fflush(stdout);
      return 1;
    }

  /* sets pointers to the allocated memory */
  products = (product_data *)main_memory;
  count_memory = MyGrids[0].total_local_size * sizeof(product_data);
  ALIGN_MEMORY_BLOCK( count_memory );
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      kdensity[igrid] = (double *)(main_memory + count_memory);
      count_memory += MyGrids[igrid].total_local_size_fft * sizeof(double);
      ALIGN_MEMORY_BLOCK( count_memory );
    }
#ifdef TWO_LPT
  kvector_2LPT = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( count_memory );
  source_2LPT = kvector_2LPT;
#ifdef THREE_LPT
  kvector_3LPT_1 = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( count_memory );
  kvector_3LPT_2 = (double *)(main_memory + count_memory);
  count_memory += MyGrids[0].total_local_size_fft * sizeof(double);
  ALIGN_MEMORY_BLOCK( count_memory );
  source_3LPT_1 = kvector_3LPT_1;
  source_3LPT_2 = kvector_3LPT_2;
#endif
#endif
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      for (i=0; i<6; i++)
	{
	  second_derivatives[igrid][i] = (double *)(main_memory + count_memory);
	  count_memory += MyGrids[igrid].total_local_size * sizeof(double);	  
	  ALIGN_MEMORY_BLOCK( count_memory );
	}
      for (i=0; i<3; i++)
	first_derivatives[igrid][i] = second_derivatives[igrid][i];
      density[igrid] = second_derivatives[igrid][0];
    }
//  for (igrid=0; igrid<Ngrids; igrid++)
//    {
//      seedtable[igrid] = (unsigned int*)(main_memory + count_memory);
//      count_memory += MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);
//    }
//
  /* allocates fft vectors */
  for (igrid=0; igrid<Ngrids; igrid++)
    {
      rvector_fft[igrid] = pfft_alloc_real(MyGrids[igrid].total_local_size_fft);
      cvector_fft[igrid] = pfft_alloc_complex(MyGrids[igrid].total_local_size_fft/2);
      if (rvector_fft[igrid] == 0x0 || cvector_fft[igrid] == 0x0)
	{
	  printf("ERROR on taks %d: I cannot allocate memory for fft vectors\n", ThisTask);
	  fflush(stdout);
	  return 1;
	}
    }  

  if (!ThisTask)
    printf("[%s] Memory has been successfully allocated\n",fdate());

#ifdef VERBOSE

  last=main_memory+memory.first_allocated;
#ifndef TWO_LPT
  bsize=3+3*Ngrids;
#else
  bsize=4+3*Ngrids;
#endif
  bsize*=2;
  bcast=(size_t*)malloc(bsize*sizeof(size_t));

  n=0;
  bcast[n++]=(size_t)((void*)kdensity[0]-(void*)products);
  for (igrid=1; igrid<Ngrids; igrid++)
    bcast[n++]=(size_t)((void*)kdensity[igrid]-(void*)kdensity[igrid-1]);
#ifdef TWO_LPT
  bcast[n++]=(size_t)((void*)kvector_2LPT-(void*)second_derivatives[Ngrids-1][0]);
  bcast[n++]=(size_t)((void*)seedtable[0]-(void*)kvector_2LPT);
#else
  bcast[n++]=(size_t)((void*)seedtable[0]-(void*)second_derivatives[Ngrids-1][0]);
#endif
  bcast[n++]=(size_t)((void*)second_derivatives[0][0]-(void*)kdensity[Ngrids-1]);
  for (igrid=1; igrid<Ngrids; igrid++)
    bcast[n++]=(size_t)((void*)second_derivatives[igrid][0]-(void*)second_derivatives[igrid-1][0]);
  /* for (igrid=1; igrid<Ngrids; igrid++) */
  /*   bcast[n++]=(size_t)((void*)seedtable[igrid]-(void*)seedtable[igrid-1]); */
  bcast[n++]=0;
  bcast[n++]=(size_t)(last-(void*)seedtable[Ngrids-1]);
  bcast[n++]=(size_t)(last-(void*)products);
  bcast[n++]=memory.fft;

  bcast[n++]=memory.prods;
  for (igrid=0; igrid<Ngrids; igrid++)
    bcast[n++]=MyGrids[igrid].total_local_size_fft * sizeof(double);
  for (igrid=0; igrid<Ngrids; igrid++)
    bcast[n++]=6*MyGrids[igrid].total_local_size * sizeof(double);
#ifdef TWO_LPT
#ifndef THREE_LPT
  bcast[n++]=MyGrids[0].total_local_size_fft * sizeof(double);
#else
  bcast[n++]=3*MyGrids[0].total_local_size_fft * sizeof(double);
#endif
#endif
  for (igrid=0; igrid<Ngrids; igrid++)
    bcast[n++]=MyGrids[igrid].GSglobal[_x_] * MyGrids[igrid].GSglobal[_y_] * sizeof(unsigned int);
  bcast[n++]=count_memory;
  bcast[n++]=memory.all;


  if (!ThisTask)
    {
      printf("Map of memory allocations for fmax (in byte):\n");
      printf("        TASK ");
      printf("    PRODUCTS");
      for (igrid=0; igrid<Ngrids; igrid++)
	printf("    KDENS G.%d",igrid);
#ifdef TWO_LPT
      printf("   KVEC 2LPT ");
#ifdef THREE_LPT
      printf("   KVEC 3LPT ");
#endif
#endif
      for (igrid=0; igrid<Ngrids; igrid++)
	printf("    DERIV G.%d",igrid);
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
#endif

  if (!ThisTask)
    {
      printf("\n");  
      fflush(stdout);
    }  return 0;
}


int deallocate_fft_vectors(int ThisGrid)
{

  pfft_free(cvector_fft[ThisGrid]);
  pfft_free(rvector_fft[ThisGrid]);

  return 0;
}

int reallocate_memory_for_fragmentation()
{

  size_t count_memory, start_memory;
  void *tmp;
#ifdef VERBOSE
  size_t *bcast;
  int bsize, s, n;
  MPI_Status status;
  void *last;
#endif

#ifdef SNAPSHOT
  /* sets zacc to -1 before redistributing it */
  for (int i=0; i<MyGrids[0].total_local_size; i++)
    products[i].zacc=-1;
#endif

#ifndef RECOMPUTE_DISPLACEMENTS

  /* in this case pointers to kdensity and kvector_* are set to zero
     and the corresponding memory is freed; otherwise they are kept */
  kdensity=0x0;
#ifdef TWO_LPT
  kvector_2LPT=0x0;
#ifdef THREE_LPT
  kvector_3LPT_1=0x0;
  kvector_3LPT_2=0x0;
#endif
#endif

#endif

#ifdef TWO_LPT
  source_2LPT=0x0;
#ifdef THREE_LPT
  source_3LPT_1=0x0;
  source_3LPT_2=0x0;
#endif
#endif

  /* these pointers are set to zero in all cases */
  density=0x0;
  first_derivatives=0x0;
  second_derivatives=0x0;
  //seedtable=0x0;

  /* products in the sub-box are stored in the memory that has been freed */

  /* enlarge memory allocation if needed */
  if (memory.frag_allocated > memory.first_allocated)
    {
      tmp=(void*)realloc(main_memory,memory.frag_allocated);
      if (tmp!=0x0)
	{
	  main_memory=tmp;
	  products = (product_data *)main_memory;
	  if (!ThisTask)
	    printf("[%s] Task 0 reallocated memory for %f Gb\n",fdate(),(double)memory.frag_allocated/GBYTE);
	}
      else
	{
	  printf("ERROR on task %d: could not reallocate memory from %zu to %zu bytes\n",ThisTask,
		 memory.first_allocated, memory.frag_allocated);
	  fflush(stdout);
	  return 1;
	}
    }

#ifndef RECOMPUTE_DISPLACEMENTS

  /* the frag structure is located just after memory products */
  frag = (product_data *)(main_memory + memory.prods);
  count_memory =  start_memory = memory.prods + memory.frag_prods;
  ALIGN_MEMORY_BLOCK( count_memory );

#else

  /* including fields_to_keep in this case */
  frag = (product_data *)(main_memory + memory.prods + memory.fields_to_keep);
  count_memory = start_memory = memory.prods + memory.fields_to_keep + memory.frag_prods;
  ALIGN_MEMORY_BLOCK( count_memory );

#endif

  groups = (group_data *)(main_memory + count_memory);
  count_memory += subbox.PredNpeaks*sizeof(group_data);
  ALIGN_MEMORY_BLOCK( count_memory );
  wheretoplace_mycat = (void*)(main_memory + count_memory);
  count_memory += subbox.PredNpeaks*sizeof(histories_data);
  ALIGN_MEMORY_BLOCK( count_memory );
#ifdef PLC
  plcgroups = (plcgroup_data *)(main_memory + count_memory);
  count_memory+=plc.Nmax*sizeof(plcgroup_data);
  ALIGN_MEMORY_BLOCK( count_memory );
#endif
  group_ID = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );
  linking_list = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );

  frag_pos = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );
  indices = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );
  indicesY = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );
  sorted_pos = (int*)(main_memory+count_memory);
  count_memory+=subbox.Nalloc*sizeof(int);
  ALIGN_MEMORY_BLOCK( count_memory );
  frag_map = (unsigned int*)(main_memory+count_memory);
  count_memory+=subbox.maplength*sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK( count_memory );
  frag_map_update = (unsigned int*)(main_memory+count_memory);
  count_memory+=subbox.maplength*sizeof(unsigned int);
  ALIGN_MEMORY_BLOCK( count_memory );

  memset((char*)frag, 0, count_memory-start_memory);

  return 0;
}
