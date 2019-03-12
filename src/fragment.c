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
#include "fragment.h"
#include <sys/types.h>
#include <sys/stat.h>

int init_fragment_1(void);
int init_fragment_2(int *);
int count_peaks();
void set_fragment_parameters(int);
void write_map(int);
int create_map(void);
void reorder(int *, int);
void reorder_nofrag(int *, int); 

//#define READ_FRAGMENT_PARAMS_FROM_FILE

// HO LEVATO LIST_OF_FRAGMENT_PARAMETERS


int ngroups;

void set_fragment_parameters(int order)
{

#ifndef TWO_LPT
  order=1;
#else
#ifndef THREE_LPT
  if (order>2)
    order=2;
#else
  if (order>3)
    order=3;
#endif
#endif

  /* Parameters for the simulation */

  f_200   = 0.171;
  switch(order)
    {
    case 1:
#ifdef USE_SIM_PARAMS
      f_m = f_a = 0.495;
      f_rm    = -0.075;
      espo    = 0.852;
      f_ra    = 0.500;
#else
      f_m = f_a = 0.505;
      f_rm    = 0.000;
      espo    = 0.820;
      f_ra    = 0.300;
#endif  
      sigmaD0 = 1.7;
      break;

    case 2:
#ifdef USE_SIM_PARAMS
      f_m = f_a = 0.475;
      f_rm    = -0.020;
      espo    = 0.780;
      f_ra    = 0.650;
#else
      f_m = f_a = 0.501;
      f_rm    = 0.052;
      espo    = 0.745;
      f_ra    = 0.334;
#endif
      sigmaD0 = 1.5;
      break;

    case 3:
#ifdef USE_SIM_PARAMS
      f_m = f_a = 0.455;
      f_rm    = 0.000;
      espo    = 0.755;
      f_ra    = 0.700;
#else
      f_m = f_a = 0.5024;
      f_rm    = 0.1475;
      espo    = 0.6852;
      f_ra    = 0.4584;
#endif
      sigmaD0 = 1.2;
      break;
 
    default:
      break;
    }

}


static int index_compare_F(const void * a,const void * b)
{
  if ( frag[*((int *)a)].Fmax == frag[*((int *)b)].Fmax )
    return 0;
  else
    if ( frag[*((int *)a)].Fmax > frag[*((int *)b)].Fmax )
      return -1;
    else
      return 1;
}


static int index_compare_P(const void * a,const void * b)
{
  if ( frag_pos[*((int *)a)] == frag_pos[*((int *)b)] )
    return 0;
  else
    if ( frag_pos[*((int *)a)] < frag_pos[*((int *)b)] )
      return -1;
    else
      return 1;
}


int fragment()
{
  int Npeaks, Ngood, tot_good, i, turn, nadd[2], nadd_all[2];
  double tmp;

  /* timing */
  cputime.frag=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] Second part: fragmentation of the collapsed medium\n",fdate());

  if (init_fragment_1())
    return 1;

  cputime.group=cputime.sort=cputime.distr=0.0;
#ifdef PLC
  cputime.plc=0.0;
#endif

  for (ThisSlice=0; ThisSlice<NSlices; ThisSlice++)
    {
      if (!ThisTask)
	printf("\n[%s] Starting slice N. %d of %d\n",fdate(),ThisSlice+1,NSlices);

      for (turn=0; turn<2; turn++)
	{
	  tmp=MPI_Wtime();

	  if (!turn)
	    {
	      if (!ThisTask)
		printf("[%s] Creating map of needed particles\n",fdate());
	      if (create_map())
		return 1;
	      /* the rest of the first turn is useless if there are PBCs on all dimensions (scalar run) */
	      if (subbox.pbc_x && subbox.pbc_y && subbox.pbc_z)
		continue;
	    }
	  else if (!(subbox.pbc_x && subbox.pbc_y && subbox.pbc_z))
	    {
	      if (!ThisTask)
		printf("[%s] Updating map of needed particles\n",fdate());
	      update_map(nadd);
	      
	      MPI_Reduce(nadd, nadd_all, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	      
	      if (!ThisTask)
		{
		  printf("[%s] Requesting %d particles from the boundary layer\n",fdate(),nadd_all[0]);
		  if (nadd_all[1])
		    printf("[%s] WARNING: %d particles requested beyond the boundary layer\n",fdate(),nadd_all[1]);
		}
	    }

	  if (!ThisTask)
	    printf("[%s] Starting %s re-distribution of products\n",fdate(),(turn?"second":"first"));
	  if (distribute())
	    return 1;

	  if (turn)
	    {
	  //LEVARE
	  for (i=0; i<NTasks; i++)
	    {
	      if (i==ThisTask)
		printf("  Task %d: %10d %10d %10f %10Lu %10d %10d\n",
		       ThisTask,subbox.Ngood,subbox.Nstored,
		       (float)subbox.Nstored/(float)subbox.Ngood,
		       MyGrids[0].Ntotal,nadd[0],nadd[1]);
	      fflush(stdout);
	      MPI_Barrier(MPI_COMM_WORLD);
	    }}

	  /* sets the map of potentially loaded particles */
	  if (!turn)
	    for (i=0; i<subbox.maplength; i++)
	      frag_map[i]=frag_map_update[i];
	  else
	    for (i=0; i<subbox.maplength; i++)
	      frag_map[i]|=frag_map_update[i];

	  // DEBUG
	  write_map(turn);

	  tmp=MPI_Wtime()-tmp;
	  cputime.distr += tmp;

	  nadd[0]=subbox.Nstored;
	  MPI_Reduce(nadd, nadd_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	  if (!ThisTask)
	    printf("[%s] %s re-distribution of Fmax done, %d particles stored by all tasks, overhead: %f, cputime = %14.6f\n",
		   fdate(),(turn?"Second":"First"),nadd_all[0],
		   (float)nadd_all[0]/(float)MyGrids[0].Ntotal,
		   tmp);

	  tmp=MPI_Wtime();
	  if (!ThisTask)
	    printf("[%s] Starting sorting\n",fdate());

	  /* sort particles in order of descending Fmax */
	  for (i=0; i<subbox.Nstored; i++)
	    *(indices+i)=i;
	  qsort((void *)indices, subbox.Nstored, sizeof(int), index_compare_F);

	  /* this is needed to reorder */
	  for (i=0; i<subbox.Nstored; i++)
	    indicesY[indices[i]]=i;

	  /* reorder the frag data structure and frag_pos in order of descending Fmax */
	  reorder(indicesY, subbox.Nstored);

	  /* sort particles in order of ascending position */
	  qsort((void *)indices, subbox.Nstored, sizeof(int), index_compare_P);
	  /* create a vector sorted_pos with sorted positions and pointers to it */
	  for (i=0; i<subbox.Nstored; i++)
	    sorted_pos[i]=frag_pos[i];

	  /* reorder the sorted vector */
	  for (i=0; i<subbox.Nstored; i++)
	    indicesY[indices[i]]=i;
	  reorder_nofrag(indicesY,subbox.Nstored);

	  /* checks that no particle is replicated TO BE REMOVED */
	  for (i=0; i<subbox.Nstored-1; i++)
	    if (sorted_pos[i]==sorted_pos[i+1])
	      printf("REPLICA!!!! %d\n",sorted_pos[i]);

	  tmp=MPI_Wtime()-tmp;
	  cputime.sort+=tmp;
	  if (!ThisTask)
	    printf("[%s] Sorting done, total cputime = %14.6f\n",fdate(),tmp);

	  /* counts the number of peaks in the sub-volume to know the number of groups */
	  Npeaks=count_peaks(&Ngood);

	  /* the number of peaks was supposed to be at most 1/10 of the number of particles
	     (we need space for Npeaks+2 groups) */
	  if (Npeaks+2 > subbox.PredNpeaks)
	    {
	      printf("ERROR on task %d: surprisingly, the number of peaks %d exceeds Npart/10 (%d)\n",
		     ThisTask,Npeaks,subbox.PredNpeaks);
	      fflush(stdout);
	      return 1;
	    }

	  if (MPI_Reduce(&Ngood, &tot_good, 1, 
			 MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
	    {
	      printf("ERROR on task %d: init_fragment failed an MPI_Reduce\n",ThisTask);
	      fflush(stdout);
	      return 1;
	    }

	  if (!ThisTask)
	    printf("[%s] Task 0 found %d peaks (predicted: %d), in the well resolved region: %d. Total number of peaks: %d\n",
		   fdate(),Npeaks,subbox.PredNpeaks,Ngood,tot_good);

	  tmp=MPI_Wtime();

	  /* sets arrays to zero */
	  memset(group_ID,0,subbox.Nalloc*sizeof(int));
	  memset(linking_list,0,subbox.Nalloc*sizeof(int));
	  memset(groups,0,subbox.PredNpeaks*sizeof(group_data));
	  memset(wheretoplace_mycat,0,subbox.PredNpeaks*sizeof(histories_data));
#ifdef PLC
	  memset(plcgroups,0,plc.Nmax*sizeof(plcgroup_data));
#endif

	  if (!turn)
	    {
	      /* quickly generates the group catalogue */
	      if (quick_build_groups(Npeaks))
		return 1;
	    }
	  else
	    {
	      /* generates the group catalog */
	      if (build_groups(Npeaks))
		return 1;
	    }

	  cputime.group+=MPI_Wtime()-tmp;
	}

      /* next slice */
      subbox.mybox_z += subbox.nbox_z_thisslice;
      subbox.start_z = find_start( MyGrids[0].GSglobal[_z_], subbox.nbox_z_allslices, subbox.mybox_z);
      subbox.Lgrid_z = find_length(MyGrids[0].GSglobal[_z_], subbox.nbox_z_allslices, subbox.mybox_z);
      subbox.stabl_z = subbox.start_z - subbox.safe_z;

    }
  /* cpu time for fragmentation */
  cputime.frag = MPI_Wtime() - cputime.frag;

  if (!ThisTask)
    printf("[%s] Finishing fragment, total cputime = %14.6f\n",fdate(),cputime.frag);

  return 0;
}



static int compare_search(const void * a,const void * b)
{
  if ( *((int *)a) == *((int *)b) )
    return 0;
  else
    if ( *((int *)a) < *((int *)b) )
      return -1;
    else
      return 1;
}

/* reorders frag_pos and frag after they have been index-sorted */
void reorder(int *ind, int n) 
{
  int i,oldf,next,tmp;
  product_data oldF,tmpF;

  memset(&oldF,0,sizeof(product_data));
  memset(&tmpF,0,sizeof(product_data));
  for (i=0; i<n; i++)
    {
      if (ind[i]!=i)
	{
	  oldf=frag_pos[i];
	  memcpy(&oldF,&frag[i],sizeof(product_data));
	  next=ind[i];
	  while (next != i) 
	    {
	      tmp=frag_pos[next];	    
	      frag_pos[next]=oldf;
	      oldf=tmp;
	      memcpy(&tmpF,&frag[next],sizeof(product_data));
	      memcpy(&frag[next],&oldF,sizeof(product_data));
	      memcpy(&oldF,&tmpF,sizeof(product_data));
	      tmp=ind[next];
	      ind[next]=next;
	      next=tmp;
	    }
	  frag_pos[i]=oldf;
	  memcpy(&frag[i],&oldF,sizeof(product_data));
	  ind[i]=next;
	}
    } 
}

/* reorders only frag_pos after it has been index-sorted */
void reorder_nofrag(int *ind, int n) 
{
  int i,oldf,next,tmp;

  for (i=0; i<n; i++)
    {
      if (ind[i]!=i)
	{
	  oldf=sorted_pos[i];
	  next=ind[i];
	  while (next != i) 
	    {
	      tmp=sorted_pos[next];
	      sorted_pos[next]=oldf;
	      oldf=tmp;
	      tmp=ind[next];
	      ind[next]=next;
	      next=tmp;
	    }
	  sorted_pos[i]=oldf;
	  ind[i]=next;
	}
    } 
}


int init_fragment_1(void)
{

  /* Initializes all quantities for the run */

  //int i;

#ifdef READ_FRAGMENT_PARAMS_FROM_FILE

  char filename[BLENGTH];
  FILE *file;
  int n;

  if (!ThisTask)
    {
      sprintf(filename,"%s.fragment_params",params.RunFlag);
      printf("Fragment parameters will be read from file %s\n",filename);

      file=fopen(filename,"r");
      if (file==0x0)
	{
	  printf("ERROR: fragment parameter file %s not found\n",filename);
	  printf("I will use standard values for the parameters\n");
	  set_fragment_parameters(ORDER_FOR_GROUPS);
	}
      else
	{
	  n=fscanf(file,"%lf %lf %lf %lf %lf %lf", &f_m, &f_rm, &espo, &f_a, &f_ra, &f_200);
	  if (n!=6)
	    {
	      printf("ERROR: fragment parameters not read from file %s\n",filename);
	      printf("I will use standard values for the parameters\n");
	      set_fragment_parameters(ORDER_FOR_GROUPS);
	    }
	  fclose(file);
	}
    }

  MPI_Bcast(&f_m   , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&f_rm  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&espo  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&f_a   , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&f_ra  , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&f_200 , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#else

  set_fragment_parameters(ORDER_FOR_GROUPS);

#endif


  /* reallocates the memory, for the part needed to redistribution */
  if (reallocate_memory_for_fragmentation())
    return 1;

  return 0;
}


int find_location(int i,int j,int k)
{
  /* this function finds the location in the frag and frag_pos vectors of a given particle */

  int pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x;
  int *nn=bsearch((void*)&pos,(void*)sorted_pos,(size_t)subbox.Nstored,sizeof(int),compare_search);
  if (nn!=0x0)
    return nn-sorted_pos;
  else
    return -1; /* if the particle is not stored, it returns -1 */
}
    

int count_peaks(int *ngood)
{

  /* counts the number of Fmax peaks down to the final redshift
     to know the total number of groups */

  int iz,i,j,k,kk,nn,ngroups_tot,peak_cond,i1,j1,k1,pos;

  FILE *fd=fopen("peaks.dat","w");

  ngroups_tot=0;     // number of groups, group 1 is the filament group
  *ngood=0;          // number of groups out of the safety boundary

  for (iz=0; iz<subbox.Nstored; iz++)
    {

      /* position on the local box */
      i  = frag_pos[iz]%subbox.Lgwbl_x;
      kk = frag_pos[iz]/subbox.Lgwbl_x;
      j  = kk%subbox.Lgwbl_y;
      k  = kk/subbox.Lgwbl_y;

      /* avoid borders */
      if ( !subbox.pbc_x && (i==0 || i==subbox.Lgwbl_x-1) ) continue;
      if ( !subbox.pbc_y && (j==0 || j==subbox.Lgwbl_y-1) ) continue;
      if ( !subbox.pbc_z && (k==0 || k==subbox.Lgwbl_z-1) ) continue;

      /* peak condition */
      peak_cond=1;
      for (nn=0; nn<6; nn++)
	{
	  switch (nn)
	    {
	    case 0:
	      i1=( subbox.pbc_x && i==0 ? subbox.Lgwbl_x-1 : i-1 );
	      j1=j;
	      k1=k;
	      break;
	    case 1:
	      i1=( subbox.pbc_x && i==subbox.Lgwbl_x-1 ? 0 : i+1 );
	      j1=j;
	      k1=k;
	      break;
	    case 2:
	      i1=i;
	      j1=( subbox.pbc_y && j==0 ? subbox.Lgwbl_y-1 : j-1 );
	      k1=k;
	      break;
	    case 3:
	      i1=i;
	      j1=( subbox.pbc_y && j==subbox.Lgwbl_y-1 ? 0 : j+1 );
	      k1=k;
	      break;
	    case 4:
	      i1=i;
	      j1=j;
	      k1=( subbox.pbc_z && k==0 ? subbox.Lgwbl_z-1 : k-1 );
	      break;
	    case 5:
	      i1=i;
	      j1=j;
	      k1=( subbox.pbc_z && k==subbox.Lgwbl_z-1 ? 0 : k+1 );
	      break;
	    }

	  pos = find_location(i1,j1,k1);
	  if (pos>=0)
	    peak_cond &= (frag[iz].Fmax > frag[indices[pos]].Fmax);

	  if (!peak_cond)
	    break;

	}

      if (peak_cond) // && frag[iz].Fmax >= outputs.Flast)  QUESTA NON SERVE PIU`
	{
	  ngroups_tot++;
	  if ( i>=subbox.safe_x && i<subbox.Lgwbl_x-subbox.safe_x &&
	       j>=subbox.safe_y && j<subbox.Lgwbl_y-subbox.safe_y &&
	       k>=subbox.safe_z && k<subbox.Lgwbl_z-subbox.safe_z)
	    (*ngood)++;
	  fprintf(fd," %2d %2d %2d   %12.10f\n",i,j,k,frag[iz].Fmax);

	}
    }

  fclose(fd);

  return ngroups_tot;
}


int create_map()
{
  int i,j,k,i1,i2,j1,j2,k1,k2;

  memset(frag_map_update, 0, subbox.maplength*sizeof(unsigned int));
  if (!subbox.pbc_x)
    {
      i1=subbox.safe_x-1;
      i2=subbox.Lgrid_x+subbox.safe_x+1;
    }
  else
    {
      i1=0;
      i2=subbox.Lgrid_x;
    }
  if (!subbox.pbc_y)
    {
      j1=subbox.safe_y-1;
      j2=subbox.Lgrid_y+subbox.safe_y+1;
    }
  else
    {
      j1=0;
      j2=subbox.Lgrid_y;
    }
  if (!subbox.pbc_z)
    {
      k1=subbox.safe_z-1;
      k2=subbox.Lgrid_z+subbox.safe_z+1;
    }
  else
    {
      k1=0;
      k2=subbox.Lgrid_z;
    }

  for (i=i1; i<i2; i++)
    for (j=j1; j<j2; j++)
      for (k=k1; k<k2; k++)
	set_mapup_bit(i,j,k);

  return 0;
}


void set_mapup_bit(int i, int j, int k)
{
  /* this operates on frag_map_update */
  unsigned int pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x;
  frag_map_update[pos/UINTLEN]|=(1<<pos%UINTLEN);
}

int get_mapup_bit(unsigned int pos)
{
  /* this operates on frag_map_update */
  /* unsigned int pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x; */
  unsigned int rem = pos%UINTLEN;
  return (frag_map_update[pos/UINTLEN] & (1<<rem))>>rem;
}

int get_map_bit(int i, int j, int k)
{
  /* this operates on frag_map */
  unsigned int pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x;
  unsigned int rem = pos%UINTLEN;
  return (frag_map[pos/UINTLEN] & (1<<rem))>>rem;
}

void write_map(int turn)
{
  int i,j,k;
  unsigned int pos,rem;
  char fname[BLENGTH];
  sprintf(fname,"map_task%d_turn%d.txt",ThisTask,turn);
  FILE *fd=fopen(fname,"w"); 
  for (k=0; k<subbox.Lgwbl_z; k++)
    {
      fprintf(fd,"k=%d\n",k);
      for (j=0; j<subbox.Lgwbl_y; j++)
	{
	  for (i=0; i<subbox.Lgwbl_x; i++)
	    {
	      pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x;
	      rem = pos%UINTLEN;
	      fprintf(fd,"%1d",(frag_map[pos/UINTLEN] & (1<<rem))>>rem);
	    }
	  fprintf(fd,"\n");
	}
    }
  fclose(fd);
  sprintf(fname,"mapup_task%d_turn%d.txt",ThisTask,turn);
  fd=fopen(fname,"w"); 
  for (k=0; k<subbox.Lgwbl_z; k++)
    {
      fprintf(fd,"k=%d\n",k);
      for (j=0; j<subbox.Lgwbl_y; j++)
	{
	  for (i=0; i<subbox.Lgwbl_x; i++)
	    {
	      pos = i + (j + k*subbox.Lgwbl_y)*subbox.Lgwbl_x;
	      rem = pos%UINTLEN;
	      fprintf(fd,"%1d",(frag_map_update[pos/UINTLEN] & (1<<rem))>>rem);
	    }
	  fprintf(fd,"\n");
	}
    }
  fclose(fd);
}


/*
  int peak, ip, jp, kp, pos;
  // check: trova i picchi 
  int counter=0;
  for (i=0; i<subbox.Nstored; i++)
    {

      ip = frag_pos[i]%subbox.Lgwbl_x;
      jp =(frag_pos[i]/subbox.Lgwbl_x)%subbox.Lgwbl_y;
      kp = frag_pos[i]/subbox.Lgwbl_x/subbox.Lgwbl_y;

	peak=1;
	if (ip==0)
	  pos=find_location(subbox.Lgwbl_x-1,jp,kp);
	else
	  pos=find_location(ip-1,jp,kp);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);

	if (ip==subbox.Lgwbl_x-1)
	  pos=find_location(0,jp,kp);
	else
	  pos=find_location(ip+1,jp,kp);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);

	if (jp==0)
	  pos=find_location(ip,subbox.Lgwbl_y-1,kp);
	else
	  pos=find_location(ip,jp-1,kp);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);

	if (jp==subbox.Lgwbl_y-1)
	  pos=find_location(ip,0,kp);
	else
	  pos=find_location(ip,jp+1,kp);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);

	if (kp==0)
	  pos=find_location(ip,jp,subbox.Lgwbl_z-1);
	else
	  pos=find_location(ip,jp,kp-1);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);

	if (kp==subbox.Lgwbl_z-1)
	  pos=find_location(ip,jp,0);
	else
	  pos=find_location(ip,jp,kp+1);
	if (pos>=0)
	  peak&=(frag[indices[pos]].Fmax<frag[i].Fmax);
	
	if (peak)
	  {
	    counter++;
	  }
    }
  printf("Numero di picchi: %d\n",counter);


*/


