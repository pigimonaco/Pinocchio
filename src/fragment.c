/* ######HEADER###### */

#include "pinocchio.h"
#include <sys/types.h>
#include <sys/stat.h>

//#define DEBUG

int fragment(void);
int count_peaks(int *);
void set_fragment_parameters(int);
void write_map(int);
int create_map(void);
void reorder(int *, int);
void reorder_nofrag(int *, int); 
void sort_and_organize(void);
int init_segmentation(void);
#ifdef RECOMPUTE_DISPLACEMENTS
int compute_future_LPT_displacements(double, int);
int shift_all_displacements(void);
int recompute_group_velocities(void);
#endif


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


int fragment_driver()
{
  /* This routine can be edited to call fragmentation for many combinations 
     of fragment parameters */

  if (!ThisTask)
    printf("\n[%s] Second part: fragmentation of the collapsed medium\n",fdate());

  /* parameters are assigned in the standard way */
  set_fragment_parameters(ORDER_FOR_GROUPS);
  
  /* reallocates the memory */
  if (reallocate_memory_for_fragmentation())
    return 1;

  if (fragment())
    return 1;

  return 0;
}


int fragment()
{
  /* This function is the driver for fragmentation of collapsed medium 
     and construction of halo catalogs */

  int Npeaks, Ngood;
#ifndef CLASSIC_FRAGMENTATION
  unsigned int nadd[2];
#endif
  double BestPredPeakFactor;
  unsigned long long mynadd[2], nadd_all[2];
  double tmp;
  int mysegment;


  /* this initializes the segmentation of the fragmentation process */
  if (init_segmentation())
    return 1;

#ifdef RECOMPUTE_DISPLACEMENTS

  /* computation of displacements for the first redshift */
  if (!ThisTask)
    printf("\n[%s] Computing displacements for first redshift %f\n",fdate(),Segment.z[0]);
  if (compute_future_LPT_displacements(Segment.z[0],0))
    return 1;

  /* computation of displacements for the second redshift */
  if (!ThisTask)
    printf("\n[%s] Computing displacements for redshift %f\n",fdate(),Segment.z[1]);
  if (compute_future_LPT_displacements(Segment.z[1],1))
    return 1;

#endif

  /* timing */
  cputime.frag=MPI_Wtime();
  cputime.group=cputime.sort=cputime.distr=0.0;
#ifdef PLC
  cputime.plc=0.0;
#endif


#ifdef CLASSIC_FRAGMENTATION
  /*************************************************************/
  /* This implements the classic fragmentation, without slices */
  /*************************************************************/

  /* redistribution of products */
  if (!ThisTask)
    printf("[%s] Starting re-distribution of products\n",fdate());

  tmp=MPI_Wtime();

  if (distribute())
    return 1;

  tmp=MPI_Wtime()-tmp;
  cputime.distr += tmp;

  if (!ThisTask)
    printf("[%s] re-distribution of Fmax done, cputime = %14.6f\n", fdate(), tmp);

  /* counts the number of peaks in the sub-volume to know the number of groups */
  Npeaks=count_peaks(&Ngood);

  mynadd[0]=Ngood;
  MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      printf("[%s] Task 0 found %d peaks, %d in the well resolved region. Total number of peaks: %Ld\n",
	     fdate(),Npeaks,Ngood,nadd_all[0]);
    }

  /* the number of peaks was supposed to be at most 1/10 of the number of particles
     (we need space for Npeaks+2 groups) */
  if (Npeaks+2 > subbox.PredNpeaks)
    {
      printf("ERROR on task %d: the number of peaks %d exceeds the predicted one (%d)\n",
	     ThisTask,Npeaks,subbox.PredNpeaks);
      printf("      Please increase PredPeakFactor to at least %4.2f and restart\n",
	     (double)Npeaks / (double)MyGrids[0].ParticlesPerTask * 6.0 + 0.01);
      fflush(stdout);
      return 1;
    }

  /* Gives the minimal PredPeakFactor */
  mynadd[0]=Npeaks;
  MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  BestPredPeakFactor=(double)nadd_all[0]/(double)MyGrids[0].ParticlesPerTask*6.0;

  /* sets arrays to zero */
  memset(group_ID,0,subbox.Nalloc*sizeof(int));
  memset(linking_list,0,subbox.Nalloc*sizeof(int));
  memset(groups,0,subbox.PredNpeaks*sizeof(group_data));
  memset(wheretoplace_mycat,0,subbox.PredNpeaks*sizeof(histories_data));
#ifdef PLC
  memset(plcgroups,0,plc.Nmax*sizeof(plcgroup_data));
#endif

  tmp=MPI_Wtime();
  if (!ThisTask)
    printf("[%s] Starting sorting\n",fdate());

  for (int i=0; i<subbox.Npart; i++)
    *(indices+i)=i;
  qsort((void *)indices, subbox.Npart, sizeof(int), index_compare_F);

  tmp=MPI_Wtime()-tmp;
  cputime.sort+=tmp;
  if (!ThisTask)
    printf("[%s] Sorting done, total cputime = %14.6f\n",fdate(),tmp);

  /* full fragmentation is done segmenting the redshift interval */
  for (mysegment=0; mysegment<Segment.n; mysegment++)
    {
      Segment.mine=mysegment;

#ifdef RECOMPUTE_DISPLACEMENTS
      // tanto cambia...
#endif

      tmp=MPI_Wtime();

      if (build_groups(Npeaks,Segment.z[mysegment],(mysegment==0)))
	return 1;

      tmp=MPI_Wtime()-tmp;
      if (!ThisTask)
	printf("[%s] Fragmentation done to redshift %6.4f, cputime = %14.6f\n",
	       fdate(),Segment.z[mysegment],tmp);
      cputime.group+=tmp;
    }

#else
  /*********************************************/
  /* This implements the updated fragmentation */
  /*********************************************/
  
  int all_pbc = (subbox.pbc[_x_] && subbox.pbc[_y_] && subbox.pbc[_z_]);

  /* fragmentation is performed twice;
     in the first turn a quick fragmentation is performed to locate
     halos and determine the interesting particles; velocities are not updated;
     in the second turn the full fragmentation is performed, segmentation is applied
     and velocities are updated */
  for (int turn=0; turn<2; turn++)
    {

      if (!turn)
	{

	  /* it creates the map by setting to 1 the well-resolved region */
	  if (!ThisTask)
	    printf("[%s] Creating map of needed particles\n",fdate());

	  if (create_map())
	    return 1;

	  /* the rest of the first turn is skipped if there are PBCs in all directions */
	  if (all_pbc)
	    continue;

	}
      else if (!all_pbc)
	{

	  /* it updates the map by adding spheres around halos in the boundary layer */
	  if (!ThisTask)
	    printf("[%s] Updating map of needed particles\n",fdate());

	  update_map(nadd);

	  mynadd[0]=nadd[0];
	  mynadd[1]=nadd[1];
	  MPI_Reduce(mynadd, nadd_all, 2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

	  if (!ThisTask)
	    {
	      printf("[%s] Requesting %Ld particles from the boundary layer\n",fdate(),nadd_all[0]);
	      if (nadd_all[1])
		printf("WARNING: %Ld requested particles lie beyond the boundary layer, some halos may be inaccurate\n",nadd_all[1]);
	    }
	}

      /* redistribution of products for queried particles */
      if (!ThisTask)
	printf("[%s] Starting %s re-distribution of products\n",fdate(),(turn?"second":"first"));

      tmp=MPI_Wtime();

      if (distribute())
	return 1;

      if (subbox.Nneeded>subbox.Nalloc)
	{
	  if (turn)
	    {
	      printf("CRITICAL WARNING in Task %d: it allocated %d but received %d particles\n",
		     ThisTask, subbox.Nalloc, subbox.Nneeded);
	      printf("The required overhead is %f\n",(float)subbox.Nneeded/(float)MyGrids[0].ParticlesPerTask);
	      printf("Please increase MaxMemPerParticle by at least %d and start again\n",
		     (int)((float)(subbox.Nneeded-subbox.Nalloc)/(float)MyGrids[0].ParticlesPerTask * 
			   (sizeof(product_data) + FRAGFIELDS * sizeof(int)))+1+(int)params.MaxMemPerParticle);
	      if (params.ExitIfExtraParticles)
		return 1;
	    }
	  else
	    {
	      printf("ERROR in Task %d: the number of allocated particles (%d) is WAY too small! I need at least %d\n",
		     ThisTask, subbox.Nalloc, subbox.Nneeded);
	      printf("The required overhead is %f\n",(float)subbox.Nneeded/(float)MyGrids[0].ParticlesPerTask);
	      printf("Please increase MaxMemPerParticle to at least %d and start again\n",
		     (int)((float)(subbox.Nneeded-subbox.Nalloc)/(float)MyGrids[0].ParticlesPerTask * 
			   (sizeof(product_data) + FRAGFIELDS * sizeof(int)))+1+(int)params.MaxMemPerParticle);
	      return 1;
	    }
	}

      tmp=MPI_Wtime()-tmp;
      cputime.distr += tmp;

      mynadd[0]=subbox.Nstored;
      MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

      if (!ThisTask)
	printf("[%s] %s re-distribution of Fmax done, %Ld particles stored by all tasks, average overhead: %f, cputime = %14.6f\n",
	       fdate(), (turn?"Second":"First"), nadd_all[0], 
	       (float)nadd_all[0]/(float)MyGrids[0].Ntotal, tmp);

      mynadd[0]=subbox.Nstored;
      MPI_Reduce(mynadd, nadd_all  , 1, MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
      MPI_Reduce(mynadd, nadd_all+1, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

      if (!ThisTask)
	printf("[%s] Smallest and largest overhead: %f, %f\n",fdate(),
	       (float)nadd_all[0]/(float)MyGrids[0].ParticlesPerTask,
	       (float)nadd_all[1]/(float)MyGrids[0].ParticlesPerTask);


      /* sets or updates the map of potentially loaded particles */
      if (!turn)
	for (int i=0; i<subbox.maplength; i++)
	  frag_map[i]=frag_map_update[i];
      else
	for (int i=0; i<subbox.maplength; i++)
	  frag_map[i]|=frag_map_update[i];

#ifdef DEBUG
      write_map(turn);
#endif

      /* this part sorts particles and reorganizes index arrays */
      sort_and_organize();

      /* counts the number of peaks in the sub-volume to know the number of groups */
      Npeaks=count_peaks(&Ngood);

      mynadd[0]=Ngood;
      MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

      if (!ThisTask)
	{
	  printf("[%s] Task 0 found %d peaks, %d in the well resolved region. Total number of peaks: %Ld\n",
		 fdate(),Npeaks,Ngood,nadd_all[0]);
	}

       /* the number of peaks was supposed to be at most 1/10 of the number of particles
	 (we need space for Npeaks+2 groups) */
      if (Npeaks+2 > subbox.PredNpeaks)
	{
	  printf("ERROR on task %d: the number of peaks %d exceeds the predicted one (%d)\n",
		 ThisTask,Npeaks,subbox.PredNpeaks);
	  printf("      Please increase PredPeakFactor and restart\n");
	  fflush(stdout);
	  return 1;
	}

      if (turn)
	{
	  /* Gives the minimal PredPeakFactor */
	  mynadd[0]=Npeaks;
	  MPI_Reduce(mynadd, nadd_all, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);

	  BestPredPeakFactor=(double)nadd_all[0]/(double)MyGrids[0].ParticlesPerTask*6.0;
	}


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

	  /* quick fragmentation is done in a shot */
	  tmp=MPI_Wtime();

	  /* quickly generates the group catalogue */
	  if (quick_build_groups(Npeaks))
	  return 1;

	  tmp=MPI_Wtime()-tmp;
	  cputime.group+=tmp;
	  if (!ThisTask)
	    printf("[%s] Quick fragmentation done, cputime = %14.6f\n",fdate(),tmp);

	}
      else
	{
	  /* full fragmentation is done segmenting the redshift interval */
	  for (mysegment=0; mysegment<Segment.n; mysegment++)
	    {
	      Segment.mine=mysegment;


#ifdef RECOMPUTE_DISPLACEMENTS
	      /* The first segment uses only the first computed velocity
		 The second segment has the two velocities already loaded */
	      if (mysegment>=2)
		{
		  /* computation of displacements for the next redshift */
		  if (!ThisTask)
		    printf("\n[%s] Computing displacements for redshift %f\n",fdate(),Segment.z[mysegment]);
		  if (shift_all_displacements())
		    return 1;

		  if (compute_future_LPT_displacements(Segment.z[mysegment],1))
		    return 1;

		  tmp=MPI_Wtime();
		  if (!ThisTask)
		    printf("[%s] Starting re-distribution of products\n",fdate());

		  full_map(); // QUESTA E` DA SCRIVERE (forse)

		  if (distribute())
		    return 1;

		  tmp=MPI_Wtime()-tmp;
		  cputime.distr += tmp;
		  if (!ThisTask)
		    printf("[%s] Re-distribution of Fmax products done, cputime = %14.6f\n",fdate(),tmp);

		  if (recompute_group_velocities())  // RIVEDERE
		    return 1;
		}
#endif


	      tmp=MPI_Wtime();
      
	      //scrivi(); //LEVARE

	      if (build_groups(Npeaks,Segment.z[mysegment],(mysegment==0)))
		return 1;

	      tmp=MPI_Wtime()-tmp;
	      if (!ThisTask)
		printf("[%s] Fragmentation done to redshift %6.4f, cputime = %14.6f\n",
		       fdate(),Segment.z[mysegment],tmp);
	      cputime.group+=tmp;
	    }
	}
    }

#endif

  /* cpu time for fragmentation */
  cputime.frag = MPI_Wtime() - cputime.frag;

  if (!ThisTask)
    printf("[%s] Finishing fragment, total cputime = %14.6f\n",fdate(),cputime.frag);

#ifdef SNAPSHOT
  if (params.WriteTimelessSnapshot)
    {
      /* redistribution of products */
      if (!ThisTask)
	printf("[%s] Starting to distribute back accretion redshifts\n",fdate());
      tmp=MPI_Wtime();

      if (distribute_back())
	return 1;

      tmp=MPI_Wtime()-tmp;
      cputime.distr += tmp;

      if (!ThisTask)
	printf("[%s] back-distribution of zacc done, cputime = %14.6f\n", fdate(), tmp);

      if (write_timeless_snapshot())
	return 1;
    }
#endif

  if (!ThisTask)
    {
      printf("\n");
      printf("[%s] Minimal memory requirements:",fdate());
      printf("[%s] The PredPeakFactor parameter could have been %5.2f in place of %5.2f\n",
	     fdate(),BestPredPeakFactor,params.PredPeakFactor);
    }

  return 0;
}


void sort_and_organize(void)
{
  double tmp;
  int i;

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

  tmp=MPI_Wtime()-tmp;
  cputime.sort+=tmp;
  if (!ThisTask)
    printf("[%s] Sorting done, total cputime = %14.6f\n",fdate(),tmp);

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


int find_location(int i,int j,int k)
{
  /* this function finds the location in the frag and frag_pos vectors of a given particle */

  int pos = COORD_TO_INDEX(i,j,k,subbox.Lgwbl);
    //i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_];
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

  int iz,i,j,k,nn,ngroups_tot,peak_cond,i1,j1,k1;

#ifdef DEBUG
  FILE *fd=fopen("peaks.dat","w");
#endif

  ngroups_tot=0;     // number of groups, group 1 is the filament group
  *ngood=0;          // number of groups out of the safety boundary

  for (iz=0; iz<subbox.Nstored; iz++)
    {
      /* position on the local box */
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,i,j,k,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],i,j,k,subbox.Lgwbl);
#endif

/* #ifdef CLASSIC_FRAGMENTATION */
/*       i  = iz%subbox.Lgwbl[_x_]; */
/*       kk = iz/subbox.Lgwbl[_x_]; */
/* #else */
/*       i  = frag_pos[iz]%subbox.Lgwbl[_x_]; */
/*       kk = frag_pos[iz]/subbox.Lgwbl[_x_]; */
/* #endif */
/*       j  = kk%subbox.Lgwbl[_y_]; */
/*       k  = kk/subbox.Lgwbl[_y_]; */

      /* avoid borders */
      if ( !subbox.pbc[_x_] && (i==0 || i==subbox.Lgwbl[_x_]-1) ) continue;
      if ( !subbox.pbc[_y_] && (j==0 || j==subbox.Lgwbl[_y_]-1) ) continue;
      if ( !subbox.pbc[_z_] && (k==0 || k==subbox.Lgwbl[_z_]-1) ) continue;

      /* peak condition */
      peak_cond=1;
      for (nn=0; nn<6; nn++)
	{
	  switch (nn)
	    {
	    case 0:
	      i1=( subbox.pbc[_x_] && i==0 ? subbox.Lgwbl[_x_]-1 : i-1 );
	      j1=j;
	      k1=k;
	      break;
	    case 1:
	      i1=( subbox.pbc[_x_] && i==subbox.Lgwbl[_x_]-1 ? 0 : i+1 );
	      j1=j;
	      k1=k;
	      break;
	    case 2:
	      i1=i;
	      j1=( subbox.pbc[_y_] && j==0 ? subbox.Lgwbl[_y_]-1 : j-1 );
	      k1=k;
	      break;
	    case 3:
	      i1=i;
	      j1=( subbox.pbc[_y_] && j==subbox.Lgwbl[_y_]-1 ? 0 : j+1 );
	      k1=k;
	      break;
	    case 4:
	      i1=i;
	      j1=j;
	      k1=( subbox.pbc[_z_] && k==0 ? subbox.Lgwbl[_z_]-1 : k-1 );
	      break;
	    case 5:
	      i1=i;
	      j1=j;
	      k1=( subbox.pbc[_z_] && k==subbox.Lgwbl[_z_]-1 ? 0 : k+1 );
	      break;
	    }

#ifdef CLASSIC_FRAGMENTATION
	  peak_cond &= (frag[iz].Fmax > frag[COORD_TO_INDEX(i1,j1,k1,subbox.Lgwbl)].Fmax);
	  //i1 +j1*subbox.Lgwbl[_x_]+ k1*Lgridxy
#else
	  /* looks for the neighbouring particle in the list */
	  int pos = find_location(i1,j1,k1);
	  if (pos>=0)
	    peak_cond &= (frag[iz].Fmax > frag[indices[pos]].Fmax);
#endif

	  if (!peak_cond)
	    break;

	}

      if (peak_cond)
	{
	  ngroups_tot++;
	  if ( i>=subbox.safe[_x_] && i<subbox.Lgwbl[_x_]-subbox.safe[_x_] &&
	       j>=subbox.safe[_y_] && j<subbox.Lgwbl[_y_]-subbox.safe[_y_] &&
	       k>=subbox.safe[_z_] && k<subbox.Lgwbl[_z_]-subbox.safe[_z_])
	    (*ngood)++;
#ifdef DEBUG
	  fprintf(fd," %2d %2d %2d   %12.10f\n",i,j,k,frag[iz].Fmax);
#endif

	}
    }

#ifdef DEBUG
  fclose(fd);
#endif

  return ngroups_tot;
}


int create_map()
{
  int i,j,k,i1,i2,j1,j2,k1,k2;

  /* sets to 1 all particles in the well-resolved region plus one row for each side (without PBCs) */
  memset(frag_map_update, 0, subbox.maplength*sizeof(unsigned int));
  if (!subbox.pbc[_x_])
    {
      i1=subbox.safe[_x_]-1;
      i2=subbox.Lgrid[_x_]+subbox.safe[_x_]+1;
    }
  else
    {
      i1=0;
      i2=subbox.Lgrid[_x_];
    }
  if (!subbox.pbc[_y_])
    {
      j1=subbox.safe[_y_]-1;
      j2=subbox.Lgrid[_y_]+subbox.safe[_y_]+1;
    }
  else
    {
      j1=0;
      j2=subbox.Lgrid[_y_];
    }
  if (!subbox.pbc[_z_])
    {
      k1=subbox.safe[_z_]-1;
      k2=subbox.Lgrid[_z_]+subbox.safe[_z_]+1;
    }
  else
    {
      k1=0;
      k2=subbox.Lgrid[_z_];
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
  unsigned int pos = COORD_TO_INDEX(i,j,k,subbox.Lgwbl);
  //i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_];
  frag_map_update[pos/UINTLEN]|=(1<<pos%UINTLEN);
}

int get_mapup_bit(unsigned int pos)
{
  /* this operates on frag_map_update */
  /* unsigned int pos = i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_]; */
  unsigned int rem = pos%UINTLEN;
  return (frag_map_update[pos/UINTLEN] & (1<<rem))>>rem;
}

int get_map_bit(int i, int j, int k)
{
  /* this operates on frag_map */
  unsigned int pos = COORD_TO_INDEX(i,j,k,subbox.Lgwbl);
  //i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_];
  unsigned int rem = pos%UINTLEN;
  return (frag_map[pos/UINTLEN] & (1<<rem))>>rem;
}

void write_map(int turn)
{
  int i,j,k;
  unsigned int pos,rem;
  char fname[LBLENGTH];
  sprintf(fname,"map_task%d_turn%d.txt",ThisTask,turn);
  FILE *fd=fopen(fname,"w");
  /* this loop is not following memory... */
  for (k=0; k<subbox.Lgwbl[_z_]; k++)
    {
      fprintf(fd,"k=%d\n",k);
      for (j=0; j<subbox.Lgwbl[_y_]; j++)
	{
	  for (i=0; i<subbox.Lgwbl[_x_]; i++)
	    {
	      pos = COORD_TO_INDEX(i,j,k,subbox.Lgwbl); 
	      //i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_];
	      rem = pos%UINTLEN;
	      fprintf(fd,"%1d",(frag_map[pos/UINTLEN] & (1<<rem))>>rem);
	    }
	  fprintf(fd,"\n");
	}
    }
  fclose(fd);
  sprintf(fname,"mapup_task%d_turn%d.txt",ThisTask,turn);
  fd=fopen(fname,"w"); 
  for (k=0; k<subbox.Lgwbl[_z_]; k++)
    {
      fprintf(fd,"k=%d\n",k);
      for (j=0; j<subbox.Lgwbl[_y_]; j++)
	{
	  for (i=0; i<subbox.Lgwbl[_x_]; i++)
	    {
	      pos =  COORD_TO_INDEX(i,j,k,subbox.Lgwbl); 
	      //i + (j + k*subbox.Lgwbl[_y_])*subbox.Lgwbl[_x_];
	      rem = pos%UINTLEN;
	      fprintf(fd,"%1d",(frag_map_update[pos/UINTLEN] & (1<<rem))>>rem);
	    }
	  fprintf(fd,"\n");
	}
    }
  fclose(fd);
}



int init_segmentation(void)
{
  /* Fragmentation is segmented into redshift bins, 
     that may be a unique one from infinity to the last redshit */

#ifdef RECOMPUTE_DISPLACEMENTS
  /* here the segmentation of the fragmentation process is defined */
  /* for now, segmentation is taken from the outputs file */
  Segment.n=outputs.n;
  for (int i=0; i<outputs.n; i++)
    Segment.z[i]=outputs.z[i];
  Segment.mine=1;
  return 0;
#else
  Segment.n=1;
  Segment.z[0]=outputs.zlast;
  Segment.mine=0;
  return 0;
#endif

}



#ifdef RECOMPUTE_DISPLACEMENTS
// DA RIVEDERE

int compute_future_LPT_displacements(double z, int after)
{

  int ia,local_x,local_y,local_z,index;
  double cputmp;
  char tag[15];
  sprintf(tag,"_%6.4f",z);

#ifdef SCALE_DEPENDENT
  ScaleDep.order=1;
  ScaleDep.redshift=z;
#endif

  cputmp=MPI_Wtime();

  /* Zeldovich displacements */
  if (!ThisTask)
    printf("[%s] Computing future LPT displacements  -- %s\n",fdate());

  for (ia=1; ia<=3; ia++)
    {
      if (!ThisTask)
	printf("[%s] Computing 1st derivative of ZEL source: %d\n",fdate(),ia);

      write_in_cvector(0, kdensity[0]);

      if (compute_derivative(0,ia,0))
	return 1;

      if (after)
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal[_x_]) * (local_y + local_z * MyGrids[0].GSlocal[_y_]);
		  products[index].Vel_after[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal[_x_] + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal[_y_]));
		}
	  if (params.WriteVmax)
	    if (write_product(13+ia,tag))
	      return 1;
	}
      else
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal[_x_]) * (local_y + local_z * MyGrids[0].GSlocal[_y_]);
		  products[index].Vel[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal[_x_] + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal[_y_]));
		}
	  if (params.WriteVmax)
	    if (write_product(ia,tag))
	      return 1;
	}
    }

#ifdef TWO_LPT
  /* 2LPT displacements */

#ifdef SCALE_DEPENDENT
  ScaleDep.order=2;
#endif

  for (ia=1; ia<=3; ia++)
    {
      if (!ThisTask)
	printf("[%s] Computing 1st derivative of 2LPT source: %d\n",fdate(),ia);

      write_in_cvector(0, kvector_2LPT);

      if (compute_derivative(0,ia,0))
	return 1;

      if (after)
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal[_x_]) * (local_y + local_z * MyGrids[0].GSlocal[_y_]);
		  products[index].Vel_2LPT_after[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal[_x_] + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal[_y_]));
		}
	  if (params.WriteVmax)
	    if (write_product(16+ia,tag))
	      return 1;
	}
      else
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal[_x_]) * (local_y + local_z * MyGrids[0].GSlocal[_y_]);
		  products[index].Vel_2LPT[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal[_x_] + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal[_y_]));
		}
	  if (params.WriteVmax)
	    if (write_product(4+ia,tag))
	      return 1;
	}
    }
#endif

  cputmp=MPI_Wtime()-cputmp;
  if (!ThisTask)
    printf("[%s] Done computing velocities, cpu time = %f s\n",fdate(),cputmp);
  cputime.vel+=cputmp;

  return 0;
}


int shift_all_displacements()
{
  /* Shifts all Vel_after to Vel, at all orders */

  int i, ia;
  /* This shifts Vel_after to Vel in the fft space */
  for (i=0; i<MyGrids[0].total_local_size; i++)
    for (ia=0; ia<3; ia++)
      {
	products[i].Vel[ia]=products[i].Vel_after[ia];
#ifdef TWO_LPT
	products[i].Vel_2LPT[ia]=products[i].Vel_2LPT_after[ia];
#ifdef THREE_LPT
	products[i].Vel_3LPT_1[ia]=products[i].Vel_3LPT_1_after[ia];
	products[i].Vel_3LPT_2[ia]=products[i].Vel_3LPT_2_after[ia];
#endif
#endif
      }
  return 0;
}


int recompute_group_velocities()
{
  /* Recompute average displacements of group velocities */
  int next,npart,i,ia;
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point > 0)
      {
	for (ia=0; ia<3; ia++)
	  groups[i].Vel[ia] = groups[i].Vel_after[ia] =
#ifdef TWO_LPT
	    groups[i].Vel_2LPT[ia] = groups[i].Vel_2LPT_after[ia] =
#ifdef THREE_LPT
	    groups[i].Vel_3LPT_1[ia] = groups[i].Vel_3LPT_1_after[ia] =
	    groups[i].Vel_3LPT_2[ia] = groups[i].Vel_3LPT_2_after[ia] =
#endif
#endif
	    0.0;

	next=groups[i].point;
	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    for (ia=0; ia<3; ia++)
	      {
		groups[i].Vel[ia]+=frag[next].Vel[ia];
		groups[i].Vel_after[ia]+=frag[next].Vel_after[ia];
#ifdef TWO_LPT
		groups[i].Vel_2LPT[ia]+=frag[next].Vel_2LPT[ia];
		groups[i].Vel_2LPT_after[ia]+=frag[next].Vel_2LPT_after[ia];
#ifdef THREE_LPT
		groups[i].Vel_3LPT_1[ia]+=frag[next].Vel_3LPT_1[ia];
		groups[i].Vel_3LPT_1_after[ia]+=frag[next].Vel_3LPT_1_after[ia];
		groups[i].Vel_3LPT_2[ia]+=frag[next].Vel_3LPT_2[ia];
		groups[i].Vel_3LPT_2_after[ia]+=frag[next].Vel_3LPT_2_after[ia];
#endif
#endif
	      }
	    next=linking_list[next];
	  }
	for (ia=0; ia<3; ia++)
	  {
	    groups[i].Vel[ia] /= (double)groups[i].Mass;
	    groups[i].Vel_after[ia] /= (double)groups[i].Mass;
#ifdef TWO_LPT
	    groups[i].Vel_2LPT[ia] /= (double)groups[i].Mass;
	    groups[i].Vel_2LPT_after[ia] /= (double)groups[i].Mass;
#ifdef THREE_LPT
	    groups[i].Vel_3LPT_1[ia] /= (double)groups[i].Mass;
	    groups[i].Vel_3LPT_1_after[ia] /= (double)groups[i].Mass;
	    groups[i].Vel_3LPT_2[ia] /= (double)groups[i].Mass;
	    groups[i].Vel_3LPT_2_after[ia] /= (double)groups[i].Mass;
#endif
#endif
	  }
      }

  return 0;
}


#endif


double myz;
gsl_function Function;

double Integrand_MF(double logm, void *param)
{
  double m=exp(logm);
  return m * AnalyticMassFunction(m,myz);
}

double compute_Nhalos_in_PLC(double z1, double z2)
{
  /* analytic prediction of the number of halos in the PLC */

  double MinMass    = log(params.ParticleMass*params.MinHaloMass);
  double delta_z    = 0.01;
  double number     = 0;
  double solidangle = (1-cos( (params.PLCAperture>90. ? 90. : params.PLCAperture)  /180.*PI) )*2.*PI;
  
  double result, error, upper, lower;
  gsl_function Function;
  
  Function.function = &Integrand_MF;
  lower=z1;
  do
    {
      upper=lower+delta_z;
      if (upper>z2)
	upper=z2;
      myz = 0.5*(upper+lower);

      gsl_integration_qags(&Function, MinMass, 37.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);

      number += result * solidangle * (pow(ComovingDistance(upper),3.) -
				       pow(ComovingDistance(lower),3.)) /3.;
      lower+=delta_z;
    } 
  while (upper<z2);

  return number;
}

int size2Mb(double *s)
{
  *s /= MBYTE;
  int gb=0;
  if (*s>1024.)
    {
      *s /= 1024.;
      gb=1;
    }
  return gb;
}

int estimate_file_size(void)
{

  /* estimates the size of output files */

  /* only Task 0 needs to work here */
  if (ThisTask)
    return 0;

  /* this works for binary output */
  if (params.CatalogInAscii)
    return 0;

  double       result, error, total=0.0, size, number;
  int          gb;
  gsl_function Function;
  
  printf("ESTIMATED STORAGE REQUIREMENTS:\n");

  double MinMass=log(params.ParticleMass*params.MinHaloMass);
  Function.function = &Integrand_MF;
  for (int iout=0; iout<outputs.n; iout++)
    {

      myz=outputs.z[iout];
      gsl_integration_qags(&Function, MinMass, 37.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
      number = result * pow(params.BoxSize_htrue,3.);
      size = number * sizeof(catalog_data);
      total+=size;
      gb=size2Mb(&size);
      printf("catalog, z=%6.4f, number of halos: %d, size: %f %s",outputs.z[iout],(int)number, size,(gb?"Gbyte":"Mbyte"));
      if (params.NumFiles>1)
	printf(" - each file will have a size of %f %s",size/(double)params.NumFiles,(gb?"Gbyte":"Mbyte"));
      printf("\n");
    }

  size = number * sizeof(catalog_data) * 1.4;
  total += size;
  gb=size2Mb(&size);
  printf("order-of-magnitude size of histories file: %f %s",1.4*size,(gb?"Gbyte":"Mbyte"));
  if (params.NumFiles>1)
    printf(" - each file will have a size of %f %s",size/(double)params.NumFiles,(gb?"Gbyte":"Mbyte"));
  printf("\n");

#ifdef PLC
  
  if (params.StartingzForPLC>0.)
    {
      number = compute_Nhalos_in_PLC(params.LastzForPLC, params.StartingzForPLC);
      size = number * sizeof(plc_write_data);
      total+=size;
      gb=size2Mb(&size);
      printf("past light cone, number of halos: %d, size: %f %s",(int)number, size,(gb?"Gbyte":"Mbyte"));
      if (params.NumFiles>1)
	printf(" - each file will have a size of %f %s",size/(double)params.NumFiles,(gb?"Gbyte":"Mbyte"));
      printf("\n");

    }
#endif

#ifdef LONGIDS
  double IDsize = (double)MyGrids[0].Ntotal * 8.;
#else
  double IDsize = (double)MyGrids[0].Ntotal * 4.;
#endif

  int nblo=3;
#ifndef TWO_LPT
  int nvel=3;
  nblo+=1;
#else
#ifndef THREE_LPT
  int nvel=6;
  nblo+=2;
#else
  int nvel=12;
  nblo+=4;
#endif
#endif

#ifdef ADD_RMAX_TO_SNAPSHOT
  nblo+=1;
#endif


#ifdef RECOMPUTE_DENSITY
  implementare...;
#endif

  if (params.WriteDensity)
    {
      size=268. + IDsize + 6. + (double)MyGrids[0].Ntotal * sizeof(float) + 6. + 2*40 + 6.;
      total+=size;
      gb=size2Mb(&size);
      printf("density snapshot size: %f %s", size, (gb?"Gbyte":"Mbyte"));
      if (params.NumFiles>1)
	printf(" - each file will have a size of %f %s",size/(double)params.NumFiles,(gb?"Gbyte":"Mbyte"));
      printf("\n");
    }

  if (params.WriteTimelessSnapshot)
    {
      size=268. + IDsize + 6. + (nvel+2) * ((double)MyGrids[0].Ntotal * sizeof(float) + 6.) + nblo*40 + 6.;
#ifdef ADD_RMAX_TO_SNAPSHOT
      size+=((double)MyGrids[0].Ntotal * sizeof(float) + 6.);
#endif
      total+=size;
      gb=size2Mb(&size);
      printf("timeless snapshot size: %f %s", size, (gb?"Gbyte":"Mbyte"));
      if (params.NumFiles>1)
	printf(" - each file will have a size of %f %s",size/(double)params.NumFiles,(gb?"Gbyte":"Mbyte"));
      printf("\n");
    }

  gb=size2Mb(&total);

  printf("Total storage: %f %s\n",total,(gb?"Gbyte":"Mbyte"));

  return 0;
}


