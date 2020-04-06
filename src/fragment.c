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
#ifdef RECOMPUTE_DISPLACEMENTS
int compute_future_LPT_displacements(double, int);
int init_segmentation();
int shift_all_displacements();
int recompute_group_velocities();
#endif

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

/* NB the list of fragment parameters option has been removed, 
   it should be implemented by repeated calls of fragment().
   One should then check what is to be reinitialized after each call. */

int fragment()
{
  /* This function is the driver for fragmentation of collapsed medium 
     and construction of halo catalogs */

  int Npeaks;
  double tmp;

  /* timing */

  if (!ThisTask)
    printf("\n[%s] Second part: fragmentation of the collapsed medium\n",fdate());


#ifdef RECOMPUTE_DISPLACEMENTS
  int mysegment;
  if (init_segmentation())
    return 1;

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

  cputime.frag=MPI_Wtime();

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

      if (init_fragment_2(&Npeaks))
	return 1;


      /* Generates the group catalogue */
      /* The call is segmented into consecutive calls */
#ifdef RECOMPUTE_DISPLACEMENTS
      for (mysegment=0; mysegment<Segment.n; mysegment++)
	{
	  Segment.mine=mysegment;

	  /* The first segment uses only the first computed velocity
	     The second segment has the two velocities already loaded */
	  if (mysegment>=2)
	    {
	      /* computation of displacements for the next redshift */
	      if (!ThisTask)
		printf("\n[%s] Computing displacements for redshift %f\n",fdate(),Segment.z[mysegment]);
	      if (shift_all_displacements())
		return 1;

	      if (compute_future_LPT_displacements(Segment.z[mysegment],(mysegment>0)))
		return 1;

	      tmp=MPI_Wtime();
	      if (!ThisTask)
		printf("[%s] Starting re-distribution of products\n",fdate());

	      if (distribute())
		return 1;

	      tmp=MPI_Wtime()-tmp;
	      cputime.distr += tmp;
	      if (!ThisTask)
		printf("[%s] Re-distribution of Fmax products done, cputime = %14.6f\n",fdate(),tmp);

	      if (recompute_group_velocities())
		return 1;
	    }

	  tmp=MPI_Wtime();
      
	  if (build_groups(Npeaks,Segment.z[mysegment]))
	    return 1;

	  cputime.group+=MPI_Wtime()-tmp;	  
	}
#else
      tmp=MPI_Wtime();

      if (build_groups(Npeaks,outputs.zlast))
	return 1;

      cputime.group+=MPI_Wtime()-tmp;
#endif


      /* next slice */
      subbox.mybox_z += subbox.nbox_z_thisslice;
      subbox.start_z = find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);
      subbox.Lgrid_z = find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);
      subbox.stabl_z = subbox.start_z - subbox.safe_z;

    }
  /* cpu time for fragmentation */
  cputime.frag = MPI_Wtime() - cputime.frag;

  if (!ThisTask)
    printf("[%s] Finishing fragment, total cputime = %14.6f\n",fdate(),cputime.frag);

  return 0;
}


static int index_compare(const void * a,const void * b)
{
  if ( frag[*((int *)a)].Fmax == frag[*((int *)b)].Fmax )
    return 0;
  else
    if ( frag[*((int *)a)].Fmax > frag[*((int *)b)].Fmax )
      return -1;
    else
      return 1;
}


int init_fragment_1(void)
{
  /* Initializes all quantities for the run */
  set_fragment_parameters(ORDER_FOR_GROUPS);

  /* reallocates the memory, for the part needed to redistribution */
  if (reallocate_memory_for_fragmentation_1())
    return 1;

  return 0;
}


int init_fragment_2(int *Npeaks)
{
  int Ngood, tot_good;
  double tmp;

  tmp=MPI_Wtime();

  /* this routine redistributes the products from planes to sub-boxes */
  if (!ThisTask)
    printf("[%s] Starting re-distribution of products\n",fdate());
  if (distribute())
    return 1;

  tmp=MPI_Wtime()-tmp;
  cputime.distr += tmp;
  if (!ThisTask)
    printf("[%s] Re-distribution of Fmax products done, cputime = %14.6f\n",fdate(),tmp);

  /* counts the number of peaks in the sub-volume to know the number of groups */
  *Npeaks=count_peaks(&Ngood);

  if (MPI_Reduce(&Ngood, &tot_good, 1, 
		 MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: init_fragment failed an MPI_Reduce\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  if (!ThisTask)
    printf("[%s] Task 0 found %d peaks, in the well resolved region: %d. Total number of peaks: %d\n",
	   fdate(),*Npeaks,Ngood,tot_good);

  if (reallocate_memory_for_fragmentation_2(*Npeaks))
    return 1;

  tmp=MPI_Wtime();
  if (!ThisTask)
    printf("[%s] Starting sorting\n",fdate());

  qsort((void *)indices, subbox.Npart, sizeof(int), index_compare);

  tmp=MPI_Wtime()-tmp;
  cputime.sort+=tmp;
  if (!ThisTask)
    printf("[%s] Sorting done, total cputime = %14.6f\n",fdate(),tmp);

  return 0;
}


int count_peaks(int *ngood)
{

  /* counts the number of Fmax peaks down to the final redshift
     to know the total number of groups */

  int iz,i,j,k,kk,nn,ngroups_tot,Lgridxy,peak_cond,i1,j1,k1;

  ngroups_tot=0;     /* number of groups, group 1 is the filament group */
  *ngood=0;          /* number of groups out of the safety boundary */
  Lgridxy = subbox.Lgwbl_x * subbox.Lgwbl_y;

  for (iz=0; iz<subbox.Npart; iz++)
    {

      /* position on the local box */
      k=iz/Lgridxy;
      kk=iz-k*Lgridxy;
      j=kk/subbox.Lgwbl_x;
      i=kk-j*subbox.Lgwbl_x;

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

	  peak_cond &= (frag[iz].Fmax > frag[i1 +j1*subbox.Lgwbl_x + k1*Lgridxy].Fmax);

	}

      if (peak_cond && frag[iz].Fmax >= outputs.Flast)
	{
	  ngroups_tot++;
	  if ( i>=subbox.safe_x && i<subbox.Lgwbl_x-subbox.safe_x &&
	       j>=subbox.safe_y && j<subbox.Lgwbl_y-subbox.safe_y &&
	       k>=subbox.safe_z && k<subbox.Lgwbl_z-subbox.safe_z)
	    (*ngood)++;

	}
    }

  return ngroups_tot;
}


#ifdef RECOMPUTE_DISPLACEMENTS
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
    printf("[%s] Computing future LPT displacements  -- %s\n",fdate(), tag); // levare tag

  for (ia=1; ia<=3; ia++)
    {
      if (!ThisTask)
	printf("[%s] Computing 1st derivative of ZEL source: %d\n",fdate(),ia);

      write_in_cvector(0, kdensity[0]);

      if (compute_derivative(0,ia,0))
	return 1;

      if (after)
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal_x) * (local_y + local_z * MyGrids[0].GSlocal_y);
		  products[index].Vel_after[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal_y));
		}
	  if (params.WriteVmax)
	    if (write_product(13+ia,tag))
	      return 1;
	}
      else
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal_x) * (local_y + local_z * MyGrids[0].GSlocal_y);
		  products[index].Vel[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal_y));
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
	  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal_x) * (local_y + local_z * MyGrids[0].GSlocal_y);
		  products[index].Vel_2LPT_after[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal_y));
		}
	  if (params.WriteVmax)
	    if (write_product(16+ia,tag))
	      return 1;
	}
      else
	{
	  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
	    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
	      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
		{
		  index = local_x + (MyGrids[0].GSlocal_x) * (local_y + local_z * MyGrids[0].GSlocal_y);
		  products[index].Vel_2LPT[ia-1] = 
		    *(rvector_fft[0] + local_x + (MyGrids[0].GSlocal_x + MyGrids[0].off) * (local_y + local_z * MyGrids[0].GSlocal_y));
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


int init_segmentation(void)
{

  /* here the segmentation of the fragmentation process is defined */
  /* for now, segmentation is taken from the outputs file */
  Segment.n=outputs.n;
  for (int i=0; i<outputs.n; i++)
    Segment.z[i]=outputs.z[i];
  Segment.mine=1;
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

// QUESTO DOVREBBE ESSERE INUTILE
/*   /\* This shifts Vel_after to Vel in the subbox space *\/ */
/*   for (i=0; i<subbox.Npart; i++) */
/*     for (ia=0; ia<3; ia++) */
/*       frag[i].Vel[ia]=frag[i].Vel_after[ia]; */
/* #ifdef TWO_LPT */
/*   frag[i].Vel_2LPT[ia]=frag[i].Vel_2LPT_after[ia]; */
/* #ifdef THREE_LPT */
/*   frag[i].Vel_3LPT_1[ia]=frag[i].Vel_3LPT_1_after[ia]; */
/*   frag[i].Vel_3LPT_2[ia]=frag[i].Vel_3LPT_2_after[ia]; */
/* #endif */
/* #endif */


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
