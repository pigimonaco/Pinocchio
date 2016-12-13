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

//#define READ_FRAGMENT_PARAMS_FROM_FILE

#if defined(LIST_OF_FRAGMENT_PARAMETERS) && defined(READ_FRAGMENT_PARAMS_FROM_FILE)
#error Trying to compile with LIST_OF_FRAGMENT_PARAMETERS and READ_FRAGMENT_PARAMS_FROM_FILE
#endif


int ngroups;

#ifdef LIST_OF_FRAGMENT_PARAMETERS
  static int first_time=1;
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



int fragment()
{
  int Npeaks;
  double tmp;

  /* timing */
  cputime.frag=MPI_Wtime();

  if (!ThisTask)
    printf("[%s] Second part: fragmentation of the collapsed medium\n",fdate());

  if (init_fragment_1())
    return 1;

#ifdef LIST_OF_FRAGMENT_PARAMETERS
  static FILE *listfile;
  int nn, more_params;
  double tobcast[6];
  char buffer[BLENGTH],MyRunFlag[BLENGTH],check[BLENGTH];
  static int save_mybox_z;
  struct stat dr;

  save_mybox_z = subbox.mybox_z;
  strcpy(MyRunFlag,params.RunFlag);
  do
    {
      /* this is the loop on fragmentation parameters */

      if (!ThisTask)
	{
	  more_params=0;
	  tobcast[0]=-99.0;
	  if (first_time)
	    {
	      listfile=fopen("list_of_fragment_parameters.txt","r");
	      if (listfile==0x0)
		printf("\n[%s] list_of_fragment_parameters.txt not found, using standard parameters\n",fdate());
	      else
		printf("\n[%s] list_of_fragment_parameters.txt found\n",fdate());

	    }
	  if (listfile!=0x0)
	    {
	      while (!more_params && !feof(listfile))
		{
		  fgets(buffer,BLENGTH,listfile);
		  if (!feof(listfile) && buffer[0]!='#' && buffer[0]!='%' &&
		      (nn=sscanf(buffer,"%lf %lf %lf %lf %lf %lf",&f_m,&f_rm,&espo,&f_a,&f_ra,&f_200)) == 6)
		    {
		      sprintf(params.RunFlag,"%s.fm%5.3f-rm%5.3f-em%5.3f-fa%5.3f-ra%5.3f-ea%5.3f",MyRunFlag,f_m,f_rm,espo,f_a,f_ra,f_200);
		      sprintf(check,"pinocchio.%6.4f.%s.mf.out",outputs.z[outputs.n-1],params.RunFlag);
		      if (stat(check,&dr))
			{
			  more_params=1;
			  tobcast[0]=f_m;
			  tobcast[1]=f_rm;
			  tobcast[2]=espo;
			  tobcast[3]=f_a;
			  tobcast[4]=f_ra;
			  tobcast[5]=f_200;
			}
		      else
			printf("This combination of parameters has already been processed: %s\n",params.RunFlag);
		    }
		}
	    }
	}
      MPI_Bcast(tobcast, 6, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (ThisTask)
	{
	  more_params=(tobcast[0]>=0.0);
	  if (more_params)
	    {
	      f_m   =tobcast[0];
	      f_rm  =tobcast[1];
	      espo  =tobcast[2];
	      f_a   =tobcast[3];
	      f_ra  =tobcast[4];
	      f_200 =tobcast[5];
	    }
	}

      if (first_time || more_params)
	{

	  if (!ThisTask)
	    {
	      printf("\n");
	      printf("Fragmentation parameters for this fragmentation run:\n");
	      printf("f_m     = %f\n",f_m);
	      printf("f_rm    = %f\n",f_rm);
	      printf("espo    = %f\n",espo);
	      printf("f_a     = %f\n",f_a);
	      printf("f_ra    = %f\n",f_ra);
	      printf("f_200   = %f\n",f_200);
	      printf("variance on the grid: %f\n",Smoothing.TrueVariance[Smoothing.Nsmooth-1]);
	      printf("\n");
	    }
#endif

	  cputime.group=cputime.sort=cputime.distr=0.0;
#ifdef PLC
	  cputime.plc=0.0;
#endif

#ifdef LIST_OF_FRAGMENT_PARAMETERS
	  subbox.mybox_z = save_mybox_z;
	  subbox.start_z = find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);
	  subbox.Lgrid_z = find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);
	  subbox.stabl_z = subbox.start_z - subbox.safe_z;
#endif

	  for (ThisSlice=0; ThisSlice<NSlices; ThisSlice++)
	    {
	      if (!ThisTask)
		printf("\n[%s] Starting slice N. %d of %d\n",fdate(),ThisSlice+1,NSlices);

	      if (init_fragment_2(&Npeaks))
		return 1;

	      tmp=MPI_Wtime();

	      /* Generates the group catalogue */
	      if (build_groups(Npeaks))
		return 1;

	      cputime.group+=MPI_Wtime()-tmp;

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

#ifdef LIST_OF_FRAGMENT_PARAMETERS
	}
      first_time=0;
    }
  while (more_params);
#endif
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

#ifndef LIST_OF_FRAGMENT_PARAMETERS
  /* writes used parameters */
  if (!ThisTask)
    {
      printf("Fragmentation parameters:\n");
      printf("f      = %f\n",f_m);
      printf("f_rm   = %f\n",f_rm);
      printf("espo   = %f\n",espo);
      //printf("f_a    = %f\n",f_a);
      printf("f_ra   = %f\n",f_ra);
      //printf("f_200  = %f\n",f_200);
      printf("variance on the grid: %f\n",Smoothing.TrueVariance[Smoothing.Nsmooth-1]);
      printf("\n");
    }
#endif

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
#ifdef LIST_OF_FRAGMENT_PARAMETERS
  if (first_time || NSlices>1)
    {
#endif
      if (!ThisTask)
	printf("[%s] Starting re-distribution of products\n",fdate());
      if (distribute())
	return 1;
#ifdef LIST_OF_FRAGMENT_PARAMETERS
    }
#endif

  tmp=MPI_Wtime()-tmp;
  cputime.distr += tmp;
  if (!ThisTask)
    printf("[%s] Re-distribution of Fmax done, cputime = %14.6f\n",fdate(),tmp);

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

  ngroups_tot=0;     // number of groups, group 1 is the filament group
  *ngood=0;          // number of groups out of the safety boundary
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
