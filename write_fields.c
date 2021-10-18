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
#include <sys/types.h>
#include <sys/stat.h>

#define FORTRAN
#define ASCII_WRITE //Cubic_Galileon

int check_DataDir(void);

int MyGrid;

int write_fields(void)
{

  char tag[1];
  tag[0]='\0';
  MyGrid=0;

  if (params.WriteFmax)
    if (write_product(0,tag))
      return 1;

  if (params.WriteRmax)
    if (write_product(4,tag))
      return 1;

#ifndef RECOMPUTE_DISPLACEMENTS
  if (params.WriteVmax)
    {
      if (write_product(1,tag))
	return 1;
      if (write_product(2,tag))
	return 1;
      if (write_product(3,tag))
	return 1;
#ifdef TWO_LPT
      if (write_product(5,tag))
	return 1;
      if (write_product(6,tag))
	return 1;
      if (write_product(7,tag))
	return 1;
#ifdef THREE_LPT
      if (write_product(8,tag))
	return 1;
      if (write_product(9,tag))
	return 1;
      if (write_product(10,tag))
	return 1;
      if (write_product(11,tag))
	return 1;
      if (write_product(12,tag))
	return 1;
      if (write_product(13,tag))
	return 1;
#endif
#endif
    }
#endif
  return 0;

}

int write_density(int ThisGrid)
{
  char tag[1];
  tag[0]='\0';
  MyGrid=ThisGrid;
  if (write_product(-1,tag))
    return 1;
  return 0;
}


int write_product(int flag, char *tag)
{
  int ix, iy, local_z, global_z;
  size_t planesize;
  int *belongs_local, *belongs;
  char filename[LBLENGTH];
  PRODFLOAT *plane;
  FILE *file;
#if defined(FORTRAN) && !defined(ASCII_WRITE)
  int idummy;
#endif
#ifdef ASCII_WRITE
  int last;
#endif
  MPI_Status status;

  /* NB: this must be changed for pencil-distributed memory */

  if (check_DataDir())
    return 1;

  /* allocation of plane for communication */
  planesize=MyGrids[MyGrid].GSglobal_x*MyGrids[MyGrid].GSglobal_y;
  if ( (plane=(PRODFLOAT*)malloc(planesize*sizeof(PRODFLOAT))) == 0x0 )
    {
      printf("ERROR on task %d: could not allocate plane for writing product %d\n",ThisTask,flag);
      fflush(stdout);
      return 1;
    }

  /* map of tasks -> planes */
  belongs        = (int*)malloc(MyGrids[MyGrid].GSglobal_z * sizeof(int));
  belongs_local  = (int*)malloc(MyGrids[MyGrid].GSglobal_z * sizeof(int));

  if (belongs==0x0 || belongs_local==0x0)
    {
      printf("ERROR on task %d: could not allocate some memory\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  for (global_z=0; global_z<MyGrids[MyGrid].GSglobal_z; global_z++)
    {
      belongs[global_z]=0;

      if (global_z >= MyGrids[MyGrid].GSstart_z && global_z < MyGrids[MyGrid].GSlocal_z + MyGrids[MyGrid].GSstart_z)
	belongs_local[global_z]=ThisTask;
      else
	belongs_local[global_z]=0;
    }

  if (MPI_Reduce(belongs_local, belongs, 
		 MyGrids[MyGrid].GSglobal_z, MPI_INT, MPI_SUM, 0, 
		 MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: MPI_Reduce failed in write_product\n",ThisTask);
      fflush(stdout);
      return 1;
    }
  
  if (MPI_Bcast(belongs, MyGrids[MyGrid].GSglobal_z, 
		MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: MPI_Bcast failed in write_product\n",ThisTask);
      fflush(stdout);
      return 1;
    }


  /* file opening and header */
  if (!ThisTask)
    {
      /* Task 0 opens the file */
      switch (flag)
	{
	case -1:
	  sprintf(filename,"Data/density.grid%d%s.%s.d",MyGrid,tag,params.RunFlag);
	  break;
	case 0:
	  sprintf(filename,"Data/Fmax%s.%s.d",tag,params.RunFlag);
	  break;
	case 1:
	  sprintf(filename,"Data/Velx%s.%s.d",tag,params.RunFlag);
	  break;
	case 2:
	  sprintf(filename,"Data/Vely%s.%s.d",tag,params.RunFlag);
	  break;
	case 3:
	  sprintf(filename,"Data/Velz%s.%s.d",tag,params.RunFlag);
	  break;
	case 4:
	  sprintf(filename,"Data/Rmax%s.%s.d",tag,params.RunFlag);
	  break;
#ifdef TWO_LPT
	case 5:
	  sprintf(filename,"Data/Velx_2LPT%s.%s.d",tag,params.RunFlag);
	  break;
	case 6:
	  sprintf(filename,"Data/Vely_2LPT%s.%s.d",tag,params.RunFlag);
	  break;
	case 7:
	  sprintf(filename,"Data/Velz_2LPT%s.%s.d",tag,params.RunFlag);
	  break;
#ifdef THREE_LPT
	case 8:
	  sprintf(filename,"Data/Velx_3LPT_1%s.%s.d",tag,params.RunFlag);
	  break;
	case 9:
	  sprintf(filename,"Data/Vely_3LPT_1%s.%s.d",tag,params.RunFlag);
	  break;
	case 10:
	  sprintf(filename,"Data/Velz_3LPT_1%s.%s.d",tag,params.RunFlag);
	  break;
	case 11:
	  sprintf(filename,"Data/Velx_3LPT_2%s.%s.d",tag,params.RunFlag);
	  break;
	case 12:
	  sprintf(filename,"Data/Vely_3LPT_2%s.%s.d",tag,params.RunFlag);
	  break;
	case 13:
	  sprintf(filename,"Data/Velz_3LPT_2%s.%s.d",tag,params.RunFlag);
	  break;
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
	case 14:
	  sprintf(filename,"Data/Velx_after%s.%s.d",tag,params.RunFlag);
	  break;
	case 15:
	  sprintf(filename,"Data/Vely_after%s.%s.d",tag,params.RunFlag);
	  break;
	case 16:
	  sprintf(filename,"Data/Velz_after%s.%s.d",tag,params.RunFlag);
	  break;
#ifdef TWO_LPT
	case 17:
	  sprintf(filename,"Data/Velx_2LPT_after%s.%s.d",tag,params.RunFlag);
	  break;
	case 18:
	  sprintf(filename,"Data/Vely_2LPT_after%s.%s.d",tag,params.RunFlag);
	  break;
	case 19:
	  sprintf(filename,"Data/Velz_2LPT_after%s.%s.d",tag,params.RunFlag);
	  break;
#endif
#endif
	default:
	  return 1;
	  break;
	}

      if ( (file=fopen(filename,"w")) == 0x0 )
	{
	  printf("ERROR on task %d: cannot open file %s\n",ThisTask, filename);
	  fflush(stdout);
	  return 1;
	}

      printf("Writing field in file %s\n",filename);

      /* Task 0 writes the header */
#ifdef ASCII_WRITE
      fprintf(file,"%d %d %d\n",MyGrids[MyGrid].GSglobal_x,MyGrids[MyGrid].GSglobal_y,MyGrids[MyGrid].GSglobal_z);
#else
#ifdef FORTRAN
      idummy=3*sizeof(int);
      fwrite(&idummy,sizeof(int),1,file);
#endif
      fwrite(&MyGrids[MyGrid].GSglobal_x,sizeof(int),1,file);
      fwrite(&MyGrids[MyGrid].GSglobal_y,sizeof(int),1,file);
      fwrite(&MyGrids[MyGrid].GSglobal_z,sizeof(int),1,file);
#ifdef FORTRAN
      fwrite(&idummy,sizeof(int),1,file);
#endif
#endif
    }

  /* loop on planes */
  for (global_z=0; global_z<MyGrids[MyGrid].GSglobal_z; global_z++)
    {
      local_z=global_z-MyGrids[MyGrid].GSstart_z;
      
      if (!ThisTask)
	{
	  /* Task 0 writes global_z on the output file */
#ifdef ASCII_WRITE
	  fprintf(file,"%d\n",global_z);
#else
#ifdef FORTRAN
	  idummy=sizeof(int);
	  fwrite(&idummy,sizeof(int),1,file);
#endif
	  fwrite(&global_z,sizeof(int),1,file);
#ifdef FORTRAN
	  fwrite(&idummy,sizeof(int),1,file);
#endif
#endif
	}

      /* the task that possesses the plane writes the wanted product on the plane vector */
      if (ThisTask==belongs[global_z])
	{
	  switch(flag)
	    {
	    case -1:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=density[MyGrid][ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)];
	      break;
	    case 0:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Fmax;
	      break;
	    case 1:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel[0];
	      break;
	    case 2:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel[1];
	      break;
	    case 3:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel[2];
	      break;
	    case 4:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=(float)products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Rmax;
	      break;
#ifdef TWO_LPT
	    case 5:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT[0];
	      break;
	    case 6:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT[1];
	      break;
	    case 7:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT[2];
	      break;
#ifdef THREE_LPT
	    case 8:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_1[0];
	      break;
	    case 9:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_1[1];
	      break;
	    case 10:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_1[2];
	      break;
	    case 11:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_2[0];
	      break;
	    case 12:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_2[1];
	      break;
	    case 13:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_3LPT_2[2];
	      break;
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
	    case 14:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_after[0];
	      break;
	    case 15:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_after[1];
	      break;
	    case 16:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_after[2];
	      break;
#ifdef TWO_LPT
	    case 17:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT_after[0];
	      break;
	    case 18:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT_after[1];
	      break;
	    case 19:
	      for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
		for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
		  plane[ix + iy*MyGrids[MyGrid].GSglobal_x]=products[ix + (MyGrids[MyGrid].GSlocal_x) * (iy + local_z * MyGrids[MyGrid].GSlocal_y)].Vel_2LPT_after[2];
	      break;
#endif
#endif
	    default:
	      return 1;
	      break;
	    }
	}

      /* if the plane does not belong to task 0 then it must be communicated */
      if (belongs[global_z])
	{
	  if (ThisTask==belongs[global_z])
	    {
	      /* Task belongs[iz] sends the plane to task 0 */
	      MPI_Send(plane, planesize, MPI_PRODFLOAT, 0, 0, MPI_COMM_WORLD);
	      
	    }
	  else if (!ThisTask)
	    {
	      /* task 0 receives the info from task belongs[global_z] and writes it on the file */
	      MPI_Recv(plane, planesize, MPI_PRODFLOAT, belongs[global_z], 0, MPI_COMM_WORLD, &status);
	    }
	}

      if (!ThisTask)
	{

#ifdef ASCII_WRITE
	  /* writes the plane */
	  for (iy=0; iy<MyGrids[MyGrid].GSglobal_y; iy++)
	    for (ix=0; ix<MyGrids[MyGrid].GSglobal_x; ix++)
	      {
		fprintf(file," %13.6g",plane[ix + iy*MyGrids[MyGrid].GSglobal_x]);
		last=1;
		if (!((1+ix + iy*MyGrids[MyGrid].GSglobal_x)%8))
		  {
 		    fprintf(file,"\n");
		    last=0;
		  }
	      }
	  if (last)
	    fprintf(file,"\n");
#else
#ifdef FORTRAN
	  idummy=planesize*sizeof(PRODFLOAT);
	  fwrite(&idummy,sizeof(int),1,file);
#endif
	  fwrite(plane,sizeof(PRODFLOAT),planesize,file);
#ifdef FORTRAN
	  fwrite(&idummy,sizeof(int),1,file);
#endif
#endif
	}
    }

  MPI_Barrier(MPI_COMM_WORLD);

  if (!ThisTask)
    {
      fclose(file);
      printf("Written file %s\n",filename);
    }

  free(belongs_local);
  free(belongs);
  free(plane);

  return 0;
}

int check_DataDir()
{
  struct stat dr;

  if (!ThisTask)
    {
      if(stat(params.DataDir,&dr))
	{
	  printf("Creating directory %s\n",params.DataDir);
	  if (mkdir(params.DataDir,0755))
	    {
	      printf("ERROR IN CREATING DIRECTORY %s (task 0)\n",params.DataDir);
	      fflush(stdout);
	      return 1;
	    }
	}
    }

  return 0;
}
