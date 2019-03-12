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
#include <sys/types.h>
#include <sys/stat.h>

#define FORTRAN
//#define ASCII_WRITE

#define GRID MyGrids[ThisGrid]

int check_DataDir(void);
int write_product(int, int);
int write_fmax(int);
int MyGrid;

int write_fields(void)
{

  MyGrid=0;

  if (params.WriteFmax)
    write_fmax(0);
    /* if (write_product(0, MyGrid)) */
    /*   return 1; */

  if (params.WriteRmax)
    if (write_product(4, MyGrid))
      return 1;

  if (params.WriteVmax)
    {
      if (write_product(1, MyGrid))
	return 1;
      if (write_product(2, MyGrid))
	return 1;
      if (write_product(3, MyGrid))
	return 1;
#ifdef TWO_LPT
      if (write_product(5, MyGrid))
	return 1;
      if (write_product(6, MyGrid))
	return 1;
      if (write_product(7, MyGrid))
	return 1;
#ifdef THREE_LPT
      if (write_product(8, MyGrid))
	return 1;
      if (write_product(9, MyGrid))
	return 1;
      if (write_product(10, MyGrid))
	return 1;
      if (write_product(11, MyGrid))
	return 1;
      if (write_product(12, MyGrid))
	return 1;
      if (write_product(13, MyGrid))
	return 1;
#endif
#endif
    }
  return 0;

}

int write_density(int ThisGrid)
{
  /* MyGrid=ThisGrid; */
  /* if (write_product(-1)) */
  /*   return 1; */
  /* return 0; */

  int       T, i, j, k;
  unsigned  index, index_i, index_j;
  unsigned  lindex, lindex_i, lindex_j;
  char      filename[200];
  FILE     *file;


  sprintf(filename, "partial_density.%4d", ThisTask);

  for(T = 0; T < NTasks; T++)
    {
      if(ThisTask == T)
	{
	  file = fopen(filename, "w");

	  for(i = 0; i < GRID.GSlocal[_x_]; i++)
	    {
	      lindex_i = i * GRID.GSlocal[_y_];
	      index_i  = (i + GRID.GSstart[_x_]) * GRID.GSglobal[_y_];
	      
	      for(j = 0; j < GRID.GSlocal[_y_]; j++)
		{
		  lindex_j = (lindex_i + j) * GRID.GSlocal[_z_];
		  index_j  = (index_i + j + GRID.GSstart[_y_]) * GRID.GSglobal[_z_];
		  
		  for(k = 0; k < GRID.GSlocal[_z_]; k++)
		    {
		      lindex = lindex_j + k;
		      index = index_j + (k + GRID.GSstart[_z_]);
		      
 		      fprintf(file, "%u %g\n", index, density[ThisGrid][lindex]);
		    }
		}
	    }
	  
	  fclose(file);
	}
      else
	MPI_Barrier(MPI_COMM_WORLD);
    }

  if(ThisTask == 0)
    {
      int err;
      err = system("cat partial_density* | sort -k1 -n > my.density");
      err = system("rm -f partial_*");
      if(err != 0)
	dprintf(VXX, 0, "something got wrong while writing density files\n");
    }

  return 0 ;
}


int write_fmax(int ThisGrid)
{
  /* MyGrid=ThisGrid; */
  /* if (write_product(-1)) */
  /*   return 1; */
  /* return 0; */

  int       T, i, j, k;
  unsigned  index, index_i, index_j, index_reverse;
  unsigned  lindex, lindex_i, lindex_j;
  char      filename[200];
  FILE     *file;

  unsigned local[3], start[3], C[3];
  
  /* C[_x_] = _y_; */
  /* C[_y_] = _x_; */
  /* C[_z_] = _z_; */

  C[_x_] = _x_;
  C[_y_] = _y_;
  C[_z_] = _z_;
  
  for(i = 0; i < 3; i++)
    {
      local[i] = GRID.GSlocal[C[i]];
      start[i] = GRID.GSstart[C[i]];
    }


  sprintf(filename, "partial_fmax.%4d", ThisTask);

  for(T = 0; T < NTasks; T++)
    {
      if(ThisTask == T)
	{
	  file = fopen(filename, "w");

	  for(i = 0; i < GRID.GSlocal[_x_]; i++)
	    {
	      lindex_i = i * GRID.GSlocal[_y_];
	      index_i  = (i + GRID.GSstart[_x_]) * GRID.GSglobal[_y_];
	      
	      for(j = 0; j < GRID.GSlocal[_y_]; j++)
		{
		  lindex_j = (lindex_i + j) * GRID.GSlocal[_z_];
		  index_j  = (index_i + j + GRID.GSstart[_y_]) * GRID.GSglobal[_z_];
		  
		  for(k = 0; k < GRID.GSlocal[_z_]; k++)
		    {
		      lindex = lindex_j + k;
		      index  = index_j + k + GRID.GSstart[_z_];
		      
		      index_reverse = ((j+start[_y_]) * GRID.GSglobal[_x_] + (i + start[_x_])) * GRID.GSglobal[_z_] + k + start[_z_];

 		      fprintf(file, "%u %g\n", index_reverse, products[lindex].Fmax);
		    }
		}
	    }
	  
	  fclose(file);
	}
      else
	MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if(ThisTask == 0)
    {
      int err;
      err = system("cat partial_fmax* | sort -k1 -n > my.fmax");
      err = system("rm -f partial_*");
      if(err != 0)
	dprintf(VXX, 0, "something got wrong while writing density files\n");
    }

  return 0 ;
}


int write_product(int flag, int ThisGrid)
{
  int     C, ix, iy, local_z, global_z;
  size_t  cubesize;
  char    filename[BLENGTH];
  double *cube;
  FILE   *file;
#ifdef FORTRAN
  int     idummy;
#endif
#ifdef ASCII_WRITE
  int     last;
#endif
  MPI_Status status;

  /* NB: this must be changed for pencil-distributed memory */

  if (check_DataDir())
    return 1;

  /* allocation of plane for communication */
  int xsize, ysize, zsize;
  xsize = GRID.GSlocal[_x_];
  ysize = GRID.GSlocal[_y_];
  zsize = GRID.GSlocal[_z_];
  cubesize= xsize * ysize * zsize;
  
  while ( (cube=(double*)malloc(cubesize*sizeof(double))) == 0x0 )
    {
      if(xsize > 2)
	xsize *=  0.8;
      else if(ysize > 2)
	ysize *= 0.8;
      else if(zsize > 2)
	zsize *= 0.8;

      if(xsize * ysize * zsize == 1)
	{
	  dprintf(VXERR, ThisTask,
		  "ERROR on task %d: could not allocate memory"
		  " to communicate field (down to %d x %d x %d size!)\n",
		  ThisTask, xsize, ysize, zsize);
	  fflush(stdout);
	  return 1;
	}

      cubesize = xsize * ysize * zsize;
    }

  free(cube);
  
  /* file opening and header */
  if (!ThisTask)
    {
      /* Task 0 opens the file */
      switch (flag)
	{
	case -1:
	  sprintf(filename,"Data/density.grid%d.%s.d",MyGrid,params.RunFlag);
	  break;
	case 0:
	  sprintf(filename,"Data/Fmax.%s.d",params.RunFlag);
	  break;
	case 1:
	  sprintf(filename,"Data/Velx.%s.d",params.RunFlag);
	  break;
	case 2:
	  sprintf(filename,"Data/Vely.%s.d",params.RunFlag);
	  break;
	case 3:
	  sprintf(filename,"Data/Velz.%s.d",params.RunFlag);
	  break;
	case 4:
	  sprintf(filename,"Data/Rmax.%s.d",params.RunFlag);
	  break;
#ifdef TWO_LPT
	case 5:
	  sprintf(filename,"Data/Velx_2LPT.%s.d",params.RunFlag);
	  break;
	case 6:
	  sprintf(filename,"Data/Vely_2LPT.%s.d",params.RunFlag);
	  break;
	case 7:
	  sprintf(filename,"Data/Velz_2LPT.%s.d",params.RunFlag);
	  break;
#ifdef THREE_LPT
	case 8:
	  sprintf(filename,"Data/Velx_3LPT_1.%s.d",params.RunFlag);
	  break;
	case 9:
	  sprintf(filename,"Data/Vely_3LPT_1.%s.d",params.RunFlag);
	  break;
	case 10:
	  sprintf(filename,"Data/Velz_3LPT_1.%s.d",params.RunFlag);
	  break;
	case 11:
	  sprintf(filename,"Data/Velx_3LPT_2.%s.d",params.RunFlag);
	  break;
	case 12:
	  sprintf(filename,"Data/Vely_3LPT_2.%s.d",params.RunFlag);
	  break;
	case 13:
	  sprintf(filename,"Data/Velz_3LPT_2.%s.d",params.RunFlag);
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
      fprintf(file,"%ld %ld %ld\n", MyGrids[MyGrid].GSglobal[_x_], MyGrids[MyGrid].GSglobal[_y_], MyGrids[MyGrid].GSglobal[_z_]);
#else
#ifdef FORTRAN
      idummy=3*sizeof(int);
      fwrite(&idummy,sizeof(int),1,file);
#endif
      fwrite(&MyGrids[MyGrid].GSglobal[_x_],sizeof(int),1,file);
      fwrite(&MyGrids[MyGrid].GSglobal[_y_],sizeof(int),1,file);
      fwrite(&MyGrids[MyGrid].GSglobal[_z_],sizeof(int),1,file);
#ifdef FORTRAN
      fwrite(&idummy,sizeof(int),1,file);
#endif
#endif
    }

/*   /\* loop on cubes (i.e. tasks) *\/ */
/*   for(C = 0; C < NTasks; C++)     */
/*     { */
      
/*       if (!ThisTask) */
/* 	{ */
/* 	  /\* Task 0 writes global_z on the output file *\/ */
/* #ifdef ASCII_WRITE */
/* 	  fprintf(file,"%d\n",global_z); */
/* #else */
/* #ifdef FORTRAN */
/* 	  idummy=sizeof(int); */
/* 	  fwrite(&idummy,sizeof(int),1,file); */
/* #endif */
/* 	  fwrite(&global_z,sizeof(int),1,file); */
/* #ifdef FORTRAN */
/* 	  fwrite(&idummy,sizeof(int),1,file); */
/* #endif */
/* #endif */
/* 	} */

/*       the task that possesses the plane writes the wanted product on the plane vector */
/*       if (ThisTask==belongs[global_z]) */
/* 	{ */
/* 	  switch(flag) */
/* 	    { */
/* 	    case -1: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=density[MyGrid][ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])]; */
/* 	      break; */
/* 	    case 0: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Fmax; */
/* 	      break; */
/* 	    case 1: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax[0]; */
/* 	      break; */
/* 	    case 2: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax[1]; */
/* 	      break; */
/* 	    case 3: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax[2]; */
/* 	      break; */
/* 	    case 4: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=(double)products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Rmax; */
/* 	      break; */
/* #ifdef TWO_LPT */
/* 	    case 5: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_2LPT[0]; */
/* 	      break; */
/* 	    case 6: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_2LPT[1]; */
/* 	      break; */
/* 	    case 7: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_2LPT[2]; */
/* 	      break; */
/* #ifdef THREE_LPT */
/* 	    case 8: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_1[0]; */
/* 	      break; */
/* 	    case 9: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_1[1]; */
/* 	      break; */
/* 	    case 10: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_1[2]; */
/* 	      break; */
/* 	    case 11: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_2[0]; */
/* 	      break; */
/* 	    case 12: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_2[1]; */
/* 	      break; */
/* 	    case 13: */
/* 	      for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 		for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 		  plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]=products[ix + (MyGrids[MyGrid].GSlocal[_x_]) * (iy + local_z * MyGrids[MyGrid].GSlocal[_y_])].Vmax_3LPT_2[2]; */
/* 	      break; */
/* #endif */
/* #endif */
/* 	    default: */
/* 	      return 1; */
/* 	      break; */
/* 	    } */
/* 	} */

/*       /\* if the plane does not belong to task 0 then it must be communicated *\/ */
/*       if (belongs[global_z]) */
/*       	{ */
/*       	  if (ThisTask==belongs[global_z]) */
/*       	    { */
/*       	      /\* Task belongs[iz] sends the plane to task 0 *\/ */
/*       	      MPI_Send(plane, planesize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); */
	      
/*       	    } */
/*       	  else if (!ThisTask) */
/*       	    { */
/*       	      /\* task 0 receives the info from task belongs[global_z] and writes it on the file *\/ */
/*       	      MPI_Recv(plane, planesize, MPI_DOUBLE, belongs[global_z], 0, MPI_COMM_WORLD, &status); */
/*       	    } */
/*       	} */

/*       if (!ThisTask) */
/* 	{ */
/* #ifdef ASCII_WRITE */
/* 	  /\* writes the plane *\/ */
/* 	  for (iy=0; iy<MyGrids[MyGrid].GSglobal[_y_]; iy++) */
/* 	    for (ix=0; ix<MyGrids[MyGrid].GSglobal[_x_]; ix++) */
/* 	      { */
/* 		fprintf(file," %3d %3d %3d %13.6g\n", */
/* 			ix,iy,global_z, */
/* 			plane[ix + iy*MyGrids[MyGrid].GSglobal[_x_]]); */
/* 		/\* last=1; *\/ */
/* 		/\* if (!((1+ix + iy*MyGrids[MyGrid].GSglobal[_x_])%8)) *\/ */
/* 		/\*   { *\/ */
/* 		/\*     fprintf(file,"\n"); *\/ */
/* 		/\*     last=0; *\/ */
/* 		/\*   } *\/ */
/* 	      } */
/* 	  /\* if (last) *\/ */
/* 	  /\*   fprintf(file,"\n"); *\/ */
/* #else */
/* #ifdef FORTRAN */
/* 	  idummy=planesize*sizeof(double); */
/* 	  fwrite(&idummy,sizeof(int),1,file); */
/* #endif */
/* 	  fwrite(plane,sizeof(double),planesize,file); */
/* #ifdef FORTRAN */
/* 	  fwrite(&idummy,sizeof(int),1,file); */
/* #endif */
/* #endif */
/* 	} */
/*     } */

  MPI_Barrier(MPI_COMM_WORLD);

  if (!ThisTask)
    {
      fclose(file);
      printf("Written file %s\n",filename);
    }


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
