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

void WriteBlockName(FILE *,int, char*);
void my_strcpy(char *,char *, int);

//#define LONGIDS
#define WRITE_INFO_BLOCK
//#define POS_IN_KPC
#define FORCE_PARTICLE_NUMBER      /* DO NOT ACTIVATE THESE */
//#define ONLY_LPT_DISPLACEMENTS     /* TWO OPTIONS TOGETHER */
#define WRITE_FMAX_TO_SNAPSHOT

#if defined(FORCE_PARTICLE_NUMBER) && defined(ONLY_LPT_DISPLACEMENTS)
#error Trying to compile with FORCE_PARTICLE_NUMBER and ONLY_LPT_DISPLACEMENTS together
#endif

#ifdef LONGIDS
#define MYIDTYPE unsigned long long int
#else
#define MYIDTYPE unsigned int
#endif

typedef struct
{
  unsigned NPart[6];
  double   Mass[6];
  double   Time;
  double   RedShift;
  int      flag_sfr;
  int      flag_feedback;
  unsigned NPartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  int      flag_stellarage;
  int      flag_metals;
  unsigned npartTotalHighWord[6];
  int      flag_entropy_instead_u;
  int      flag_metalcooling;
  int      flag_stellarevolution;
  char     fill[52];  /* filling to 256 Bytes */
} SnapshotHeader_data;

#ifdef WRITE_INFO_BLOCK
typedef struct
{
  char name[4], type[8];
  int ndim, active[6];
} InfoBlock_data;
#endif


/* In this routine each task outputs all the particles belonging to good groups.
   If FORCE_PARTICLE_NUMBER is not active, some particles will be duplicated 
   but if FORCE_PARTICLE_NUMBER is active, some groups will have 
   fewer particles that their putative mass */
int write_snapshot(int iout)
{
  /* writes displacements of all particles in a GADGET format 
     particles in halos are distributed around the center of mass
     as NFW spheres with virial velocity distribution */

  int NTasksPerFile,collector,itask,next,ThisFile,index,npart,i,j,dummy,myNpart,NInFile,
    good_particle;
  int SGrid[3],Lgridxy,ibox,jbox,kbox,kk,global_x,global_y,global_z, part;
  unsigned long long int myNtotal, longdummy, NPartTot;
  double vfact, rvir, MyPos[3], *MyVel, Dz, conc, rnd, nfwfac, area, xrnd, yrnd, 
    probfunc, u, theta, sigma;
  MyVel=MyPos;
  char filename[LBLENGTH],wn[5],wt[9];
  FILE *file;
  SnapshotHeader_data Header;
  MYIDTYPE *ID,*GRID;
  MPI_Status status;
  typedef struct
  {
    float axis[3];
  } AuxStruct;
  
  AuxStruct *Pos,*Vel;
#ifdef WRITE_FMAX_TO_SNAPSHOT
  float *FMAX;
  int *RMAX;
#endif

#ifdef WRITE_INFO_BLOCK
#ifndef ONLY_LPT_DISPLACEMENTS
#ifndef WRITE_FMAX_TO_SNAPSHOT
#define NBLOCKS 4
#else
#define NBLOCKS 6
#endif
#else
#define NBLOCKS 3
#endif
  InfoBlock_data InfoBlock[NBLOCKS];
#endif
  
#ifdef ROTATE_BOX
  static int rot[3]={1,2,0};
#else
  static int rot[3]={0,1,2};
#endif


#ifdef ONLY_LPT_DISPLACEMENTS
  if (!ThisTask)
    {
      printf("directive ONLY_LPT_DISPLACEMENTS is present\n");
      printf("the code will use LPT displacements for all particles\n");
    }
#endif

#ifdef FORCE_PARTICLE_NUMBER
  if (!ThisTask)
    {
      printf("directive FORCE_PARTICLE_NUMBER is present\n");
      printf("the code will conserve the number of particles and will not check for group consistency\n");
    }
#endif



  /* Snapshots are written only if fragmentation is done in one slice */
  if (NSlices>1)
    return 0;

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  SGrid[0]=(double)subbox.stabl_x;
  SGrid[1]=(double)subbox.stabl_y;
  SGrid[2]=(double)subbox.stabl_z;
  NPartTot = 
    (unsigned long long int)MyGrids[0].GSglobal_x * 
    (unsigned long long int)MyGrids[0].GSglobal_y * 
    (unsigned long long int)MyGrids[0].GSglobal_z;
  Lgridxy = subbox.Lgwbl_x * subbox.Lgwbl_y;

  Dz = GrowingMode(outputs.z[iout],params.k_for_GM);

  /* Each task counts the number of particles that will be output */
  /* first loop is on good particles that are not in groups */
  for (i=0, myNpart=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	myNpart++;
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      myNpart+=groups[i].Mass;
#endif


  /* Task 0 counts the total number of particles that will be contained in the snapshot */
  longdummy=(unsigned long long)myNpart;
  MPI_Reduce(&longdummy, &myNtotal, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    printf("[%s] The snapshot will contain a total of %Ld particles (true total number: %Ld, duplication: %f)\n",
	   fdate(), myNtotal, NPartTot, (double)(myNtotal-NPartTot)/(double)NPartTot);

  /* collector task computes the total number of particles in the file */
  for (next=0; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (!next)
	{
	  if (ThisTask==collector)  
	    NInFile=myNpart;
	}
      else
	{
	  if (ThisTask==collector)
	    {
	      npart=0;
	      MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	      NInFile+=npart;
	    }
	  else if (ThisTask==itask)
	    MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	}
    }

  /* The collector task opens the file and writes the header */
  if (ThisTask==collector)
    {
      if (params.NumFiles>1)
	sprintf(filename,"pinocchio.%6.4f.%s.snapshot.out.%d",
		outputs.z[iout],params.RunFlag,ThisFile);
      else
	sprintf(filename,"pinocchio.%6.4f.%s.snapshot.out",
		outputs.z[iout],params.RunFlag);

      if (!ThisTask)
	printf("[%s] Task 0 will write %d particles in snapshot file %s\n",fdate(),NInFile,filename);

      if ( (file=fopen(filename,"w"))==0x0)
	{
	  printf("Error on Task 0: cannot open file %s\n",filename);
	  return 1;
	}

      /* writes information on the snapshot header */
      memset(&Header, 0, sizeof(SnapshotHeader_data));

      Header.NPart[1]=NInFile;
      Header.Mass[1]=params.ParticleMass*params.Hubble100*1e-10;
      Header.NPartTotal[1]=(unsigned int)( (myNtotal<<32) >>32 );
      Header.npartTotalHighWord[1]=(unsigned int)(myNtotal>>32);
      Header.Time=1./(1+outputs.z[iout]);
      Header.RedShift=outputs.z[iout];
      Header.flag_sfr=0;
      Header.flag_feedback=0;
      Header.flag_cooling=0;
      Header.num_files=params.NumFiles;
      Header.BoxSize=params.BoxSize_h100;
#ifdef POS_IN_KPC
      Header.BoxSize*=1000.;
#endif
      Header.Omega0=params.Omega0;
      Header.OmegaLambda=params.OmegaLambda;
      Header.HubbleParam=params.Hubble100;
      Header.flag_stellarage=0;
      Header.flag_metals=0;
      Header.flag_entropy_instead_u=0;
      Header.flag_metalcooling=0;
      Header.flag_stellarevolution=0;

      dummy=sizeof(SnapshotHeader_data);
      WriteBlockName(file,dummy,"HEAD");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&Header, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);

#ifdef WRITE_INFO_BLOCK
      /* info block */
      memset(&InfoBlock, 0, NBLOCKS*sizeof(InfoBlock_data));

      int tb=0;
      my_strcpy(InfoBlock[tb].name,"ID  ",4);
#ifdef LONGIDS
      my_strcpy(InfoBlock[tb].type,"LLONG   ",8);
#else
      my_strcpy(InfoBlock[tb].type,"LONG    ",8);
#endif
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"POS ",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"VEL ",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

#ifndef ONLY_LPT_DISPLACEMENTS
      my_strcpy(InfoBlock[tb].name,"GRID",4);
#ifdef LONGIDS
      my_strcpy(InfoBlock[tb].type,"LLONG   ",8);
#else
      my_strcpy(InfoBlock[tb].type,"LONG    ",8);
#endif
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

#ifdef WRITE_FMAX_TO_SNAPSHOT
      my_strcpy(InfoBlock[tb].name,"FMAX",4);
      my_strcpy(InfoBlock[tb].type,"FLOAT   ",8);
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"RMAX",4);
      my_strcpy(InfoBlock[tb].type,"LONG    ",8);
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=1;
      if (!ThisTask)
 	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
     ++tb;

#endif
#endif

      dummy=NBLOCKS*sizeof(InfoBlock_data);
      WriteBlockName(file,dummy,"INFO");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&InfoBlock, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);
#endif

    }

  /* each task builds its catalogue: IDs */
  ID = (MYIDTYPE*)malloc(myNpart * sizeof(MYIDTYPE));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global_x = ibox + subbox.stabl_x;
	  if (global_x<0) global_x+=MyGrids[0].GSglobal_x;
	  if (global_x>=MyGrids[0].GSglobal_x) global_x-=MyGrids[0].GSglobal_x;
	  global_y = jbox + subbox.stabl_y;
	  if (global_y<0) global_y+=MyGrids[0].GSglobal_y;
	  if (global_y>=MyGrids[0].GSglobal_y) global_y-=MyGrids[0].GSglobal_y;
	  global_z = kbox + subbox.stabl_z;
	  if (global_z<0) global_z+=MyGrids[0].GSglobal_z;
	  if (global_z>=MyGrids[0].GSglobal_z) global_z-=MyGrids[0].GSglobal_z;
#ifdef ROTATE_BOX
	  ID[index]=1+(MYIDTYPE)global_x + 
	    ( (MYIDTYPE)global_z + 
	      (MYIDTYPE)global_y * (MYIDTYPE)MyGrids[0].GSglobal_z ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#else
	  ID[index]=1+(MYIDTYPE)global_x + 
	    ( (MYIDTYPE)global_y + 
	      (MYIDTYPE)global_z * (MYIDTYPE)MyGrids[0].GSglobal_y ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#endif
	  ++index;
	}
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	next=groups[i].point;
	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    kbox=next/Lgridxy;
	    kk=next-kbox*Lgridxy;
	    jbox=kk/subbox.Lgwbl_x;
	    ibox=kk-jbox*subbox.Lgwbl_x;
	    /* particle coordinates in the box, imposing PBCs */
	    global_x = ibox + subbox.stabl_x;
	    if (global_x<0) global_x+=MyGrids[0].GSglobal_x;
	    if (global_x>=MyGrids[0].GSglobal_x) global_x-=MyGrids[0].GSglobal_x;
	    global_y = jbox + subbox.stabl_y;
	    if (global_y<0) global_y+=MyGrids[0].GSglobal_y;
	    if (global_y>=MyGrids[0].GSglobal_y) global_y-=MyGrids[0].GSglobal_y;
	    global_z = kbox + subbox.stabl_z;
	    if (global_z<0) global_z+=MyGrids[0].GSglobal_z;
	    if (global_z>=MyGrids[0].GSglobal_z) global_z-=MyGrids[0].GSglobal_z;
#ifdef ROTATE_BOX
	    ID[index]=1+(MYIDTYPE)global_x + 
	      ( (MYIDTYPE)global_z + 
		(MYIDTYPE)global_y * (MYIDTYPE)MyGrids[0].GSglobal_z ) * 
	      (MYIDTYPE)MyGrids[0].GSglobal_x;
#else
	    ID[index]=1+(MYIDTYPE)global_x + 
	      ( (MYIDTYPE)global_y + 
		(MYIDTYPE)global_z * (MYIDTYPE)MyGrids[0].GSglobal_y ) * 
	      (MYIDTYPE)MyGrids[0].GSglobal_x;
#endif
	    ++index;
	    next=linking_list[next];
	  }
      }

#endif

  /* writing of IDs */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(MYIDTYPE);
      WriteBlockName(file,dummy,"ID  ");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(ID, sizeof(MYIDTYPE), myNpart, file); 
      free(ID);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  ID=(MYIDTYPE*)malloc(npart * sizeof(MYIDTYPE));
	  MPI_Recv(ID, npart*sizeof(MYIDTYPE), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(ID, sizeof(MYIDTYPE), npart, file);
	  free(ID);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(ID, myNpart*sizeof(MYIDTYPE), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(ID);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


  /* each task builds its catalogue: Pos */
  Pos = (AuxStruct *)calloc(myNpart , sizeof(AuxStruct));

  /* the first loop is on particles outside groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global_x = ibox + subbox.stabl_x;
	  if (global_x<0) global_x+=MyGrids[0].GSglobal_x;
	  if (global_x>=MyGrids[0].GSglobal_x) global_x-=MyGrids[0].GSglobal_x;
	  global_y = jbox + subbox.stabl_y;
	  if (global_y<0) global_y+=MyGrids[0].GSglobal_y;
	  if (global_y>=MyGrids[0].GSglobal_y) global_y-=MyGrids[0].GSglobal_y;
	  global_z = kbox + subbox.stabl_z;
	  if (global_z<0) global_z+=MyGrids[0].GSglobal_z;
	  if (global_z>=MyGrids[0].GSglobal_z) global_z-=MyGrids[0].GSglobal_z;

	  /* particles outside groups (uncollapsed and filament particles)
	     are displaced with LPT */
#ifdef FORCE_PARTICLE_NUMBER
	  if (group_ID[i]<=FILAMENT)
	    {
#endif
	      set_point(global_x,global_y,global_z,i,outputs.F[iout],&obj);
	      for (j=0; j<3; j++)
		{
		  Pos[index].axis[j]=(float)(q2x(rot[j], &obj, ORDER_FOR_CATALOG) * 
					     params.InterPartDist * params.Hubble100);
		  if (Pos[index].axis[j] >= params.BoxSize_h100)
		    Pos[index].axis[j] -= (float)params.BoxSize_h100;
		  if (Pos[index].axis[j] < 0.0) 
		    Pos[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
		  Pos[index].axis[j] *= 1000.;
#endif
		}
#ifdef FORCE_PARTICLE_NUMBER
	    }
	  else
	    {
	      set_obj(group_ID[i], outputs.F[iout], &obj);
	      for (j=0; j<3; j++)
		MyPos[j]=(float)((q2x(rot[j], &obj, ORDER_FOR_CATALOG) + SGrid[rot[j]]) *
				 params.InterPartDist * params.Hubble100);
	      /* particles in the group are distributed as NFW */
	      /* Concentration taken from Bhattacharya, et al. 2013 */
	      conc = pow(Dz,0.54) * 5.9
		* pow( (1.12*pow(groups[group_ID[i]].Mass*params.ParticleMass*params.Hubble100 
				 /5.e13,0.3) + 0.53)/Dz , -0.35);
	      rvir = pow(0.01 * GRAVITY * groups[group_ID[i]].Mass*params.ParticleMass /
		   pow(Hubble(outputs.z[iout]),2.0), 1./3.) * (1. + outputs.z[iout]) * params.Hubble100;

	      part=0;
	      do
		{
		  rnd = gsl_rng_uniform(random_generator);
		  nfwfac = log(1.+conc)-conc/(1.+conc);
		  area = 1.1*conc/(4.*nfwfac);
		  xrnd = rnd*area/(1.1*conc/(4.*nfwfac));
		  rnd = gsl_rng_uniform(random_generator);
		  yrnd = rnd*1.1*conc*xrnd/(4.*nfwfac);
		  probfunc = conc*conc*xrnd/pow(1.+conc*xrnd,2)/nfwfac;
		  if (yrnd <= probfunc)
		    {
		      /* From: http://mathworld.wolfram.com/SpherePointPicking.html  */
		      u = -1.+2.*gsl_rng_uniform(random_generator);
		      theta = 2.*PI*gsl_rng_uniform(random_generator);

		      Pos[index].axis[0] = MyPos[0] + (float)(xrnd * rvir * sqrt(1.-u*u)*cos(theta));
		      Pos[index].axis[1] = MyPos[1] + (float)(xrnd * rvir * sqrt(1.-u*u)*sin(theta));
		      Pos[index].axis[2] = MyPos[2] + (float)(xrnd * rvir * u);

		      for (j=0; j<3; j++)
			{
			  if (Pos[index].axis[j] >= params.BoxSize_h100) 
			    Pos[index].axis[j] -= (float)params.BoxSize_h100;
			  if (Pos[index].axis[j] < 0.0) 
			    Pos[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
			  Pos[index].axis[j] *= 1000.;
#endif
			}
		      part=1;
		    }
		}
	      while (part==0);

	    }
#endif
	  ++index;

	}
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	/* the center of mass of groups is displaced with LPT */
	set_obj(i, outputs.F[iout], &obj);
	for (j=0; j<3; j++)
	  MyPos[j]=(float)((q2x(rot[j], &obj, ORDER_FOR_CATALOG) + SGrid[rot[j]]) *
			   params.InterPartDist * params.Hubble100);

	/* virial radius and concentration of the group */
	/* Concentration taken from Bhattacharya, et al. 2013 */
	conc = pow(Dz,0.54) * 5.9
	  * pow( (1.12*pow(groups[i].Mass*params.ParticleMass*params.Hubble100 /5.e13,0.3)+ 0.53)/Dz ,
		 -0.35);
	rvir = pow(0.01 * GRAVITY * groups[i].Mass*params.ParticleMass /
		   pow(Hubble(outputs.z[iout]),2.0), 1./3.) * (1. + outputs.z[iout]) * params.Hubble100;

	/* particles in the group are distributed as NFW */
	part=0;
	do
	  {
	    rnd = gsl_rng_uniform(random_generator);
	    nfwfac = log(1.+conc)-conc/(1.+conc);
	    area = 1.1*conc/(4.*nfwfac);
	    xrnd = rnd*area/(1.1*conc/(4.*nfwfac));
	    rnd = gsl_rng_uniform(random_generator);
	    yrnd = rnd*1.1*conc*xrnd/(4.*nfwfac);
	    probfunc = conc*conc*xrnd/pow(1.+conc*xrnd,2)/nfwfac;
	    if (yrnd <= probfunc)
	      {
		/* From: http://mathworld.wolfram.com/SpherePointPicking.html  */
		u = -1.+2.*gsl_rng_uniform(random_generator);
		theta = 2.*PI*gsl_rng_uniform(random_generator);

		Pos[index].axis[0] = MyPos[0] + (float)(xrnd * rvir * sqrt(1.-u*u)*cos(theta));
		Pos[index].axis[1] = MyPos[1] + (float)(xrnd * rvir * sqrt(1.-u*u)*sin(theta));
		Pos[index].axis[2] = MyPos[2] + (float)(xrnd * rvir * u);

		for (j=0; j<3; j++)
		  {
		    if (Pos[index].axis[j] >= params.BoxSize_h100) 
		      Pos[index].axis[j] -= (float)params.BoxSize_h100;
		    if (Pos[index].axis[j] < 0.0) 
		      Pos[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
		    Pos[index].axis[j] *= 1000.;
#endif
		  }
		index++;
		part++;
	      }
	  }
	while (part<groups[i].Mass);
      }
#endif

  /* writing of Pos */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"POS ");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Pos, sizeof(AuxStruct), myNpart, file);
      free(Pos);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Pos = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Pos, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Pos, sizeof(AuxStruct), npart, file);
	  free(Pos);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Pos, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Pos);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


  /* each task builds its catalogue: Vel */
  Vel = (AuxStruct *)malloc(myNpart * sizeof(AuxStruct));
  vfact=sqrt(1.+outputs.z[iout]);
  /* the first loop is on particles outside groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global_x = ibox + subbox.stabl_x;
	  if (global_x<0) global_x+=MyGrids[0].GSglobal_x;
	  if (global_x>=MyGrids[0].GSglobal_x) global_x-=MyGrids[0].GSglobal_x;
	  global_y = jbox + subbox.stabl_y;
	  if (global_y<0) global_y+=MyGrids[0].GSglobal_y;
	  if (global_y>=MyGrids[0].GSglobal_y) global_y-=MyGrids[0].GSglobal_y;
	  global_z = kbox + subbox.stabl_z;
	  if (global_z<0) global_z+=MyGrids[0].GSglobal_z;
	  if (global_z>=MyGrids[0].GSglobal_z) global_z-=MyGrids[0].GSglobal_z;
	  
	  /* particles outside groups (uncollapsed and filament particles)
	     are displaced with LPT */
#ifdef FORCE_PARTICLE_NUMBER
	  if (group_ID[i]<=FILAMENT)
	    {
#endif
	      set_point(global_x,global_y,global_z,i,outputs.F[iout],&obj);
	      for (j=0; j<3; j++)
		/* GADGET format requires velocities to be divided by sqrt(a) */
		Vel[index].axis[j]=(float)(vel(rot[j], &obj)*vfact);
#ifdef FORCE_PARTICLE_NUMBER
	    }
	  else
	    {
	      /* the center of mass of groups is displaced with LPT */
	      set_obj(group_ID[i], outputs.F[iout], &obj);
	      set_obj_vel(group_ID[i], outputs.F[iout], &obj);
	      for (j=0; j<3; j++)
		MyVel[j]=(float)(vel(rot[j], &obj)*vfact);

	      sigma = sqrt(GRAVITY * groups[group_ID[i]].Mass * params.ParticleMass / 3./ rvir);

	      for (j=0; j<3; j++)
		Vel[index].axis[j]=MyVel[j] + gsl_ran_gaussian(random_generator, sigma) * vfact;
	    }
#endif
	  ++index;

	}
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	/* the center of mass of groups is displaced with LPT */
	set_obj(i, outputs.F[iout], &obj);
	set_obj_vel(i, outputs.F[iout], &obj);
	for (j=0; j<3; j++)
	  MyVel[j]=(float)(vel(rot[j], &obj)*vfact);

	/* /\* virial radius and concentration of the group *\/ */
	/* Dz = GrowingMode(outputs.z[iout],params.k_for_GM); */
	/* /\* this is the virial radius in physical true kpc *\/ */
	/* rvir = pow(0.01 * GRAVITY * groups[i].Mass*params.ParticleMass / */
	/* 	   pow(Hubble(outputs.z[iout]),2.0), 1./3.);   */
	sigma = sqrt(GRAVITY * groups[i].Mass * params.ParticleMass / 3./ rvir);

	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    for (j=0; j<3; j++)
	      Vel[index].axis[j]=MyVel[j] + gsl_ran_gaussian(random_generator, sigma) * vfact;
	    ++index;
	  }
      }
#endif

  /* writing of Vel */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"VEL ");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Vel, sizeof(AuxStruct), myNpart, file);
      free(Vel);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Vel, sizeof(AuxStruct), npart, file);
	  free(Vel);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Vel, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Vel);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

#ifndef ONLY_LPT_DISPLACEMENTS
  /* each task builds its catalogue: GRID */
  GRID = (MYIDTYPE*)malloc(myNpart * sizeof(MYIDTYPE));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
#ifndef FORCE_PARTICLE_NUMBER
      if (group_ID[i]<=FILAMENT && good_particle)
	GRID[index++]=group_ID[i];
#else
      if (good_particle)
	{
	  if (group_ID[i]<=FILAMENT)
	    GRID[index++]=group_ID[i];
	  else
	    GRID[index++]=groups[group_ID[i]].name;
	}
#endif

    }

#ifndef FORCE_PARTICLE_NUMBER
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	next=groups[i].point;
	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    GRID[index++]=groups[i].name;
	    next=linking_list[next];
	  }
      }
#endif

  /* writing of GRIDs */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(MYIDTYPE);
      WriteBlockName(file,dummy,"GRID");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(GRID, sizeof(MYIDTYPE), myNpart, file); 
      free(GRID);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  GRID=(MYIDTYPE*)malloc(npart * sizeof(MYIDTYPE));
	  MPI_Recv(GRID, npart*sizeof(MYIDTYPE), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(GRID, sizeof(MYIDTYPE), npart, file);
	  free(GRID);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(GRID, myNpart*sizeof(MYIDTYPE), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(GRID);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

#ifdef WRITE_FMAX_TO_SNAPSHOT

  /* each task builds its catalogue: FMAX */
  FMAX = (float*)malloc(myNpart * sizeof(float));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	FMAX[index++]=frag[i].Fmax;
    }

#ifndef FORCE_PARTICLE_NUMBER
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	next=groups[i].point;
	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    FMAX[index++]=frag[next].Fmax;
	    next=linking_list[next];
	  }
      }
#endif

  /* writing of FMAX */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(float);
      WriteBlockName(file,dummy,"FMAX");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(FMAX, sizeof(float), myNpart, file); 
      free(FMAX);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  FMAX=(float*)malloc(npart * sizeof(float));
	  MPI_Recv(FMAX, npart*sizeof(float), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(FMAX, sizeof(float), npart, file);
	  free(FMAX);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(FMAX, myNpart*sizeof(float), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(FMAX);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);



  /* each task builds its catalogue: RMAX */
  RMAX = (int*)malloc(myNpart * sizeof(int));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	RMAX[index++]=frag[i].Rmax;
    }

#ifndef FORCE_PARTICLE_NUMBER
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	next=groups[i].point;
	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    RMAX[index++]=frag[next].Rmax;
	    next=linking_list[next];
	  }
      }
#endif

  /* writing of RMAX */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(int);
      WriteBlockName(file,dummy,"RMAX");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(RMAX, sizeof(int), myNpart, file); 
      free(RMAX);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  RMAX=(int*)malloc(npart * sizeof(int));
	  MPI_Recv(RMAX, npart*sizeof(int), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(RMAX, sizeof(int), npart, file);
	  free(RMAX);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(RMAX, myNpart*sizeof(int), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(RMAX);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


#endif
#endif


  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  return 0;
}



#ifdef TIMELESS_SNAPSHOT
/* In this routine a snapshot is written that contains all information
   needed to reconstruct at the post-processing level the full matter
   density field */
int write_timeless_snapshot(void)
{
/* NB this has fewer options the other code, because it is thought to 
   be used for post-processing only */

  int NTasksPerFile,collector,itask,next,ThisFile,index,npart,i,j,dummy,myNpart,NInFile,
    good_particle;
  int Lgridxy,ibox,jbox,kbox,kk,global_x,global_y,global_z;
  unsigned long long int NPartTot;
  char filename[LBLENGTH],wn[5],wt[9];
  FILE *file;
  SnapshotHeader_data Header;
  MYIDTYPE *ID;
  MPI_Status status;
  typedef struct
  {
    float axis[3];
  } AuxStruct;
  
  AuxStruct *Vel;
  float *FMAX,*ZACC;
  int *RMAX;

/* 1 ID
   2 VZEL
   3 V2
   4 V3_1
   5 V3_2
   6 FMAX
   7 RMAX
   8 ZACC */

#define NBLOCKS_T 8
  InfoBlock_data InfoBlock[NBLOCKS_T];
  
#ifdef ROTATE_BOX
  static int rot[3]={1,2,0};
#else
  static int rot[3]={0,1,2};
#endif

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  NPartTot = 
    (unsigned long long int)MyGrids[0].GSglobal_x * 
    (unsigned long long int)MyGrids[0].GSglobal_y * 
    (unsigned long long int)MyGrids[0].GSglobal_z;
  Lgridxy = subbox.Lgwbl_x * subbox.Lgwbl_y;


  /* this is the total number of particles for the task */
  for (i=0, myNpart=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	myNpart++;
    }

  /* collector task computes the total number of particles in the file */
  for (next=0; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (!next)
	{
	  if (ThisTask==collector)  
	    NInFile=myNpart;
	}
      else
	{
	  if (ThisTask==collector)
	    {
	      npart=0;
	      MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	      NInFile+=npart;
	    }
	  else if (ThisTask==itask)
	    MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	}
    }


  /* The collector task opens the file and writes the header */
  if (ThisTask==collector)
    {
      if (NSlices==1)
	{
	  if (params.NumFiles>1)
	    sprintf(filename,"pinocchio.%s.t_snapshot.out.%d",params.RunFlag,ThisFile);
	  else
	    sprintf(filename,"pinocchio.%s.t_snapshot.out",params.RunFlag);
	}
      else
	{
	  sprintf(filename,"pinocchio.%s.t_snapshot.out.%d",params.RunFlag,ThisFile+ThisSlice*params.NumFiles);
	}
	    

      if (!ThisTask)
	printf("[%s] Task 0 will write %d particles in snapshot file %s\n",fdate(),NInFile,filename);

      if ( (file=fopen(filename,"w"))==0x0)
	{
	  printf("Error on Task 0: cannot open file %s\n",filename);
	  return 1;
	}

      /* writes information on the snapshot header */
      memset(&Header, 0, sizeof(SnapshotHeader_data));

      Header.NPart[1]=NInFile;
      Header.Mass[1]=params.ParticleMass*params.Hubble100*1e-10;
      Header.NPartTotal[1]=(unsigned int)( (NPartTot<<32) >>32 );
      Header.npartTotalHighWord[1]=(unsigned int)(NPartTot>>32);
      Header.Time=1.;
      Header.RedShift=0.0;
      Header.flag_sfr=0;
      Header.flag_feedback=0;
      Header.flag_cooling=0;
      Header.num_files=params.NumFiles;
      Header.BoxSize=params.BoxSize_h100;
#ifdef POS_IN_KPC
      Header.BoxSize*=1000.;
#endif
      Header.Omega0=params.Omega0;
      Header.OmegaLambda=params.OmegaLambda;
      Header.HubbleParam=params.Hubble100;
      Header.flag_stellarage=0;
      Header.flag_metals=0;
      Header.flag_entropy_instead_u=0;
      Header.flag_metalcooling=0;
      Header.flag_stellarevolution=0;

      dummy=sizeof(SnapshotHeader_data);
      WriteBlockName(file,dummy,"HEAD");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&Header, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);

      /* info block */
      memset(&InfoBlock, 0, NBLOCKS_T*sizeof(InfoBlock_data));
      int tb=0;
      my_strcpy(InfoBlock[tb].name,"ID  ",4);
#ifdef LONGIDS
      my_strcpy(InfoBlock[tb].type,"LLONG   ",8);
#else
      my_strcpy(InfoBlock[tb].type,"LONG    ",8);
#endif
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)	
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"VZEL",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

#ifdef TWO_LPT
      my_strcpy(InfoBlock[tb].name,"V2  ",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

#ifdef THREE_LPT
      my_strcpy(InfoBlock[tb].name,"V3_1",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"V3_2",4);
      my_strcpy(InfoBlock[tb].type,"FLOATN  ",8);
      InfoBlock[tb].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;
#endif
#endif

      my_strcpy(InfoBlock[tb].name,"FMAX",4);
      my_strcpy(InfoBlock[tb].type,"FLOAT   ",8);
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"RMAX",4);
      my_strcpy(InfoBlock[tb].type,"LONG    ",8);
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      my_strcpy(InfoBlock[tb].name,"ZACC",4);
      my_strcpy(InfoBlock[tb].type,"FLOAT   ",8);
      InfoBlock[tb].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[tb].active[i]=0;
      InfoBlock[tb].active[1]=1;
      if (!ThisTask)
	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[9]='\0';
	  printf("name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 /* InfoBlock[tb].name, */
		 /* InfoBlock[tb].type, */
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
      ++tb;

      dummy=NBLOCKS_T*sizeof(InfoBlock_data);
      WriteBlockName(file,dummy,"INFO");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&InfoBlock, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);

    }


  /* ****************** ID ****************** */

  /* each task builds its catalogue: IDs */
  ID = (MYIDTYPE*)malloc(myNpart * sizeof(MYIDTYPE));
  /* loop on good particles */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global_x = ibox + subbox.stabl_x;
	  if (global_x<0) global_x+=MyGrids[0].GSglobal_x;
	  if (global_x>=MyGrids[0].GSglobal_x) global_x-=MyGrids[0].GSglobal_x;
	  global_y = jbox + subbox.stabl_y;
	  if (global_y<0) global_y+=MyGrids[0].GSglobal_y;
	  if (global_y>=MyGrids[0].GSglobal_y) global_y-=MyGrids[0].GSglobal_y;
	  global_z = kbox + subbox.stabl_z;
	  if (global_z<0) global_z+=MyGrids[0].GSglobal_z;
	  if (global_z>=MyGrids[0].GSglobal_z) global_z-=MyGrids[0].GSglobal_z;
#ifdef ROTATE_BOX
	  ID[index]=1+(MYIDTYPE)global_x + 
	    ( (MYIDTYPE)global_z + 
	      (MYIDTYPE)global_y * (MYIDTYPE)MyGrids[0].GSglobal_z ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#else
	  ID[index]=1+(MYIDTYPE)global_x + 
	    ( (MYIDTYPE)global_y + 
	      (MYIDTYPE)global_z * (MYIDTYPE)MyGrids[0].GSglobal_y ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#endif
	  ++index;
	}
    }

  /* writing of IDs */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(MYIDTYPE);
      WriteBlockName(file,dummy,"ID  ");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(ID, sizeof(MYIDTYPE), myNpart, file); 
      free(ID);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  ID=(MYIDTYPE*)malloc(npart * sizeof(MYIDTYPE));
	  MPI_Recv(ID, npart*sizeof(MYIDTYPE), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(ID, sizeof(MYIDTYPE), npart, file);
	  free(ID);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(ID, myNpart*sizeof(MYIDTYPE), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(ID);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

  /* ****************** VZEL ****************** */

  /* each task builds its catalogue: Vzel */
  Vel = (AuxStruct *)malloc(myNpart * sizeof(AuxStruct));
  /* loop on good particles */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	{
	  for (j=0; j<3; j++)
	    Vel[index].axis[j]=(float)(frag[i].Vel[rot[j]]);
	  ++index;
	}
    }

  /* writing of Vel */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"VZEL");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Vel, sizeof(AuxStruct), myNpart, file);
      free(Vel);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Vel, sizeof(AuxStruct), npart, file);
	  free(Vel);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Vel, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Vel);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

#ifdef TWO_LPT
  /* ****************** V2 ****************** */

  /* each task builds its catalogue: Vzel */
  Vel = (AuxStruct *)malloc(myNpart * sizeof(AuxStruct));
  /* loop on good particles */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	{
	  for (j=0; j<3; j++)
	    Vel[index].axis[j]=(float)(frag[i].Vel_2LPT[rot[j]]);
	  ++index;
	}
    }

  /* writing of V2 */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"V2  ");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Vel, sizeof(AuxStruct), myNpart, file);
      free(Vel);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Vel, sizeof(AuxStruct), npart, file);
	  free(Vel);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Vel, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Vel);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

#ifdef THREE_LPT
  /* ****************** V3_1 ****************** */

  /* each task builds its catalogue: Vzel */
  Vel = (AuxStruct *)malloc(myNpart * sizeof(AuxStruct));
  /* loop on good particles */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	{
	  for (j=0; j<3; j++)
	    Vel[index].axis[j]=(float)(frag[i].Vel_3LPT_1[rot[j]]);
	  ++index;
	}
    }

  /* writing of V3_1 */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"V3_1");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Vel, sizeof(AuxStruct), myNpart, file);
      free(Vel);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Vel, sizeof(AuxStruct), npart, file);
	  free(Vel);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Vel, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Vel);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


  /* ****************** V3_2 ****************** */

  /* each task builds its catalogue: Vzel */
  Vel = (AuxStruct *)malloc(myNpart * sizeof(AuxStruct));
  /* loop on good particles */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	{
	  for (j=0; j<3; j++)
	    Vel[index].axis[j]=(float)(frag[i].Vel_3LPT_2[rot[j]]);
	  ++index;
	}
    }

  /* writing of V3_1 */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(AuxStruct);
      WriteBlockName(file,dummy,"V3_2");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(Vel, sizeof(AuxStruct), myNpart, file);
      free(Vel);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	  MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(Vel, sizeof(AuxStruct), npart, file);
	  free(Vel);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(Vel, sizeof(AuxStruct)*myNpart, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(Vel);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

#endif
#endif
  /* ****************** FMAX ****************** */

  /* each task builds its catalogue: FMAX */
  FMAX = (float*)malloc(myNpart * sizeof(float));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	FMAX[index++]=frag[i].Fmax;
    }

  /* writing of FMAX */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(float);
      WriteBlockName(file,dummy,"FMAX");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(FMAX, sizeof(float), myNpart, file); 
      free(FMAX);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  FMAX=(float*)malloc(npart * sizeof(float));
	  MPI_Recv(FMAX, npart*sizeof(float), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(FMAX, sizeof(float), npart, file);
	  free(FMAX);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(FMAX, myNpart*sizeof(float), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(FMAX);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


  /* ****************** RMAX ****************** */

  /* each task builds its catalogue: RMAX */
  RMAX = (int*)malloc(myNpart * sizeof(int));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	RMAX[index++]=frag[i].Rmax;
    }

  /* writing of RMAX */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(int);
      WriteBlockName(file,dummy,"RMAX");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(RMAX, sizeof(int), myNpart, file); 
      free(RMAX);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  RMAX=(int*)malloc(npart * sizeof(int));
	  MPI_Recv(RMAX, npart*sizeof(int), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(RMAX, sizeof(int), npart, file);
	  free(RMAX);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(RMAX, myNpart*sizeof(int), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(RMAX);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

  /* ****************** ZACC ****************** */

  /* each task builds its catalogue: ZACC */
  ZACC = (float*)malloc(myNpart * sizeof(float));
  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl_x;
      ibox=kk-jbox*subbox.Lgwbl_x;
      good_particle = ( ibox>=subbox.safe_x && ibox<subbox.Lgwbl_x-subbox.safe_x && 
			jbox>=subbox.safe_y && jbox<subbox.Lgwbl_y-subbox.safe_y && 
			kbox>=subbox.safe_z && kbox<subbox.Lgwbl_z-subbox.safe_z );
      if (good_particle)
	ZACC[index++]=frag[i].zacc;
    }

  /* writing of ZACC */
  if (ThisTask==collector)
    {
      dummy=NInFile*sizeof(float);
      WriteBlockName(file,dummy,"ZACC");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(ZACC, sizeof(float), myNpart, file); 
      free(ZACC);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  ZACC=(float*)malloc(npart * sizeof(float));
	  MPI_Recv(ZACC, npart*sizeof(float), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(ZACC, sizeof(float), npart, file);
	  free(ZACC);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(ZACC, myNpart*sizeof(float), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	  free(ZACC);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);


  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  return 0;
}
#endif


int write_LPT_snapshot(double redshift)
{
  // TODO: growth rate obtained directly in k-space

  /* writes displacements of all particles in a GADGET format */

  int NTasksPerFile,collector,itask,next,ThisFile,index,x,y,z,NInFile,npart,i,dummy;
  int SGrid[3],Q[3];  /* LGrid[3],GGrid[3]; */
  unsigned long long int NPartTot;
  double vf, G;
  char filename[LBLENGTH];
  FILE *file;
  SnapshotHeader_data Header;
  MYIDTYPE *ID;
  MPI_Status status;
  typedef struct 
  {
    float axis[3];
  } AuxStruct;
  
  AuxStruct *Pos,*Vel;

#ifdef WRITE_INFO_BLOCK
  InfoBlock_data InfoBlock[3];
#endif

#ifdef TWO_LPT
  double vf2,G2;
#ifdef THREE_LPT
  //double G31,G32;
#endif
#endif
#ifdef ROTATE_BOX
  static int rot[3]={1,2,0};
#else
  static int rot[3]={0,1,2};
#endif

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  /* GADGET format requires velocities to be divided by sqrt(a) */
  G=GrowingMode(redshift,params.k_for_GM);
  vf = params.InterPartDist*fomega(redshift,params.k_for_GM)*Hubble(redshift)/sqrt(1.+redshift) * G;
#ifdef TWO_LPT
  G2  = GrowingMode_2LPT(redshift,params.k_for_GM);
  vf2 = params.InterPartDist*fomega_2LPT(redshift,params.k_for_GM)*Hubble(redshift)/sqrt(1.+redshift) * G2;
#ifdef THREE_LPT
  //G31  = GrowingMode_3LPT_1(redshift,params.k_for_GM);
  //G32  = GrowingMode_3LPT_2(redshift,params.k_for_GM);
#endif
#endif
  SGrid[0]=(double)MyGrids[0].GSstart_x;
  SGrid[1]=(double)MyGrids[0].GSstart_y;
  SGrid[2]=(double)MyGrids[0].GSstart_z;

  NPartTot = (unsigned long long int)MyGrids[0].GSglobal_x * 
    (unsigned long long int)MyGrids[0].GSglobal_y * 
    (unsigned long long int)MyGrids[0].GSglobal_z;

  /* collector task computes the total number of particles in the file */
  for (next=0; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (!next)
	{
	  if (ThisTask==collector)  
	    NInFile=MyGrids[0].total_local_size;
	}
      else
	{
	  if (ThisTask==collector)
	    {
	      npart=0;
	      MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	      NInFile+=npart;
	    }
	  else if (ThisTask==itask)
	    MPI_Send(&MyGrids[0].total_local_size, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	}
    }
	      

  /* The collector task opens the file and writes the header */
  if (ThisTask==collector)
    {
      if (params.NumFiles>1)
	sprintf(filename,"pinocchioICs.%6.4f.%s.%d",
		redshift,params.RunFlag,ThisFile);
      else
	sprintf(filename,"pinocchioICs.%6.4f.%s",
		redshift,params.RunFlag);

      if (!ThisTask)
	printf("[%s] Writing file %s\n",fdate(),filename);

      if ( (file=fopen(filename,"w"))==0x0)
	{
	  printf("Error on Task 0: cannot open file %s\n",filename);
	  return 1;
	}

      /* writes information on the snapshot header */
      memset(&Header, 0, sizeof(SnapshotHeader_data));

      Header.NPart[1]=NInFile;
      Header.Mass[1]=params.ParticleMass*params.Hubble100*1e-10;
      Header.NPartTotal[1]=(unsigned int)( (NPartTot<<32) >>32 );
      Header.npartTotalHighWord[1]=(unsigned int)(NPartTot>>32);
      Header.Time=1./(1+redshift);
      Header.RedShift=redshift;
      Header.flag_sfr=0;
      Header.flag_feedback=0;
      Header.flag_cooling=0;
      Header.num_files=params.NumFiles;
      Header.BoxSize=params.BoxSize_h100;
      Header.Omega0=params.Omega0;
      Header.OmegaLambda=params.OmegaLambda;
      Header.HubbleParam=params.Hubble100;
      Header.flag_stellarage=0;
      Header.flag_metals=0;
      Header.flag_entropy_instead_u=0;
      Header.flag_metalcooling=0;
      Header.flag_stellarevolution=0;

      dummy=sizeof(SnapshotHeader_data);
      WriteBlockName(file,dummy,"HEAD");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&Header, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);
#ifdef WRITE_INFO_BLOCK
      /* info block */
      memset(&InfoBlock, 0, 3*sizeof(InfoBlock_data));

      my_strcpy(InfoBlock[0].name,"POS ",4);
      my_strcpy(InfoBlock[0].type,"FLOATN  ",8);
      InfoBlock[0].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[0].active[i]=1;

      my_strcpy(InfoBlock[1].name,"VEL ",4);
      my_strcpy(InfoBlock[1].type,"FLOATN  ",8);
      InfoBlock[1].ndim=3;
      for (i=0; i<6; i++)
	InfoBlock[1].active[i]=1;

      my_strcpy(InfoBlock[2].name,"ID  ",4);
#ifdef LONGIDS
      my_strcpy(InfoBlock[2].type,"LLONG   ",8);
#else
      my_strcpy(InfoBlock[2].type,"LONG    ",8);
#endif
      InfoBlock[2].ndim=1;
      for (i=0; i<6; i++)
	InfoBlock[2].active[i]=1;

      dummy=3*sizeof(InfoBlock_data);
      WriteBlockName(file,dummy,"INFO");
      fwrite(&dummy, sizeof(dummy), 1, file);
      fwrite(&InfoBlock, dummy, 1, file);
      fwrite(&dummy, sizeof(dummy), 1, file);
#endif

    }


  /* each task builds its catalogue: Pos */
  if (MyGrids[0].total_local_size)
    {
      Pos = (AuxStruct *)malloc(MyGrids[0].total_local_size * sizeof(AuxStruct));

      for (z=0; z<MyGrids[0].GSlocal_z; z++)
	{
	  Q[2]=z;
	  for (y=0; y<MyGrids[0].GSlocal_y; y++)
	    {
	      Q[1]=y;
	      for (x=0; x<MyGrids[0].GSlocal_x; x++)
		{
		  Q[0]=x;

		  index = x + MyGrids[0].GSlocal_x *(y + z* MyGrids[0].GSlocal_y);

		  for (i=0; i<3; i++)
		    {

		      Pos[index].axis[i] = (float)(( (double)(Q[rot[i]]+SGrid[rot[i]]) + G*VEL_for_displ[rot[i]][index] )
						   * params.InterPartDist * params.Hubble100)
#ifdef TWO_LPT
			+ (float)(G2 * VEL2_for_displ[rot[i]][index] * params.InterPartDist * params.Hubble100)
#ifdef THREE_LPT

#endif
#endif
			;

		      if (Pos[index].axis[i] < 0) Pos[index].axis[i] += (float)params.BoxSize_h100;
		      if (Pos[index].axis[i] >= params.BoxSize_h100) Pos[index].axis[i] -= (float)params.BoxSize_h100;
		    }
		}
	    }
	}
    }

  /* writing of Pos */
  if (ThisTask==collector)
    {
      if (MyGrids[0].total_local_size)
	{
	  dummy=NInFile*sizeof(AuxStruct);
	  WriteBlockName(file,dummy,"POS ");
	  fwrite(&dummy, sizeof(dummy), 1, file);
	  fwrite(Pos, sizeof(AuxStruct), MyGrids[0].total_local_size, file);
	  free(Pos);
	}
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  if (npart)
	    {
	      Pos = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	      MPI_Recv(Pos, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	      fwrite(Pos, sizeof(AuxStruct), npart, file);
	      free(Pos);
	    }
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&MyGrids[0].total_local_size, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  if (MyGrids[0].total_local_size)
	    {
	      MPI_Send(Pos, sizeof(AuxStruct)*MyGrids[0].total_local_size, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	      free(Pos);
	    }
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy, sizeof(dummy), 1, file);

  /* each task builds its catalogue: Vel */
  if (MyGrids[0].total_local_size)
    {
      Vel = (AuxStruct *)malloc(MyGrids[0].total_local_size * sizeof(AuxStruct));

      for (z=0; z<MyGrids[0].GSlocal_z; z++)
	for (y=0; y<MyGrids[0].GSlocal_y; y++)
	  for (x=0; x<MyGrids[0].GSlocal_x; x++)
	    {
	      index = x + MyGrids[0].GSlocal_x *(y + z* MyGrids[0].GSlocal_y);

	      for (i=0; i<3; i++)
		{
 		  Vel[index].axis[i] = (float)(vf*VEL_for_displ[rot[i]][index])
#ifdef TWO_LPT
		    + (float)(vf2*VEL2_for_displ[rot[i]][index])
#ifdef THREE_LPT

#endif
#endif
		  ;
		}
	    }
    }

  /* writing of Vel */
  if (ThisTask==collector)
    {
      if (MyGrids[0].total_local_size)
	{
	  dummy=NInFile*sizeof(AuxStruct);
	  WriteBlockName(file,dummy,"VEL ");
	  fwrite(&dummy, sizeof(dummy), 1, file);
	  fwrite(Vel, sizeof(AuxStruct), MyGrids[0].total_local_size, file);
	  free(Vel);
	}
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  if (npart)
	    {
	      Vel = (AuxStruct *)malloc(npart * sizeof(AuxStruct));
	      MPI_Recv(Vel, npart*sizeof(AuxStruct), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	      fwrite(Vel, sizeof(AuxStruct), npart, file);
	      free(Vel);
	    }
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&MyGrids[0].total_local_size, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  if (MyGrids[0].total_local_size)
	    {
	      MPI_Send(Vel, sizeof(AuxStruct)*MyGrids[0].total_local_size, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	      free(Vel);
	    }
	}
    }

  /* collector task closes the block and the file */
  if (ThisTask==collector)
    {
      fwrite(&dummy, sizeof(dummy), 1, file);
    }


  /* each task builds its catalogue: IDs */
  if (MyGrids[0].total_local_size)
    {
      ID = (MYIDTYPE*)malloc(MyGrids[0].total_local_size * sizeof(MYIDTYPE));

      for (z=0; z<MyGrids[0].GSlocal_z; z++)
	for (y=0; y<MyGrids[0].GSlocal_y; y++)
	  for (x=0; x<MyGrids[0].GSlocal_x; x++)
	    {
	      index = x + MyGrids[0].GSlocal_x *(y + z* MyGrids[0].GSlocal_y);

#ifdef ROTATE_BOX
	  ID[index]=1+(MYIDTYPE)(x+MyGrids[0].GSstart_x) + 
	    ( (MYIDTYPE)(z+MyGrids[0].GSstart_z) + 
	      (MYIDTYPE)(y+MyGrids[0].GSstart_y) * (MYIDTYPE)MyGrids[0].GSglobal_z ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#else
	  ID[index]=1+(MYIDTYPE)(x+MyGrids[0].GSstart_x) + 
	    ( (MYIDTYPE)(y+MyGrids[0].GSstart_y) + 
	      (MYIDTYPE)(z+MyGrids[0].GSstart_z) * (MYIDTYPE)MyGrids[0].GSglobal_y ) * 
	    (MYIDTYPE)MyGrids[0].GSglobal_x;
#endif
	  /*
#ifdef ROTATE_BOX
	      ID[index] =  1 + (MYIDTYPE)(x+MyGrids[0].GSstart_x) + MyGrids[0].GSglobal_x * 
		( (MYIDTYPE)(z+MyGrids[0].GSstart_z) +
		  (MYIDTYPE)(y+MyGrids[0].GSstart_y) * MyGrids[0].GSglobal_z);
#else
	      ID[index] =  1 + (MYIDTYPE)(z+MyGrids[0].GSstart_z) + MyGrids[0].GSglobal_z *
		( (MYIDTYPE)(y+MyGrids[0].GSstart_y) +
		  (MYIDTYPE)(x+MyGrids[0].GSstart_x) * MyGrids[0].GSglobal_y );
#endif
	  */
	    }
    }

  /* writing of IDs */
  if (ThisTask==collector)
    {
      if (MyGrids[0].total_local_size)
	{
	  dummy=NInFile*sizeof(MYIDTYPE);
	  WriteBlockName(file,dummy,"ID  ");
	  fwrite(&dummy, sizeof(dummy), 1, file);
	  fwrite(ID, sizeof(MYIDTYPE), MyGrids[0].total_local_size, file); 
	  free(ID);
	}
    }


  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  if (npart)
	    {
	      ID=(MYIDTYPE*)malloc(npart * sizeof(MYIDTYPE));
	      MPI_Recv(ID, npart*sizeof(MYIDTYPE), MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	      fwrite(ID, sizeof(MYIDTYPE), npart, file);
	      free(ID);
	    }
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&MyGrids[0].total_local_size, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  if (MyGrids[0].total_local_size)
	    {
	      MPI_Send(ID, MyGrids[0].total_local_size*sizeof(MYIDTYPE), MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	      free(ID);
	    }
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    {
      fwrite(&dummy, sizeof(dummy), 1, file);
      fclose(file);
    }

  return 0;
}





void WriteBlockName(FILE *fd,int dummy, char* LABEL)
{
int dummy2;

 dummy2=8;
 fwrite(&dummy2,sizeof(dummy2),1,fd);
 fwrite(LABEL,4,1,fd);
 dummy2=dummy+8;
 fwrite(&dummy2,sizeof(dummy2),1,fd);
 dummy2=8;
 fwrite(&dummy2,sizeof(dummy2),1,fd);
} 


void my_strcpy(char *to, char *from, int n)
{
  int i;
  for (i=0; i<n; i++)
    *(to+i)=*(from+i);
}
