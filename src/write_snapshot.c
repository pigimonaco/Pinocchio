/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco, Tom Theuns, Giuliano Taffoni, Marius Lepinzan, 
 Chiara Moretti, Luca Tornatore, David Goz, Tiago Castro
 Copyright (C) 2025
 
 github: https://github.com/pigimonaco/Pinocchio
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

//#define POS_IN_KPC   /* with this directive on, positions will be in kpc/h */


/* This file contains functions to handle gadget-like snapshots
   (format 2 with INFO block) with information on all particles 
   Possible snapshots are:
   - timeless snapshot -- gives ID and LPT fields for all particles, plus Fmax, Rmax (if required) 
     and the ZACC time at which the particle enters a halo
   - LPT snapshot -- gives ID, position and velocity of each particle, according
     to LPT at a given order
   - density -- it writes the ID and the linear density field
*/
  

#ifdef SNAPSHOT

static unsigned long long largest32 = (unsigned)1<<31;

#ifdef LONGIDS
#define MYIDTYPE unsigned long long int
#else
#define MYIDTYPE unsigned int
#endif

/* gadget header */
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

/* information for the INFO block */
typedef struct
{
  char name[4], type[8];
  int ndim, active[6];
  size_t sizeof_type;
  void *data;
} Block_data;

Block_data *InfoBlock;
int NBlocks, NextBlock;

/* for vector blocks */
typedef struct
{
  float axis[3];
} AuxStruct;

MPI_Status status;

void WriteBlockName(FILE *,unsigned long long, char*);
void my_strcpy(char *,char *, int);
int write_header();
int write_block(Block_data);
void free_block(Block_data);
int add_to_info(Block_data);
int write_info_block(void);
int initialize_ID(Block_data*);
int initialize_FMAX(Block_data*);
int initialize_RMAX(Block_data*);
int initialize_ZEL(Block_data*);
#ifdef TWO_LPT
int initialize_2LPT(Block_data*);
#ifdef THREE_LPT
int initialize_3LPT_1(Block_data*);
int initialize_3LPT_2(Block_data*);
#endif
#endif
int initialize_ZACC(Block_data*);
int initialize_GRUP(Block_data*);
int initialize_POS(Block_data*);
int initialize_VEL(Block_data*);
int initialize_density(int, Block_data*);
void set_point_timedep(double);

char filename[LBLENGTH];
FILE *file;
int NTasksPerFile,collector,ThisFile,NPartInFile,myNpart,myiout;
int *Npart_array;

/* this is the redshift at which an LPT snapshot is written */
int myiout;

/* this is used to shift particles to the final position */
//pos_data myobj;

double redshift;

int write_LPT_snapshot()
{
  /* write positions of all particles obtained with LPT */

  Block_data block;
  //myiout=0;
  //set_point_timedep(outputs.z[myiout]);
  redshift = outputs.z[0];

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%6.4f.%s.LPT_snapshot.out",outputs.z[myiout],params.RunFlag);

  /* allocates structure to handle the INFO block */
  NBlocks=4;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of POS */
  if (initialize_POS(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_VEL(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(Npart_array);
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;

}


int write_timeless_snapshot()
{
  /* writes the fmax products as in a snapshot format */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.t_snapshot.out",params.RunFlag);

  myiout=outputs.n-1;

  /* allocates structure to handle the INFO block */
#ifdef TWO_LPT
#ifdef THREE_LPT
  NBlocks=8;
#else
  NBlocks=6;
#endif
#else
  NBlocks=4;
#endif
#ifdef ADD_RMAX_TO_SNAPSHOT
  ++NBlocks;
#endif

  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing products in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef ADD_RMAX_TO_SNAPSHOT
  /* writing of RMAX */
  if (initialize_RMAX(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);
#endif

  /* writing of FMAX */
  if (initialize_FMAX(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_ZEL(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef TWO_LPT
  /* writing of VEL2 */
  if (initialize_2LPT(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef THREE_LPT
  /* writing of VL31 */
  if (initialize_3LPT_1(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VL32 */
  if (initialize_3LPT_2(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);
#endif
#endif

  /* writing of ZACC */
  if (initialize_ZACC(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of particle GROUP_ID */
  if (initialize_GRUP(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;
}

int write_density(int ThisGrid)
{
  /* writes the fmax products as in a snapshot format */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.density%d.out",params.RunFlag,ThisGrid);

  /* allocates structure to handle the INFO block */
  NBlocks=2;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing density in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header())
    return 1;

  /* writing of IDs */
  if (initialize_ID(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of density */
  if (initialize_density(ThisGrid, &block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of INFO block */
  if (write_info_block())
    return 1;
  free(InfoBlock);

  /* collector task closes the file */
  if (ThisTask==collector)
    fclose(file);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  return 0;

}


int write_header()
{

  /* This function counts the number of particles that the task will write in the snapshot
     and writes the snapshot header.
     The type of particles depend on the type of snapshot that is being written...
  */

  int itask, dummy;
  unsigned long long int myNtotal;
  char my_filename[LBLENGTH],tag[SBLENGTH];
  int *tmp_array;
  SnapshotHeader_data Header;

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  myNpart = MyGrids[0].total_local_size;

  /* Number of particles in each file */
  Npart_array = (int*)calloc(params.NumFiles, sizeof(int));
  tmp_array   = (int*)calloc(params.NumFiles, sizeof(int));
  tmp_array[collector]=myNpart;

  /* this computes the number of particles in each file */
  MPI_Reduce(tmp_array, Npart_array, params.NumFiles, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(Npart_array, params.NumFiles, MPI_INT, 0, MPI_COMM_WORLD);
  NPartInFile=Npart_array[collector];

  /* total number of particles that will be contained in the complete snapshot */
  for (myNtotal=0, itask=0; itask<params.NumFiles; itask++)
    myNtotal+=Npart_array[itask];

  /* sets arrays to zero (probably superfluous) */
  memset(tmp_array, 0, params.NumFiles*sizeof(int));
  memset(Npart_array, 0, params.NumFiles*sizeof(int));
  tmp_array[collector]=myNpart;

  /* this computes the largest number of particles a task can have for each collector */
  MPI_Reduce(tmp_array, Npart_array, params.NumFiles, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Bcast(Npart_array, params.NumFiles, MPI_INT, 0, MPI_COMM_WORLD);
  free(tmp_array);

  // if (!ThisTask)
  //   printf("[%s] The snapshot will contain a total of %Ld particles (true total number: %Ld, duplication factor: %f percent)\n",
  // 	   fdate(), myNtotal, MyGrids[0].Ntotal, 100.*(double)(myNtotal-MyGrids[0].Ntotal)/(double)MyGrids[0].Ntotal);

  /* The collector task opens the file and writes the header */
  if (ThisTask==collector)
    {
      strcpy(my_filename,filename);
      if (params.NumFiles>1)
	{
	  sprintf(tag,".%d",ThisFile);
	  strcat(my_filename,tag);
	}

      if (!ThisTask)
	printf("[%s] Task 0 will write %d particles out of %Ld in snapshot file %s\n",fdate(), NPartInFile, MyGrids[0].Ntotal, my_filename);

      if ( (file=fopen(my_filename,"w"))==0x0)
	{
	  printf("Error on Task %d: cannot open file %s\n",ThisTask,my_filename);
	  return 1;
	}

      /* writes information on the snapshot header */
      memset(&Header, 0, sizeof(SnapshotHeader_data));

      Header.NPart[1]=NPartInFile;
      Header.Mass[1]=params.ParticleMass*params.Hubble100*1e-10;
      Header.NPartTotal[1]=(unsigned int)( (myNtotal<<32) >>32 );
      Header.npartTotalHighWord[1]=(unsigned int)(myNtotal>>32);
      Header.Time=1./(1+outputs.z[myiout]);
      Header.RedShift=outputs.z[myiout];
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

    }

  return 0;
}


int write_block(Block_data block)
{
  /* writes a specified block in the snapshot */
  int npart, next, itask;
  unsigned long dummy_4;
  unsigned long long dummy_8;

  if (ThisTask==collector)
    {
      dummy_8=NPartInFile*block.sizeof_type;
      if (dummy_8>largest32-10)
	dummy_4=0;
      else
	dummy_4=(unsigned long)dummy_8;
      WriteBlockName(file,dummy_8,block.name);
      fwrite(&dummy_4, sizeof(dummy_4), 1, file);
      fwrite(block.data, block.sizeof_type, myNpart, file);
    }

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  npart=0;
	  MPI_Recv(&npart, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  MPI_Recv(block.data, npart*block.sizeof_type, MPI_BYTE, itask, 0, MPI_COMM_WORLD, &status);
	  fwrite(block.data, block.sizeof_type, npart, file);
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&myNpart, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  MPI_Send(block.data, myNpart*block.sizeof_type, MPI_BYTE, collector, 0, MPI_COMM_WORLD);
	}
    }

  /* collector task closes the block */
  if (ThisTask==collector)
    fwrite(&dummy_4, sizeof(dummy_4), 1, file);

  return 0;
}

void free_block(Block_data block)
{

  free(block.data);

}

int add_to_info(Block_data block)
{
  /* adds a block to the structure needed to write the INFO block */
  if (ThisTask!=collector)
    return 0;

  if (NextBlock>=NBlocks)
    {
      if (!ThisTask)
	printf("ERROR: I am trying to write the %d N-th block but NBLOCKS=%d\n",
	       NextBlock+1,NBlocks);
      return 1;
    }

  my_strcpy(InfoBlock[NextBlock].name,block.name,4);
  my_strcpy(InfoBlock[NextBlock].type,block.type,8);
  InfoBlock[NextBlock].ndim=block.ndim;
  /* pinocchio snapshots have only one particle type */
  for (int i=0; i<6; i++)
    InfoBlock[NextBlock].active[i]=(i==1);

  NextBlock++;

  return 0;
}


int write_info_block(void)
{
  /* writes the INFO block */
  int dummy,tb;
  char wn[5],wt[9];

  if (ThisTask!=collector)
    return 0;

  for (tb=0; tb<NBlocks; tb++)
    {

      if (!ThisTask)
 	{
	  my_strcpy(wn,InfoBlock[tb].name,4);
	  wn[4]='\0';
	  my_strcpy(wt,InfoBlock[tb].type,8);
	  wt[8]='\0';
	  printf("Written block: name=%s, type=%s, ndim=%d, active=[%d,%d,%d,%d,%d,%d]\n",
		 wn,wt,
		 InfoBlock[tb].ndim,
		 InfoBlock[tb].active[0],
		 InfoBlock[tb].active[1],
		 InfoBlock[tb].active[2],
		 InfoBlock[tb].active[3],
		 InfoBlock[tb].active[4],
		 InfoBlock[tb].active[5]);
	}
    }

  dummy=NBlocks*40;
  WriteBlockName(file,dummy,"INFO");
  fwrite(&dummy, sizeof(dummy), 1, file);
  for (tb=0; tb<NBlocks; tb++)
    fwrite(&InfoBlock[tb], sizeof(char), 40, file);
  fwrite(&dummy, sizeof(dummy), 1, file);

  return 0;
}


int initialize_ID(Block_data *block)
{
  /* this routine initializes the ID block */

  my_strcpy(block->name,"ID  ",4);
#ifdef LONGIDS
  my_strcpy(block->type,"LLONG   ",8);
  block->sizeof_type = sizeof(unsigned long long);
#else
  my_strcpy(block->type,"LONG    ",8);
  block->sizeof_type = sizeof(unsigned int);
#endif
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(MYIDTYPE));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* loop on all particles */
  for (int local_x = 0; local_x < MyGrids[0].GSlocal[_x_]; local_x++)
    {
      int idx_x = local_x * MyGrids[0].GSlocal[_y_];
      MYIDTYPE IDX_x = (local_x + MyGrids[0].GSstart[_x_]) * MyGrids[0].GSglobal[_y_];
      for (int local_y = 0; local_y < MyGrids[0].GSlocal[_y_]; local_y++)
	{
	  int idx_y = (idx_x + local_y) * MyGrids[0].GSlocal[_z_];
	  MYIDTYPE IDX_y = (IDX_x + local_y + MyGrids[0].GSstart[_y_]) * MyGrids[0].GSglobal[_z_];
	  for (int local_z = 0; local_z < MyGrids[0].GSlocal[_z_]; local_z++)
	    {

	      /* NB this is coherent with INDEX_TO_COORD */
	      ((MYIDTYPE*)block->data)[idx_y + local_z] = 1 + IDX_y + local_z + MyGrids[0].GSstart[_z_];

	    }
	}
    }

  return 0;
}


int initialize_density(int ThisGrid, Block_data *block)
{
  /* this routine initializes the density */

  my_strcpy(block->name,"DENS",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(float);
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(float));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[ThisGrid].total_local_size; i++)
    ((float*)block->data)[i] = density[ThisGrid][i];

  return 0;
}

int initialize_FMAX(Block_data *block)
{
  /* this routine initializes the FMAX block */

  my_strcpy(block->name,"FMAX",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(float);
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(float));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    ((float*)block->data)[i] = products[i].Fmax;

  return 0;
}

int initialize_RMAX(Block_data *block)
{
  /* this routine initializes the RMAX block */

  my_strcpy(block->name,"RMAX",4);
  my_strcpy(block->type,"LONG    ",8);
  block->sizeof_type = sizeof(int);
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(int));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    ((int*)block->data)[i] = products[i].Rmax;

  return 0;
}

int initialize_ZEL(Block_data *block)
{
  /* this routine initializes the ZEL block */

  my_strcpy(block->name,"ZEL ",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      ((AuxStruct*)block->data)[i].axis[0] = products[i].Vel[0];
      ((AuxStruct*)block->data)[i].axis[1] = products[i].Vel[1];
      ((AuxStruct*)block->data)[i].axis[2] = products[i].Vel[2];
    }

  return 0;
}

#ifdef TWO_LPT
int initialize_2LPT(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"2LPT",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      ((AuxStruct*)block->data)[i].axis[0] = products[i].Vel_2LPT[0];
      ((AuxStruct*)block->data)[i].axis[1] = products[i].Vel_2LPT[1];
      ((AuxStruct*)block->data)[i].axis[2] = products[i].Vel_2LPT[2];
    }

  return 0;
}

#ifdef THREE_LPT
int initialize_3LPT_1(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"31PT",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      ((AuxStruct*)block->data)[i].axis[0] = products[i].Vel_3LPT_1[0];
      ((AuxStruct*)block->data)[i].axis[1] = products[i].Vel_3LPT_1[1];
      ((AuxStruct*)block->data)[i].axis[2] = products[i].Vel_3LPT_1[2];
    }

  return 0;
}

int initialize_3LPT_2(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"32PT",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      ((AuxStruct*)block->data)[i].axis[0] = products[i].Vel_3LPT_2[0];
      ((AuxStruct*)block->data)[i].axis[1] = products[i].Vel_3LPT_2[1];
      ((AuxStruct*)block->data)[i].axis[2] = products[i].Vel_3LPT_2[2];
    }

  return 0;
}
#endif
#endif

int initialize_ZACC(Block_data *block)
{
  /* this routine initializes the ZACC block */

  my_strcpy(block->name,"ZACC",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(float);
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(float));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    ((float*)block->data)[i] = products[i].zacc;

  return 0;
}

int initialize_GRUP(Block_data *block)
{
  /* this routine initializes the GRUOP_ID block */

  my_strcpy(block->name,"GRUP",4);
  my_strcpy(block->type,"LONG    ",8);
  block->sizeof_type = sizeof(int);
  block->ndim=1;

  /* each task builds its catalogue: IDs */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(int));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    ((int*)block->data)[i] = products[i].group_ID;

  return 0;
}

int initialize_POS(Block_data *block)
{
  /* initializes the POS block by moving all particles with LPT */
  int ibox, jbox, kbox;

  my_strcpy(block->name,"POS ",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue: Pos */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      INDEX_TO_COORD(i, ibox, jbox, kbox, MyGrids[0].GSlocal);
      int q[3]={ibox+MyGrids[0].GSstart[0],jbox+MyGrids[0].GSstart[1],kbox+MyGrids[0].GSstart[2]};

      for (int j=0; j<3; j++)
	{
	  PRODFLOAT pos = q[j] + SHIFT + products[i].Vel[j]; // * myobj.D;
#ifdef TWO_LPT
	  pos += products[i].Vel_2LPT[j]; //  * myobj.D2;
#ifdef THREE_LPT
	  pos += products[i].Vel_3LPT_1[j] // * myobj.D31 
	    + products[i].Vel_3LPT_2[j]; // * myobj.D32;
#endif
#endif
	  if (pos>=MyGrids[0].GSglobal[j]) pos-=MyGrids[0].GSglobal[j];
	  if (pos< 0.0)                    pos+=MyGrids[0].GSglobal[j];

	  ((AuxStruct*)block->data)[i].axis[j] = (float)(pos * params.InterPartDist * params.Hubble100)
#ifdef POS_IN_KPC
	    * 1000.
#endif
	    ;
	}
    }

  return 0;

}


int initialize_VEL(Block_data *block)
{

  /* initializes the VEL block by moving all particles with LPT */

  my_strcpy(block->name,"VEL ",4);
  my_strcpy(block->type,"FLOATN  ",8);
  block->sizeof_type = sizeof(AuxStruct);
  block->ndim=3;

  /* each task builds its catalogue: Vel */
  block->data = (void *)malloc(Npart_array[collector] * sizeof(AuxStruct));
  if (block->data==0x0)
    {
      printf("ERROR on Task %d: could not allocate block to be written in snapshot\n",ThisTask);
      return 1;
    }

  /* Warning: the scale for fomega should be set somehow... (very little impact) */
  double vfact=Hubble(redshift)/(1.+redshift)*params.InterPartDist * (1+redshift) * fomega(redshift, params.k_for_GM)/sqrt(1.+outputs.z[myiout]);

  /* loop on all particles */
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    {
      for (int j=0; j<3; j++)
	{
	    PRODFLOAT vv = products[i].Vel[j]; //* myobj.Dv;
#ifdef TWO_LPT
	    vv += products[i].Vel_2LPT[j];// * myobj.D2v;
#ifdef THREE_LPT
	    vv += products[i].Vel_3LPT_1[j]  // * myobj.D31v
	      +   products[i].Vel_3LPT_1[j]; // * myobj.D32v;
#endif
#endif
	  ((AuxStruct*)block->data)[i].axis[j]=(float)(vv * vfact);
	}
    }

  return 0;
}


void WriteBlockName(FILE *fd,unsigned long long dummy, char* LABEL)
{

  if (dummy>largest32-10)
    {
      unsigned long dummy_4;
      unsigned long long dummy2_8;
      dummy_4=12;
      fwrite(&dummy_4,sizeof(dummy_4),1,fd);
      fwrite(LABEL,4,1,fd);
      dummy2_8=dummy+8;
      fwrite(&dummy2_8,sizeof(dummy2_8),1,fd);
      fwrite(&dummy_4,sizeof(dummy_4),1,fd);
    }
  else
    {
      unsigned int dummy2;
      dummy2=8;
      fwrite(&dummy2,sizeof(dummy2),1,fd);
      fwrite(LABEL,4,1,fd);
      dummy2=(unsigned long)dummy+8;
      fwrite(&dummy2,sizeof(dummy2),1,fd);
      dummy2=8;
      fwrite(&dummy2,sizeof(dummy2),1,fd);
    }
} 


void my_strcpy(char *to, char *from, int n)
{
  int i;
  for (i=0; i<n; i++)
    *(to+i)=*(from+i);
}


#endif
