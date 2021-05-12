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


//#define LONGIDS
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

typedef struct
{
  char name[4], type[8];
  int ndim, active[6];
  size_t sizeof_type;
  void *data;
} Block_data;

Block_data *InfoBlock;
int NBlocks, NextBlock;

typedef struct
{
  float axis[3];
} AuxStruct;

MPI_Status status;

void WriteBlockName(FILE *,int, char*);
void my_strcpy(char *,char *, int);
int write_header(int);
int write_block(Block_data);
void free_block(Block_data);
int add_to_info(Block_data);
int write_info_block(void);
int count_particles_for_final_snapshot(void);
int initialize_ID_subvolumes(Block_data*);
int initialize_POS_withNFW(Block_data*);
int initialize_VEL_withNFW(Block_data*);
int initialize_ID_fftspace(Block_data*);
int initialize_FMAX_fftspace(Block_data*);
int initialize_RMAX_fftspace(Block_data*);
int initialize_VEL_fftspace(Block_data*);
#ifdef TWO_LPT
int initialize_VEL_2LPT_fftspace(Block_data*);
#ifdef THREE_LPT
int initialize_VEL_3LPT_1_fftspace(Block_data*);
int initialize_VEL_3LPT_2_fftspace(Block_data*);
#endif
#endif
int initialize_density(int, Block_data*);

#ifdef ROTATE_BOX
static int rot[3]={1,2,0};
#else
static int rot[3]={0,1,2};
#endif

/* In this routine each task outputs all the particles belonging to good groups.
   If FORCE_PARTICLE_NUMBER is not active, some particles will be duplicated 
   but if FORCE_PARTICLE_NUMBER is active, some groups will have 
   fewer particles that their putative mass */

char filename[LBLENGTH];
FILE *file;
int NTasksPerFile,collector,ThisFile,NPartInFile,myNpart,myiout;
int *Npart_array;

int write_snapshot(int iout)
{
  /* writes displacements of all particles in a GADGET format 
     particles in halos are distributed around the center of mass
     as NFW spheres with virial velocity distribution */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.snapshot_%03d",params.RunFlag,iout);

  int my_snap_type=0;
  myiout=iout;

  /* allocates structure to handle the INFO block */
  NBlocks=3;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header(my_snap_type))
    return 1;

  /* writing of IDs */
  if (initialize_ID_subvolumes(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of POS */
  if (initialize_POS_withNFW(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_VEL_withNFW(&block))
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


int write_products()
{
  /* writes the fmax products as in a snapshot format */

  Block_data block;

  /* Snapshot filename */
  sprintf(filename,"pinocchio.%s.products.out",params.RunFlag);

  int my_snap_type=1;

  /* allocates structure to handle the INFO block */
#ifdef TWO_LPT
#ifdef THREE_LPT
  NBlocks=7;
#else
  NBlocks=5;
#endif
#else
  NBlocks=4;
#endif
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing products in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header(my_snap_type))
    return 1;

  /* writing of IDs */
  if (initialize_ID_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of RMAX */
  if (initialize_RMAX_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of FMAX */
  if (initialize_FMAX_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VEL */
  if (initialize_VEL_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef TWO_LPT
  /* writing of VEL2 */
  if (initialize_VEL_2LPT_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

#ifdef THREE_LPT
  /* writing of VL31 */
  if (initialize_VEL_3LPT_1_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);

  /* writing of VL32 */
  if (initialize_VEL_3LPT_2_fftspace(&block))
    return 1;
  if (add_to_info(block))
    return 1;
  if (write_block(block))
    return 1;
  free_block(block);
#endif
#endif

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
  sprintf(filename,"pinocchio.%s.density.out",params.RunFlag);

  int my_snap_type=1;

  /* allocates structure to handle the INFO block */
  NBlocks=2;
  NextBlock=0;
  InfoBlock=(Block_data *)calloc(NBlocks, sizeof(Block_data));

  if (!ThisTask)
    printf("[%s] Writing density in snapshot file %s\n",fdate(),filename);

  /* this routine opens the file */
  if (write_header(my_snap_type))
    return 1;

  /* writing of IDs */
  if (initialize_ID_fftspace(&block))
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


int write_header(int snap_type)
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

  switch (snap_type)
    {
    case 0: /* particles live in the subvolume domain */
      myNpart = count_particles_for_final_snapshot(); // SEPARARE FORCE_PARTICLE_NUMBER QUI
      break;
    case 1: /* particles live in the fft domain */
      myNpart = MyGrids[0].total_local_size; // NON E` VERO!!!
      break;
    default:
      myNpart = 0;
      break;
    }

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

  if (!ThisTask)
    printf("[%s] The snapshot will contain a total of %Ld particles (true total number: %Ld, duplication factor: %f percent)\n",
	   fdate(), myNtotal, MyGrids[0].Ntotal, 100.*(double)(myNtotal-MyGrids[0].Ntotal)/(double)MyGrids[0].Ntotal);

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
	printf("[%s] Task 0 will write %d particles in snapshot file %s\n",fdate(),NPartInFile,my_filename);

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
  int dummy, npart, next, itask;

  if (ThisTask==collector)
    {
      dummy=NPartInFile*block.sizeof_type;
      WriteBlockName(file,dummy,block.name);
      fwrite(&dummy, sizeof(dummy), 1, file);
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
    fwrite(&dummy, sizeof(dummy), 1, file);

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


int count_particles_for_final_snapshot(void)
{
  /* This function counts the number of particles that will be written by ThisTask.
     These live in the subbox space, each task outputs all particles that are assigned
     to its halos, even if they live in different sub-volumes. Note that this 
     does not guarantee particle conservation.
   */

  int Lgridxy, myN, i, ibox,jbox,kbox,kk, good_particle;

  /* this is relevant when the particles live in the sub-volume domain */
  Lgridxy = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_];

  /* Each task counts the number of particles that will be output */
  /* first loop is on good particles that are not in groups */
  for (i=0, myN=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl[_x_];
      ibox=kk-jbox*subbox.Lgwbl[_x_];
      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER) // PROBABILMENTE SI PUO` TOGLIERE
	  && group_ID[i]<=FILAMENT
#endif
	  )
	myN++;
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      myN+=groups[i].Mass;
#endif

  return myN;
}


int initialize_ID_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

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
  for (int i=0; i < MyGrids[0].total_local_size; i++)
    ((float*)block->data)[i] = density[ThisGrid][i];

  return 0;
}

int initialize_FMAX_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

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

int initialize_RMAX_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

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

int initialize_VEL_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"VEL ",4);
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
int initialize_VEL_2LPT_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"VEL2",4);
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
int initialize_VEL_3LPT_1_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"VL31",4);
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
int initialize_VEL_3LPT_2_fftspace(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */

  my_strcpy(block->name,"VL32",4);
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

// DA QUI IN POI E` TUTTO DA CAPIRE COME FARE

int initialize_ID_subvolumes(Block_data *block)
{
  /* this routine initializes the ID block when particles live in sub-volume domain */
  int index, i, ibox, jbox, kbox, kk, good_particle, global[3], Lgridxy;

  Lgridxy = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_];

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

  /* first loop is on good particles that are not in groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(i,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[i],ibox,jbox,kbox,subbox.Lgwbl);
#endif

      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  ((MYIDTYPE*)block->data)[index]=1+
	    COORD_TO_INDEX((long long)((ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]),
			   (long long)((jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]),
			   (long long)((kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]),
			   MyGrids[0].GSglobal);
	  ++index;
	}
    }

#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
  /* second loop is on good groups */
  int next;
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point >= 0 && groups[i].good)
      {
	next=groups[i].point;
	for (int npart=0; npart<groups[i].Mass; npart++)
	  {
	    kbox=next/Lgridxy;
	    kk=next-kbox*Lgridxy;
	    jbox=kk/subbox.Lgwbl[_x_];
	    ibox=kk-jbox*subbox.Lgwbl[_x_];
	    /* particle coordinates in the box, imposing PBCs */
	    global[_x_] = ibox + subbox.stabl[_x_];
	    if (global[_x_]<0) global[_x_]+=MyGrids[0].GSglobal[_x_];
	    if (global[_x_]>=MyGrids[0].GSglobal[_x_]) global[_x_]-=MyGrids[0].GSglobal[_x_];
	    global[_y_] = jbox + subbox.stabl[_y_];
	    if (global[_y_]<0) global[_y_]+=MyGrids[0].GSglobal[_y_];
	    if (global[_y_]>=MyGrids[0].GSglobal[_y_]) global[_y_]-=MyGrids[0].GSglobal[_y_];
	    global[_z_] = kbox + subbox.stabl[_z_];
	    if (global[_z_]<0) global[_z_]+=MyGrids[0].GSglobal[_z_];
	    if (global[_z_]>=MyGrids[0].GSglobal[_z_]) global[_z_]-=MyGrids[0].GSglobal[_z_];
#ifdef ROTATE_BOX
	    ((MYIDTYPE*)block->data)[index]=1+(MYIDTYPE)global[_x_] + 
	      ( (MYIDTYPE)global[_z_] + 
		(MYIDTYPE)global[_y_] * (MYIDTYPE)MyGrids[0].GSglobal[_z_] ) * 
	      (MYIDTYPE)MyGrids[0].GSglobal[_x_];
#else
	    ((MYIDTYPE*)block->data)[index]=1+(MYIDTYPE)global[_x_] + 
	      ( (MYIDTYPE)global[_y_] + 
		(MYIDTYPE)global[_z_] * (MYIDTYPE)MyGrids[0].GSglobal[_y_] ) * 
	      (MYIDTYPE)MyGrids[0].GSglobal[_x_];
#endif
	    ++index;
	    next=linking_list[next];
	  }
      }

#endif

  return 0;
}


int initialize_POS_withNFW(Block_data *block)
{
  /* */
  int index, i, ibox, jbox, kbox, kk, good_particle, global[3], j, part;
  double MyPos[3], conc, rvir, rnd, nfwfac, area, xrnd, yrnd, probfunc, u, theta;

  int Lgridxy = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_];

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

  double Dz = GrowingMode(outputs.z[myiout],params.k_for_GM);

  /* the first loop is on particles outside groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl[_x_];
      ibox=kk-jbox*subbox.Lgwbl[_x_];
      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global[_x_] = ibox + subbox.stabl[_x_];
	  if (global[_x_]<0) global[_x_]+=MyGrids[0].GSglobal[_x_];
	  if (global[_x_]>=MyGrids[0].GSglobal[_x_]) global[_x_]-=MyGrids[0].GSglobal[_x_];
	  global[_y_] = jbox + subbox.stabl[_y_];
	  if (global[_y_]<0) global[_y_]+=MyGrids[0].GSglobal[_y_];
	  if (global[_y_]>=MyGrids[0].GSglobal[_y_]) global[_y_]-=MyGrids[0].GSglobal[_y_];
	  global[_z_] = kbox + subbox.stabl[_z_];
	  if (global[_z_]<0) global[_z_]+=MyGrids[0].GSglobal[_z_];
	  if (global[_z_]>=MyGrids[0].GSglobal[_z_]) global[_z_]-=MyGrids[0].GSglobal[_z_];

	  /* particles outside groups (uncollapsed and filament particles)
	     are displaced with LPT */
#ifdef FORCE_PARTICLE_NUMBER
	  if (group_ID[i]<=FILAMENT)
	    {
#endif
	      set_point(global[_x_],global[_y_],global[_z_],i,outputs.F[myiout],&obj);
	      for (j=0; j<3; j++)
		{
		  ((AuxStruct*)block->data)[index].axis[j]=
		    (float)(q2x(rot[j], &obj, ORDER_FOR_CATALOG) * 
			    params.InterPartDist * params.Hubble100);
		  if (((AuxStruct*)block->data)[index].axis[j] >= params.BoxSize_h100)
		    ((AuxStruct*)block->data)[index].axis[j] -= (float)params.BoxSize_h100;
		  if (((AuxStruct*)block->data)[index].axis[j] < 0.0) 
		    ((AuxStruct*)block->data)[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
		  ((AuxStruct*)block->data)[index].axis[j] *= 1000.;
#endif
		}
#ifdef FORCE_PARTICLE_NUMBER
	    }
	  else
	    {
	      set_obj(group_ID[i], outputs.F[myiout], &obj);
	      for (j=0; j<3; j++)
		MyPos[j]=(float)((q2x(rot[j], &obj, ORDER_FOR_CATALOG) + subbox.stabl[rot[j]]) *
				 params.InterPartDist * params.Hubble100);
	      /* particles in the group are distributed as NFW */
	      /* Concentration taken from Bhattacharya, et al. 2013 */
	      conc = pow(Dz,0.54) * 5.9
		* pow( (1.12*pow(groups[group_ID[i]].Mass*params.ParticleMass*params.Hubble100 
				 /5.e13,0.3) + 0.53)/Dz , -0.35);
	      rvir = pow(0.01 * GRAVITY * groups[group_ID[i]].Mass*params.ParticleMass /
		   pow(Hubble(outputs.z[myiout]),2.0), 1./3.) * (1. + outputs.z[myiout]) * params.Hubble100;

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

		      ((AuxStruct*)block->data)[index].axis[_x_] = MyPos[_x_] + (float)(xrnd * rvir * sqrt(1.-u*u)*cos(theta));
		      ((AuxStruct*)block->data)[index].axis[_y_] = MyPos[_y_] + (float)(xrnd * rvir * sqrt(1.-u*u)*sin(theta));
		      ((AuxStruct*)block->data)[index].axis[_z_] = MyPos[_z_] + (float)(xrnd * rvir * u);

		      for (j=0; j<3; j++)
			{
			  if (((AuxStruct*)block->data)[index].axis[j] >= params.BoxSize_h100) 
			    ((AuxStruct*)block->data)[index].axis[j] -= (float)params.BoxSize_h100;
			  if (((AuxStruct*)block->data)[index].axis[j] < 0.0) 
			    ((AuxStruct*)block->data)[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
			  ((AuxStruct*)block->data)[index].axis[j] *= 1000.;
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
	set_obj(i, outputs.F[myiout], &obj);
	for (j=0; j<3; j++)
	  MyPos[j]=(float)((q2x(rot[j], &obj, ORDER_FOR_CATALOG) + subbox.stabl[rot[j]]) *
			   params.InterPartDist * params.Hubble100);

	/* virial radius and concentration of the group */
	/* Concentration taken from Bhattacharya, et al. 2013 */
	conc = pow(Dz,0.54) * 5.9
	  * pow( (1.12*pow(groups[i].Mass*params.ParticleMass*params.Hubble100 /5.e13,0.3)+ 0.53)/Dz ,
		 -0.35);
	rvir = pow(0.01 * GRAVITY * groups[i].Mass*params.ParticleMass /
		   pow(Hubble(outputs.z[myiout]),2.0), 1./3.) * (1. + outputs.z[myiout]) * params.Hubble100;

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

		((AuxStruct*)block->data)[index].axis[_x_] = MyPos[_x_] + (float)(xrnd * rvir * sqrt(1.-u*u)*cos(theta));
		((AuxStruct*)block->data)[index].axis[_y_] = MyPos[_y_] + (float)(xrnd * rvir * sqrt(1.-u*u)*sin(theta));
		((AuxStruct*)block->data)[index].axis[_z_] = MyPos[_z_] + (float)(xrnd * rvir * u);

		for (j=0; j<3; j++)
		  {
		    if (((AuxStruct*)block->data)[index].axis[j] >= params.BoxSize_h100) 
		      ((AuxStruct*)block->data)[index].axis[j] -= (float)params.BoxSize_h100;
		    if (((AuxStruct*)block->data)[index].axis[j] < 0.0) 
		      ((AuxStruct*)block->data)[index].axis[j] += (float)params.BoxSize_h100;
#ifdef POS_IN_KPC
		    ((AuxStruct*)block->data)[index].axis[j] *= 1000.;
#endif
		  }
		index++;
		part++;
	      }
	  }
	while (part<groups[i].Mass);
      }
#endif

  return 0;

}


int initialize_VEL_withNFW(Block_data *block)
{

  int index, i, ibox, jbox, kbox, kk, good_particle, global[3], j;
  double MyVel[3], vfact, sigma, rvir;

  int Lgridxy = subbox.Lgwbl[_x_] * subbox.Lgwbl[_y_];

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

  vfact=sqrt(1.+outputs.z[myiout]);
  /* the first loop is on particles outside groups */
  for (i=0, index=0; i<subbox.Npart; i++)
    {
      /* grid coordinates from the indices (sub-box coordinates) */
      kbox=i/Lgridxy;
      kk=i-kbox*Lgridxy;
      jbox=kk/subbox.Lgwbl[_x_];
      ibox=kk-jbox*subbox.Lgwbl[_x_];
      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );
      if (good_particle
#if !defined(ONLY_LPT_DISPLACEMENTS) && !defined(FORCE_PARTICLE_NUMBER)
	  && group_ID[i]<=FILAMENT
#endif
	  )
	{
	  /* particle coordinates in the box, imposing PBCs */
	  global[_x_] = ibox + subbox.stabl[_x_];
	  if (global[_x_]<0) global[_x_]+=MyGrids[0].GSglobal[_x_];
	  if (global[_x_]>=MyGrids[0].GSglobal[_x_]) global[_x_]-=MyGrids[0].GSglobal[_x_];
	  global[_y_] = jbox + subbox.stabl[_y_];
	  if (global[_y_]<0) global[_y_]+=MyGrids[0].GSglobal[_y_];
	  if (global[_y_]>=MyGrids[0].GSglobal[_y_]) global[_y_]-=MyGrids[0].GSglobal[_y_];
	  global[_z_] = kbox + subbox.stabl[_z_];
	  if (global[_z_]<0) global[_z_]+=MyGrids[0].GSglobal[_z_];
	  if (global[_z_]>=MyGrids[0].GSglobal[_z_]) global[_z_]-=MyGrids[0].GSglobal[_z_];
	  
	  /* particles outside groups (uncollapsed and filament particles)
	     are displaced with LPT */
#ifdef FORCE_PARTICLE_NUMBER
	  if (group_ID[i]<=FILAMENT)
	    {
#endif
	      set_point(global[_x_],global[_y_],global[_z_],i,outputs.F[myiout],&obj);
	      for (j=0; j<3; j++)
		/* GADGET format requires velocities to be divided by sqrt(a) */
		((AuxStruct*)block->data)[index].axis[j]=(float)(vel(rot[j], &obj)*vfact);
#ifdef FORCE_PARTICLE_NUMBER
	    }
	  else
	    {
	      /* the center of mass of groups is displaced with LPT */
	      set_obj(group_ID[i], outputs.F[myiout], &obj);
	      set_obj_vel(group_ID[i], outputs.F[myiout], &obj);
	      for (j=0; j<3; j++)
		MyVel[j]=(float)(vel(rot[j], &obj)*vfact);

	      rvir = pow(0.01 * GRAVITY * groups[group_ID[i]].Mass*params.ParticleMass /
		   pow(Hubble(outputs.z[myiout]),2.0), 1./3.) * (1. + outputs.z[myiout]) * params.Hubble100;
	      sigma = sqrt(GRAVITY * groups[group_ID[i]].Mass * params.ParticleMass / 3./ rvir);

	      for (j=0; j<3; j++)
		((AuxStruct*)block->data)[index].axis[j]=MyVel[j] + gsl_ran_gaussian(random_generator, sigma) * vfact;
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
	set_obj(i, outputs.F[myiout], &obj);
	set_obj_vel(i, outputs.F[myiout], &obj);
	for (j=0; j<3; j++)
	  MyVel[j]=(float)(vel(rot[j], &obj)*vfact);

	/* /\* virial radius and concentration of the group *\/ */
	/* Dz = GrowingMode(outputs.z[myiout],params.k_for_GM); */
	/* /\* this is the virial radius in physical true kpc *\/ */
	/* rvir = pow(0.01 * GRAVITY * groups[i].Mass*params.ParticleMass / */
	/* 	   pow(Hubble(outputs.z[myiout]),2.0), 1./3.);   */
	sigma = sqrt(GRAVITY * groups[i].Mass * params.ParticleMass / 3./ rvir);

	for (npart=0; npart<groups[i].Mass; npart++)
	  {
	    for (j=0; j<3; j++)
	      ((AuxStruct*)block->data)[index].axis[j]=MyVel[j] + gsl_ran_gaussian(random_generator, sigma) * vfact;
	    ++index;
	  }
      }
#endif

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
