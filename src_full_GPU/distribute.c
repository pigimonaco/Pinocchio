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

#define BUFLEN 50000

//#define DEBUG

static product_data *comm_buffer;
#ifndef CLASSIC_FRAGMENTATION
static unsigned long long frag_offset;
#endif

int intersection(int *, int *, int *);
int send_data(int *, int);
int recv_data(int *, int);
int keep_data(int *, int *);
unsigned int fft_space_index(unsigned int, int *);
unsigned int subbox_space_index(unsigned int, int *);
#ifndef CLASSIC_FRAGMENTATION
int get_distmap_bit(unsigned int *, unsigned int);
void set_distmap_bit(unsigned int *, unsigned int, int);
void build_distmap(unsigned int *, int *);
void update_distmap(unsigned int *, int *);
#endif

#ifdef DEBUG
FILE *DBGFD;
#endif

int distribute(void)
{
  /* Distributes products from fft-space to sub-volumes */

  int my_fft_box[6],my_subbox[6];
  int log_ntask, bit, receiver, sender;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* this defines the box that the task possesses in the subbox space 
     (the starting coordinate may be negative) */
  my_subbox[0] = subbox.stabl[_x_];
  my_subbox[1] = subbox.stabl[_y_];
  my_subbox[2] = subbox.stabl[_z_];
  my_subbox[3] = subbox.Lgwbl[_x_];
  my_subbox[4] = subbox.Lgwbl[_y_];
  my_subbox[5] = subbox.Lgwbl[_z_];

#ifdef DEBUG
  char fname[SBLENGTH];
  sprintf(fname,"Task%d.dbg",ThisTask);
  DBGFD = fopen(fname,"a");
#endif

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* these are the already stored particles */
#ifndef CLASSIC_FRAGMENTATION
  frag_offset=subbox.Nneeded;
#endif

  /* stores the data relative to the intersection of the fft box and subbox */
  if (keep_data(my_fft_box, my_subbox))
    return 1;

  comm_buffer=(product_data*)calloc(BUFLEN , sizeof(product_data));
  if (comm_buffer==0x0)
    {
      printf("ERROR on task %d: could not allocate comm_buffer in distribute\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  /* hypercubic communication scheme */
  for (log_ntask=0; log_ntask<1000; log_ntask++)
    if (1<<log_ntask >= NTasks)
      break;

  /* loop on hypercube dimension */
  for (bit=1; bit<1<<log_ntask; bit++)
    {
    /* loop on tasks */
      for (sender=0; sender<NTasks; sender++)
	{
	  /* receiver task is computed with a bitwise xor */
	  receiver = sender ^ bit;

	  /* condition on sender and receiver */
	  if (receiver < NTasks && sender < receiver)
	    {

	      /* the communication will be done 
		 first sender -> receiver and then receiver -> sender */

	      if (ThisTask==sender)
		send_data(my_fft_box, receiver);
	      else if (ThisTask==receiver)
		{
		  if (recv_data(my_subbox, sender))
		    return 1;
		}

	      if (ThisTask==receiver)
		send_data(my_fft_box, sender);
	      else if (ThisTask==sender)
		{
		  if (recv_data(my_subbox, receiver))
		    return 1;
		}

	    }
	}
    }


  /* updates the number of stored particles */
#ifdef CLASSIC_FRAGMENTATION
  subbox.Nstored=subbox.Npart;
#else

  if (frag_offset > subbox.Nalloc)
    subbox.Nstored=subbox.Nalloc;
  else
    subbox.Nstored=frag_offset;

  subbox.Nneeded=frag_offset;

#endif

  free(comm_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
  fclose(DBGFD);
#endif

  return 0;
}


int intersection(int *fbox, int *sbox, int *ibox)
{

  unsigned int istart[3], istop[3], istart2[3], istop2[3], dim, Nint,
    stop1, stop2, this, off;
  unsigned int ISTARTx, ISTOPx, ISTARTy, ISTOPy, ISTARTz, ISTOPz, ax, ay, az;

  /* intersection of the two boxes is considered dimension by dimension */
  Nint=1;
  for (dim=0; dim<3; dim++)
    {
      /* in case, fix negative starting point */
      if (sbox[dim]<0)
	sbox[dim]+=MyGrids[0].GSglobal[dim];

      /* first stopping point for subbox is at most the global box edge */
      stop1=fbox[dim]+fbox[dim+3];
      if (sbox[dim]+sbox[dim+3] > MyGrids[0].GSglobal[dim])
	stop2=MyGrids[0].GSglobal[dim];
      else
	stop2=sbox[dim]+sbox[dim+3];

      /* intersection up to the global box edge */
      istart[dim] = ( fbox[dim] > sbox[dim] ? fbox[dim] : sbox[dim] );
      istop[dim] = ( stop1 < stop2 ? stop1 : stop2 );

      /* if the subbox goes beyond the global box edge, 
	 apply PBCs to the other segment and check the intersection */
      if ( (stop2=sbox[dim]+sbox[dim+3]) > MyGrids[0].GSglobal[dim] )
	{
	  stop2 = stop2%MyGrids[0].GSglobal[dim];
	  istart2[dim] = ( fbox[dim] > 0 ? fbox[dim] : 0);
	  istop2[dim] = ( stop1 < stop2 ? stop1 : stop2 );
	}
      else
	{
	  istart2[dim]=1;
	  istop2[dim]=0;
	}

      /* this dimension contributes 0, 1 or 2 */
      Nint *= (istart[dim] < istop[dim]) + (istart2[dim] < istop2[dim]);
    }

  /* store all intersections, looping on the two options for each dimension */
  if (Nint)
    {
      this=0;
      for (ax=0; ax<2; ax++)
	{
	  if (ax)
	    {
	      ISTARTx=istart[0];
	      ISTOPx=istop[0];
	    }
	  else
	    {
	      ISTARTx=istart2[0];
	      ISTOPx=istop2[0];
	    }
	  for (ay=0; ay<2; ay++)
	    {
	      if (ay)
		{
		  ISTARTy=istart[1];
		  ISTOPy=istop[1];
		}
	      else
		{
		  ISTARTy=istart2[1];
		  ISTOPy=istop2[1];
		}
	      for (az=0; az<2; az++)
		{
		  if (az)
		    {
		      ISTARTz=istart[2];
		      ISTOPz=istop[2];
		    }
		  else
		    {
		      ISTARTz=istart2[2];
		      ISTOPz=istop2[2];
		    }

		  if ( (ISTARTx<ISTOPx) & (ISTARTy<ISTOPy) & (ISTARTz<ISTOPz) )
		    {
		      off=this*6;
		      ibox[  off]=ISTARTx;
		      ibox[1+off]=ISTARTy;
		      ibox[2+off]=ISTARTz;
		      ibox[3+off]=ISTOPx-ISTARTx;
		      ibox[4+off]=ISTOPy-ISTARTy;
		      ibox[5+off]=ISTOPz-ISTARTz;
		      this++;
		    }
		}
	    }
	}
    }	  


#ifdef DEBUG
  fprintf(DBGFD,"INTERSECTION\n");
  fprintf(DBGFD,"Task %d, fft box:        %d %d %d   %d %d %d\n",
	 ThisTask,fbox[0],fbox[1],fbox[2],fbox[3],fbox[4],fbox[5]);
  fprintf(DBGFD,"         subvolume:      %d %d %d   %d %d %d\n",
	 sbox[0],sbox[1],sbox[2],sbox[3],sbox[4],sbox[5]);
  for (this=0; this<Nint; this++)
    {
      off=this*6;
      fprintf(DBGFD,"         intersection %d: %d %d %d   %d %d %d\n",
	     this,ibox[off],ibox[1+off],ibox[2+off],ibox[3+off],ibox[4+off],ibox[5+off]);
    }
  if (!Nint)
    fprintf(DBGFD,"        no intersections\n");
#endif

  return Nint;
}


int send_data(int *mybox, int target)
{
  /* This routine sends the content of a box to a target task.
     Communication is divided in these stages:
     1) the sender communicates the start and length of its box,
     2) the sender receives the start and the length of the needed box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the sender receives a map of needed particles
     5) the sender loops on particles and sends them in chunks of size BUFLEN
  */
  
#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
#endif
  unsigned int size, bufcount;
  int targetbox[6], interbox[48], i, Nint, box, off;
  MPI_Status status;

  MPI_Send(mybox    , 6, MPI_INT, target, 0, MPI_COMM_WORLD);
  MPI_Recv(targetbox, 6, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

  /* up to 8 intersections of the two boxes */
  Nint=intersection(mybox, targetbox, interbox);
  for (box=0; box<Nint; box++)
    {
      off=box*6;
      size = interbox[3+off] * interbox[4+off] * interbox[5+off];

#ifdef DEBUG
      fprintf(DBGFD,"Task %d will send this box: %d %d %d -- %d %d %d\n",
	     ThisTask, interbox[0+off],interbox[1+off],interbox[2+off],
	     interbox[3+off],interbox[4+off],interbox[5+off]);
#endif

#ifndef CLASSIC_FRAGMENTATION
      /* receive the map of needed particles */
      mapl = size/8 + 1;
      map = (unsigned int*)calloc(mapl, sizeof(unsigned int));

      /* the receiver has built its map and is sending it */
      MPI_Recv(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

      /* updates the map and sends it back to the target */
      update_distmap(map, interbox+off);
      MPI_Send(map, mapl, MPI_INT, target, 0, MPI_COMM_WORLD);
#endif

      /* load particles on the buffer and send them to the target */
      bufcount=0;
      for (i=0; i<size; i++)
	{
#ifndef CLASSIC_FRAGMENTATION
	  if (get_distmap_bit(map, i))
#endif
	    {
	      memcpy(&comm_buffer[bufcount],&products[fft_space_index(i, interbox+off)],sizeof(product_data));
	      //comm_buffer[bufcount] = products[fft_space_index(i, interbox+off)];
	      //*(comm_buffer + bufcount) = *(products + fft_space_index(i, interbox+off));

	      //fprintf(DBGFD," %d %d %d -- %f %f\n",i, fft_space_index(i, interbox+off),bufcount,
	      //        products[fft_space_index(i, interbox+off)].Fmax, comm_buffer[bufcount].Fmax); // LEVARE!!!

	      bufcount++;
	    }
	  if (bufcount==BUFLEN)
	    {
#ifdef DEBUG
	      fprintf(DBGFD,"...sending %d products to Task %d... %d\n",bufcount,target,(int)sizeof(product_data));
	      if (bufcount<10)
		{
		  for (int u=0; u<bufcount; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
	      else
		{
		  for (int u=0; u<5; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD," ... ");
		  for (int u=bufcount-5; u<bufcount; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
#endif
	      MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
	      MPI_Send(comm_buffer, bufcount * sizeof(product_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
	      bufcount=0;
	    }
	}

      if (bufcount)
	{
#ifdef DEBUG
	  fprintf(DBGFD,"...sending %d products to Task %d...\n",bufcount,target);
	  if (bufcount<10)
	    {
	      for (int u=0; u<bufcount; u++)
		fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
	      fprintf(DBGFD,"\n");
	    }
	  else
	    {
	      for (int u=0; u<5; u++)
		fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
	      fprintf(DBGFD," ... ");
	      for (int u=bufcount-5; u<bufcount; u++)
		fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
	      fprintf(DBGFD,"\n");
	    }
#endif
	  MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
	  MPI_Send(comm_buffer, bufcount * sizeof(product_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
	}

#ifndef CLASSIC_FRAGMENTATION
      free(map);
#endif
    }

  /* done */
  return 0;
}


int recv_data(int *mybox, int sender)
{
  /* This routine receives the content of a box from a sender task.
     Communication is divided in these stages:
     1) the target receives the start and length of sender box,
     2) the target sends the start and the length of its box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the target sends a map of needed particles
     5) the target loops on particles and receives them in chunks of size BUFLEN
  */
  
#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
  int nstore;
#endif
  unsigned int size, nsent, boff, expected, received;
  int senderbox[6], interbox[48], i, Nint, box;
  MPI_Status status;

  MPI_Recv(senderbox, 6, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
  MPI_Send(mybox    , 6, MPI_INT, sender, 0, MPI_COMM_WORLD);

  /* up to 8 intersections of the two boxes */
  Nint=intersection(senderbox, mybox, interbox);
  for (box=0; box<Nint; box++)
    {
      boff=box*6;
      size = interbox[3+boff] * interbox[4+boff] * interbox[5+boff];

#ifdef DEBUG
      fprintf(DBGFD,"Task %d will receive this box: %d %d %d -- %d %d %d\n",
	     ThisTask, interbox[0+boff],interbox[1+boff],interbox[2+boff],
	     interbox[3+boff],interbox[4+boff],interbox[5+boff]);
#endif

#ifndef CLASSIC_FRAGMENTATION
      /* send the map of needed particles */
      mapl = size/8 + 1;
      map = (unsigned int*)calloc(mapl, sizeof(unsigned int));

      /* constructs the map for the intersection and send it to the sender */
      build_distmap(map, interbox+boff);
      MPI_Send(map, mapl, MPI_INT, sender, 0, MPI_COMM_WORLD);

      /* map is nulled before getting it back */
      memset(map, 0, mapl * sizeof(unsigned int));

      /* the map is updated by the sender and sent back */
      MPI_Recv(map, mapl, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);

      /* record positions of particles that will be stored 
	 and count how many particles will be sent */
      int off=frag_offset;
      int count=frag_offset;

      for (i=0; i<size; i++)
	if (get_distmap_bit(map, i))
	  {
	    ++count;
	    if (off<subbox.Nalloc)
	      frag_pos[off++]=subbox_space_index(i, interbox+boff);
	  }

      expected = count-frag_offset;
#else
      expected=size;
#endif

      /* receive particles */
      received=0;

      if (expected)
	{

	  do
	    {

	      MPI_Recv(&nsent, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
	      MPI_Recv(comm_buffer, nsent * sizeof(product_data), MPI_BYTE, sender, 0, MPI_COMM_WORLD, &status);
#ifdef DEBUG
	      fprintf(DBGFD,"...receiving %d products from Task %d... %d\n",nsent,sender,(int)sizeof(product_data));
	      if (nsent<10)
		{
		  for (int u=0; u<nsent; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
	      else
		{
		  for (int u=0; u<5; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD," ... ");
		  for (int u=nsent-5; u<nsent; u++)
		    fprintf(DBGFD,"  %f  ",comm_buffer[u].Fmax);
		  fprintf(DBGFD,"\n");
		}
#endif

#ifdef CLASSIC_FRAGMENTATION
	      for (i=0; i<nsent; i++)
		memcpy(&frag[subbox_space_index(i+received, interbox+boff)],&comm_buffer[i],sizeof(product_data));
#else
	      if (frag_offset + nsent < subbox.Nalloc)
		nstore = nsent;
	      else
		nstore = (long long)subbox.Nalloc - (long long)frag_offset;

	      for (i=0; i<nstore; i++)
		memcpy(&frag[frag_offset + i],&comm_buffer[i],sizeof(product_data));
	      frag_offset+=nsent;
#endif
	      received+=nsent;

	    }
	  while (received<expected);

	}

#ifndef CLASSIC_FRAGMENTATION
      free(map);
#endif
    }

  /* done */
  return 0;
}


int keep_data(int *fft_box, int *sub_box)
{
  /* this routine transfers products from the fft space to the subbox space */

#ifndef CLASSIC_FRAGMENTATION
  unsigned int *map, mapl;
#endif
  int interbox[48], i, Nint, box, off;
  unsigned int size;

  Nint = intersection(fft_box, sub_box, interbox);
  for (box=0; box<Nint; box++)
    {
      off=box*6;
      size = interbox[3+off] * interbox[4+off] * interbox[5+off];

#ifdef DEBUG
      fprintf(DBGFD,"Task %d will keep this box: %d %d %d -- %d %d %d\n",
	     ThisTask, interbox[0+off],interbox[1+off],interbox[2+off],interbox[3+off],interbox[4+off],interbox[5+off]);
#endif

#ifdef CLASSIC_FRAGMENTATION
      /* copy data */
      for (i=0; i<size; i++)
	memcpy(&frag[subbox_space_index(i,interbox+off)],
	       &products[fft_space_index(i, interbox+off)],sizeof(product_data));
#else

      /* map of needed particles */
      mapl = size/8 + 1;
      map = (unsigned int*)calloc(mapl, sizeof(unsigned int));
      build_distmap(map, interbox+off);
      update_distmap(map, interbox+off);

      /* copy data */
      for (i=0; i<size; i++)
	{
	  if (get_distmap_bit(map, i))
	    {
	      if (frag_offset<subbox.Nalloc)
		{
		  frag_pos[frag_offset]=subbox_space_index(i,interbox+off);
		  memcpy(&frag[frag_offset], 
			 &products[fft_space_index(i, interbox+off)], sizeof(product_data));
		}
	      ++frag_offset;
	    }
	}
      free(map);
#endif
    }

  return 0;
}


unsigned int fft_space_index(unsigned int pos, int *box)
{
  /* here we move:
     (1) from index i to the relative position of the point in intersection,
     (2) from that to global position without PBCs,
     (3) then we impose PBCs
     (4) then we compute the position in the local FFT box
     (5) and finally we compute the local particle index 
     (these are the logical steps, formulas are more compact)
  */

  unsigned int ip, jp, kp;

  INDEX_TO_COORD(pos, ip, jp, kp, (box+3));

  return 
    COORD_TO_INDEX((ip + box[0]) % MyGrids[0].GSglobal[_x_] - MyGrids[0].GSstart[_x_],
		   (jp + box[1]) % MyGrids[0].GSglobal[_y_] - MyGrids[0].GSstart[_y_],
		   (kp + box[2]) % MyGrids[0].GSglobal[_z_] - MyGrids[0].GSstart[_z_],
		   MyGrids[0].GSlocal);

}


unsigned int subbox_space_index(unsigned int pos, int *box)
{
  /* here we go 
     (1) from index pos to relative position in intersection,
     (2) from that to global position, imposing PBCs,
     (3) then we compute the position in the local subbox, imposing PBCs again,
     (4) and finally the index */

  unsigned int ip, jp, kp;

  INDEX_TO_COORD(pos, ip, jp, kp, (box+3)); /* coords within the intersection */

  return 
    COORD_TO_INDEX(((ip + box[_x_])%MyGrids[0].GSglobal[_x_] - subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_],
		   ((jp + box[_y_])%MyGrids[0].GSglobal[_y_] - subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_],
		   ((kp + box[_z_])%MyGrids[0].GSglobal[_z_] - subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_],
		   subbox.Lgwbl);

}


#ifndef CLASSIC_FRAGMENTATION
int get_distmap_bit(unsigned int *map, unsigned int pos)
{
  /* this operates on the map allocated for the distribution */
  /* gets the bit corresponding to position pos */
  unsigned int rem = pos%UINTLEN;
  return (map[pos/UINTLEN] & (1<<rem))>>rem;
}


void set_distmap_bit(unsigned int *map, unsigned int pos, int value)
{
  /* this operates on the map allocated for the distribution */
  /* sets to 1 the bit corresponding to position pos */
  if (value)
    map[pos/UINTLEN] |= (1<<pos%UINTLEN);
  else
    map[pos/UINTLEN] &= ~(1<<pos%UINTLEN);
  
}


void build_distmap(unsigned int *map, int *box)
{
  /* the receiver builds the map used for distribution */
  unsigned int size = box[3] * box[4] * box[5];

  for (unsigned int i=0; i<size; i++)
    set_distmap_bit(map, i, get_mapup_bit(subbox_space_index(i,box)));

}


void update_distmap(unsigned int *map, int *box)
{

  unsigned int size = box[3] * box[4] * box[5];
  unsigned int bit, i;

  for (i=0; i< size; i++)
    {

      bit=get_distmap_bit(map, i);
      set_distmap_bit(map, i, (bit & (products[fft_space_index(i, box)].Fmax>=outputs.Flast) ) );

    }
}

#endif


#ifdef SNAPSHOT
/* distribution of zacc back to the FFT space */

typedef struct
{
  unsigned int pos;
  PRODFLOAT zacc;
  int group_ID;
} back_data;
back_data *back_buffer;

int keep_data_back(int *);
int send_data_back(int);
int recv_data_back(int *, int);

int distribute_back(void)
{
  /* Distributes accretion times from sub-volumes to fft-space */

  int my_fft_box[6];
  int log_ntask, bit, receiver, sender;

  /* this defines the box that the task possesses in the FFT space */
  my_fft_box[0] = MyGrids[0].GSstart[_x_];
  my_fft_box[1] = MyGrids[0].GSstart[_y_];
  my_fft_box[2] = MyGrids[0].GSstart[_z_];
  my_fft_box[3] = MyGrids[0].GSlocal[_x_];
  my_fft_box[4] = MyGrids[0].GSlocal[_y_];
  my_fft_box[5] = MyGrids[0].GSlocal[_z_];

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  /* stores the data relative to the intersection of the fft box and subbox */
  if (keep_data_back(my_fft_box))
    return 1;

  back_buffer=(back_data*)calloc(BUFLEN , sizeof(back_data));
  if (back_buffer==0x0)
    {
      printf("ERROR on task %d: could not allocate back_buffer in distribute_back\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  /* hypercubic communication scheme */
  for (log_ntask=0; log_ntask<1000; log_ntask++)
    if (1<<log_ntask >= NTasks)
      break;

  /* loop on hypercube dimension */
  for (bit=1; bit<1<<log_ntask; bit++)
    {
    /* loop on tasks */
      for (sender=0; sender<NTasks; sender++)
	{
	  /* receiver task is computed with a bitwise xor */
	  receiver = sender ^ bit;

	  /* condition on sender and receiver */
	  if (receiver < NTasks && sender < receiver)
	    {

	      /* the communication will be done 
		 first sender -> receiver and then receiver -> sender */

	      if (ThisTask==sender)
		send_data_back(receiver);
	      else if (ThisTask==receiver)
		{
		  if (recv_data_back(my_fft_box, sender))
		    return 1;
		}

	      if (ThisTask==receiver)
		send_data_back(sender);
	      else if (ThisTask==sender)
		{
		  if (recv_data_back(my_fft_box, receiver))
		    return 1;
		}

	    }
	}
    }


  free(back_buffer);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}

int keep_data_back(int *fft_box)
{
  /* this routine transfers products from the fft space to the subbox space */

  unsigned int ibox, jbox, kbox, good_particle, fftpos;

  /* find particles that are wanted by the receiver */
  for (int iz=0; iz<subbox.Nstored; iz++)
    {
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);
#endif

      /* this is still in the subbox frame */
      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      /* global box frame */
      ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_];
      jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_];
      kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_];

      if (good_particle &&
	    ibox >= fft_box[0] && ibox < fft_box[0]+fft_box[3] && 
	    jbox >= fft_box[1] && jbox < fft_box[1]+fft_box[4] &&
	    kbox >= fft_box[2] && kbox < fft_box[2]+fft_box[5])
	  {
	    /* position in the fft domain */
	    fftpos = COORD_TO_INDEX(ibox - fft_box[0], jbox - fft_box[1], kbox - fft_box[2], (fft_box+3));
	    products[fftpos].zacc = frag[iz].zacc;
      products[fftpos].group_ID = frag[iz].group_ID;
      
      // printf(" Task %d keep: particle_ID:  %d  zacc: %f   group_ID: %d\n",ThisTask,iz,frag[iz].zacc, frag[iz].group_ID);

	  }
      }

  return 0;
}

int send_data_back(int target)
{
  /* This routine sends the content of a box to a target task.
     Communication is divided in these stages:
     1) the sender communicates the start and length of its box,
     2) the sender receives the start and the length of the needed box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the sender receives a map of needed particles
     5) the sender loops on particles and sends them in chunks of size BUFLEN
  */
  
  int ibox, jbox, kbox, good_particle;
  int recv_box[6];
  MPI_Status status;

  MPI_Recv(recv_box, 6, MPI_INT, target, 0, MPI_COMM_WORLD, &status);

//  printf("BOX Task %d: [%d, %d, %d] - [%d, %d, %d]\n",ThisTask,recv_box[0],recv_box[1],recv_box[2],recv_box[3],recv_box[4],recv_box[5]); //LEVARE

  int bufcount=0;
  int sent=0;

  /* find particles that are wanted by the receiver */
  for (int iz=0; iz<subbox.Nstored; iz++)
    {
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);
#endif
      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      /* global box frame */
      ibox = (ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_];
      jbox = (jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_];
      kbox = (kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_];

      //printf("  oioiuoi %d %d %d %d %d\n" , iz, ibox, jbox, kbox, good_particle); // LEVARE

      if (good_particle &&
	    ibox >= recv_box[0] && ibox < recv_box[0]+recv_box[3] && 
	    jbox >= recv_box[1] && jbox < recv_box[1]+recv_box[4] &&
	    kbox >= recv_box[2] && kbox < recv_box[2]+recv_box[5])
	  {
	    /* position in the fft domain */
	    back_buffer[bufcount].pos = COORD_TO_INDEX(ibox - recv_box[0], jbox - recv_box[1], kbox - recv_box[2], (recv_box+3));
	    back_buffer[bufcount].zacc = frag[iz].zacc;
      back_buffer[bufcount].group_ID =frag[iz].group_ID;

	    // printf(" Task %d send:  %d  %d  %f %d\n",ThisTask,iz+1,COORD_TO_INDEX(ibox,jbox,kbox,MyGrids[0].GSglobal)+1,back_buffer[iz].zacc,back_buffer[iz].group_ID); // LEVARE
	    ++bufcount;
	    ++sent;

	    if (bufcount==BUFLEN)
	      {
		MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
		MPI_Send(back_buffer, bufcount * sizeof(back_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
		bufcount=0;
	      }		
	  }
    }

  if (bufcount)
    {
      MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
      MPI_Send(back_buffer, bufcount * sizeof(back_data), MPI_BYTE, target, 0, MPI_COMM_WORLD);
    }
  
  bufcount=0;
  MPI_Send(&bufcount, 1, MPI_INT, target, 0, MPI_COMM_WORLD);

//  printf("Task %d sent %d particles\n",ThisTask,sent); // LEVARE

  /* done */
  return 0;
}

int recv_data_back(int *mybox, int sender)
{
  /* This routine receives the content of a box from a sender task.
     Communication is divided in these stages:
     1) the target receives the start and length of sender box,
     2) the target sends the start and the length of its box,
     3) intersection of the sender and target boxes is computed
     4) if there is an intersection the target sends a map of needed particles
     5) the target loops on particles and receives them in chunks of size BUFLEN
  */
  
  unsigned int nsent, received;
  MPI_Status status;

  MPI_Send(mybox, 6, MPI_INT, sender, 0, MPI_COMM_WORLD);

  /* receive particles */
  received=0;

  do
    {
      MPI_Recv(&nsent, 1, MPI_INT, sender, 0, MPI_COMM_WORLD, &status);
      if (nsent)
	MPI_Recv(back_buffer, nsent * sizeof(back_data), MPI_BYTE, sender, 0, MPI_COMM_WORLD, &status);
      for (int iz=0; iz<nsent; iz++)
	{
	  products[back_buffer[iz].pos].zacc=back_buffer[iz].zacc;
    products[back_buffer[iz].pos].group_ID=back_buffer[iz].group_ID;

    // printf(" Task %d recv:  %d  %d  %f  %d\n",ThisTask,iz+1,back_buffer[iz].pos+1,back_buffer[iz].zacc, back_buffer[iz].group_ID);
	}
      received+=nsent;
    }
  while (nsent);

//  printf("Task %d received %d particles\n",ThisTask,received); // LEVARE

  /* done */
  return 0;
}

#endif
