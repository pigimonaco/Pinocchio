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

//#define VERBOSE

typedef struct
{
  int task,box_x,box_y,box_z,bz,gz,flag,check;  // flag si puo` levare...
} comm_struct;

static int largest_size;
static product_data *sub_plane;

int find_task(int, int, int, int, int *, int *, int *);
int send_data(comm_struct *);
int recv_data(comm_struct *);
int keep_data(comm_struct *);

int distribute(void)
{
  /*
  Distributes products from planes to sub-volumes, including boundary layers
  */

  int global_z, *belongs_local, *belongs_global;
  int i1,i2,i3,main_task,second_task,sender,box_z,second_z;

  int log_ntask,bit,receiver,ss,rr,dd,turn,i,count_send,count_recv,count_keep,isend,irecv;
  comm_struct *comm_send, *comm_recv, *comm_keep;

  /* let's synchronize the tasks here */
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

#ifdef VERBOSE
  printf("\n");
  printf("task %d, plan of sub-boxes:  %d %d %d - %d %d %d/%d - %d %d %d - %d %d %d - %d %d %d - %d %d %d - %d %d %d\n",
	 ThisTask,MyGrids[0].GSlocal_x,MyGrids[0].GSlocal_y,MyGrids[0].GSlocal_z,
	 subbox.nbox_x, subbox.nbox_y, subbox.nbox_z_thisslice, subbox.nbox_z_allslices,
	 subbox.mybox_x,subbox.mybox_y,subbox.mybox_z,
	 subbox.Lgrid_x,subbox.Lgrid_y,subbox.Lgrid_z,
	 subbox.start_x,subbox.start_y,subbox.start_z,
	 subbox.Lgwbl_x,subbox.Lgwbl_y,subbox.Lgwbl_z,
	 subbox.stabl_x,subbox.stabl_y,subbox.stabl_z);
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /* here it builds a map of which task owns which plane */
  /* NB: here it is assumed that tasks are distributed by planes */

  belongs_local  = (int*)malloc(MyGrids[0].GSglobal_z*sizeof(int));
  belongs_global = (int*)malloc(MyGrids[0].GSglobal_z*sizeof(int));
  for (global_z=0; global_z<MyGrids[0].GSglobal_z; global_z++)
    {
      belongs_global[global_z]=0;
      if (global_z >= MyGrids[0].GSstart_z && global_z < MyGrids[0].GSstart_z + MyGrids[0].GSlocal_z)
	belongs_local[global_z]=ThisTask;
      else
	belongs_local[global_z]=0;
    }  

  if (MPI_Reduce(belongs_local, belongs_global, 
		 MyGrids[0].GSglobal_z, MPI_INT, MPI_SUM, 0, 
		 MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: distribute could not perform an MPI_Reduce\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  if (MPI_Bcast(belongs_global, MyGrids[0].GSglobal_z, 
		MPI_INT, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: distribute could not perform an MPI_Bcast\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  /* loop on global planes to count the number of communications needed */
  /* i3 is the counter of sub-boxes, on the z direction, updated at the end of the do cycle */

  for (global_z = i3=count_send=count_recv=count_keep=0; global_z<MyGrids[0].GSglobal_z; global_z++)
    {
      /* updates the z-sub-box index if necessary
         NB: i3 is related to the complete box, with all slices */
      if (i3<subbox.nbox_z_allslices-1 && global_z==find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,i3+1))
        i3++;

       /* loop on all sub-planes, that is sub-boxes on the plane */
      for (i2=0; i2<subbox.nbox_y; i2++)
	for (i1=0; i1<subbox.nbox_x; i1++)
	  {

	    /* find the main and second task that will need the sub-plane */
	    main_task = find_task(global_z, i1, i2, i3, &second_task, &box_z, &second_z);

	    /* check if the task must send, receive or keep data */
	    if (belongs_global[global_z]==ThisTask)
	      {
		if (main_task==ThisTask)
		  count_keep++;
		else if (main_task>=0)
		  count_send++;

		if (second_task==ThisTask)
		  count_keep++;
		else if (second_task>=0)
		  count_send++;
	      }

	    if (main_task==ThisTask && belongs_global[global_z]!=ThisTask)
	      count_recv++;

	    if (second_task==ThisTask && belongs_global[global_z]!=ThisTask)
              count_recv++;
	  }
    }
  /* allocate buffers for sendind, receiving, and keeping */

  comm_send=(comm_struct*)malloc(count_send * sizeof(comm_struct));
  comm_recv=(comm_struct*)malloc(count_recv * sizeof(comm_struct));
  comm_keep=(comm_struct*)malloc(count_keep * sizeof(comm_struct));

  /* loop on global planes to record the needed communications */
  for (global_z = i3=count_send=count_recv=count_keep=0; global_z<MyGrids[0].GSglobal_z; global_z++)
    {
      /* updates the z-sub-box index if necessary */
      if (i3<subbox.nbox_z_allslices-1 && global_z==find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,i3+1))
        i3++;

       /* loop on all sub-planes, that is sub-boxes on the plane */
      for (i2=0; i2<subbox.nbox_y; i2++)
	for (i1=0; i1<subbox.nbox_x; i1++)
	  {

	    /* find the main and second task that will need the sub-plane */
	    main_task = find_task(global_z, i1, i2, i3, &second_task, &box_z, &second_z);

	    /* store the information */
	    if (belongs_global[global_z]==ThisTask)
	      {
		if (main_task==ThisTask)
		  {
		    comm_keep[count_keep  ].task  = ThisTask;   /* sender or receiver */
		    comm_keep[count_keep  ].box_x = i1;         /* x-index of subbox to be processed */ 
		    comm_keep[count_keep  ].box_y = i2;         /* y-index of subbox to be processed */ 
		    comm_keep[count_keep  ].box_z = i3;         /* z-index (global) of subbox to be processed */ 
		    comm_keep[count_keep  ].gz    = global_z;   /* global plane where to find data */
		    comm_keep[count_keep  ].bz    = box_z;      /* local plane where to store data */
		    comm_keep[count_keep  ].flag  = 1;
		    comm_keep[count_keep++].check = 0;
		  }
		else if (main_task>=0)
		  {
		    comm_send[count_send  ].task  = main_task;
		    comm_send[count_send  ].box_x = i1;
		    comm_send[count_send  ].box_y = i2;
		    comm_send[count_send  ].box_z = i3;
		    comm_send[count_send  ].gz    = global_z;		    
		    comm_send[count_send  ].bz    = box_z;
		    comm_send[count_send  ].flag  = 1;
		    comm_send[count_send++].check = 0;
		  }

		if (second_task==ThisTask)
		  {
		    comm_keep[count_keep  ].task  = ThisTask;
		    comm_keep[count_keep  ].box_x = i1;
		    comm_keep[count_keep  ].box_y = i2;
		    comm_keep[count_keep  ].box_z = i3;
		    comm_keep[count_keep  ].gz    = global_z;
		    comm_keep[count_keep  ].bz    = second_z;
		    comm_keep[count_keep  ].flag  = 2;
		    comm_keep[count_keep++].check = 0;
		  }
		else if (second_task>=0)
		  {
		    comm_send[count_send  ].task  = second_task;
		    comm_send[count_send  ].box_x = i1;
		    comm_send[count_send  ].box_y = i2;
		    comm_send[count_send  ].box_z = i3;
		    comm_send[count_send  ].gz    = global_z;
		    comm_send[count_send  ].bz    = second_z;
		    comm_send[count_send  ].flag  = 2;
		    comm_send[count_send++].check = 0;
		  }
	      }

	    if (main_task==ThisTask && belongs_global[global_z]!=ThisTask)
	      {
		comm_recv[count_recv  ].task  = belongs_global[global_z];
		comm_recv[count_recv  ].box_x = i1;
		comm_recv[count_recv  ].box_y = i2;
		comm_recv[count_recv  ].box_z = i3;
		comm_recv[count_recv  ].gz    = global_z;
		comm_recv[count_recv  ].bz    = box_z;
		comm_recv[count_recv  ].flag  = 1;
		comm_recv[count_recv++].check = 0;
	      }

	    if (second_task==ThisTask && belongs_global[global_z]!=ThisTask)
	      {
		comm_recv[count_recv  ].task  = belongs_global[global_z];
		comm_recv[count_recv  ].box_x = i1;
		comm_recv[count_recv  ].box_y = i2;
		comm_recv[count_recv  ].box_z = i3;
		comm_recv[count_recv  ].gz    = global_z;
		comm_recv[count_recv  ].bz    = second_z;
		comm_recv[count_recv  ].flag  = 2;
		comm_recv[count_recv++].check = 0;
	      }
	  }
    }

#ifdef VERBOSE
  int itask;
  for (itask=0; itask<NTasks; itask++)
    {
      if (ThisTask==itask)
	{
	  printf("Task %d\n",ThisTask);
	  printf("keep_data: %d\n",count_keep);
	  for (i=0; i<count_keep; i++)
	    printf("task=%d, box_x=%d, box_y=%d, box_z=%d, global_z=%d, local_z=%d, flag=%d\n",
		   comm_keep[i].task,
		   comm_keep[i].box_x,
		   comm_keep[i].box_y,
		   comm_keep[i].box_z,
		   comm_keep[i].gz,
		   comm_keep[i].bz,
		   comm_keep[i].flag);
	  printf("send_data: %d\n",count_send);
	  for (i=0; i<count_send; i++)
	    printf("task=%d, box_x=%d, box_y=%d, box_z=%d, global_z=%d, local_z=%d, flag=%d\n",
		   comm_send[i].task,
		   comm_send[i].box_x,
		   comm_send[i].box_y,
		   comm_send[i].box_z,
		   comm_send[i].gz,
		   comm_send[i].bz,
		   comm_send[i].flag);
	  printf("recv_data: %d\n",count_recv);
	  for (i=0; i<count_recv; i++)
	    printf("task=%d, box_x=%d, box_y=%d, box_z=%d, global_z=%d, local_z=%d, flag=%d\n",
		   comm_recv[i].task,
		   comm_recv[i].box_x,
		   comm_recv[i].box_y,
		   comm_recv[i].box_z,
		   comm_recv[i].gz,
		   comm_recv[i].bz,
		   comm_recv[i].flag);
	}
    }
#endif

  largest_size=(find_length(MyGrids[0].GSglobal_x,subbox.nbox_x,0) + 2.*subbox.safe_x) *
    (find_length(MyGrids[0].GSglobal_y,subbox.nbox_y,0) + 2.*subbox.safe_y);
  sub_plane=(product_data*)malloc(largest_size * sizeof(product_data));
  if (sub_plane==0x0)
    {
      printf("ERROR on task %d: could not allocate sub_plane in distribute\n",ThisTask);
      fflush (stdout);
      return 1;
    }

  /* assign the bunches of memory that should not be communicated */
  for (i=0; i<count_keep; i++)
    {
      if (keep_data(comm_keep+i))
	return 1;
      comm_keep[i].check=1;
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

	      /* the communication will be done sender -> receiver and receiver -> sender */
	      ss=sender;
	      rr=receiver;
#ifdef VERBOSE
	      if (!ThisTask)
		printf("COMMUNICATION %d %d\n",ss,rr);
#endif
	      for (turn=0; turn<2; turn++)
		{
		  if (ThisTask==ss)
		    {
		      for (isend=0; isend<count_send; isend++)
			if (comm_send[isend].task==rr)
			  {
			    if (comm_send[isend].check) 
			      printf("CASINO SEND %d %d %d %d\n",ThisTask,rr,ss,isend);
			    if (send_data(comm_send+isend))
			      return 1;
			    comm_send[isend].check=1;
			  }
		    }
		  else if (ThisTask==rr)
		    {
		      for (irecv=0; irecv<count_recv ; irecv++)
			if (comm_recv[irecv].task==ss)
			  {
			    if (comm_recv[irecv].check) 
			      printf("CASINO RECV %d %d %d %d\n",ThisTask,rr,ss,isend);
			    if (recv_data(comm_recv+irecv))
			      return 1;
			    comm_recv[irecv].check=1;
			  }
		    }		

		  /* exchanges sender and receiver */
		  dd=rr;
		  rr=ss;
		  ss=dd;
		}

	    }
	}
    }


  /* check that all wanted communications have been performed */
  for (i=0; i<count_keep; i++)
    if (!comm_keep[i].check) 
      printf("CASINO FINALE KEEP %d %d %d\n",ThisTask,i,comm_keep[i].task);
  for (i=0; i<count_send; i++)
    if (!comm_send[i].check) 
      printf("CASINO FINALE SEND %d %d %d\n",ThisTask,i,comm_send[i].task);
  for (i=0; i<count_recv; i++)
    if (!comm_recv[i].check) 
      printf("CASINO FINALE RECV %d %d %d\n",ThisTask,i,comm_recv[i].task);

  free(sub_plane);

  free(comm_keep);
  free(comm_recv);
  free(comm_send);

  free(belongs_global);
  free(belongs_local);

  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


int find_task(int global_z, int i1, int i2, int i3, int *second_task, int *box_z, int *second_z)
{
  /* given the z-plane, find the main and secondary task
     to which the sub-plane must be communicated */

  int main_task,i3bis,TL;

  /* sub-planes are distributed to tasks first in z, then in y and x */
  if (i3>=subbox.nbox_z_thisslice*ThisSlice && i3<subbox.nbox_z_thisslice*(ThisSlice+1))
    main_task=subbox.nbox_z_thisslice*(i1*subbox.nbox_y+i2)+(i3-subbox.nbox_z_thisslice*ThisSlice);
  else
    main_task=-1;

  /* this is the local z variable in the subbox-system WITHOUT the boundary layer */
  *box_z=global_z-find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,i3);
  TL=find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,i3);

  /* the second task to communicate the sub_plane to is decided using box_z */
  if (*box_z<subbox.safe_z)            // lower border, second task is the previous in z
    {
      i3bis=i3-1;
      if (i3bis<0)
	i3bis+=subbox.nbox_z_allslices;
      if (i3bis>=subbox.nbox_z_thisslice*ThisSlice && i3bis<subbox.nbox_z_thisslice*(ThisSlice+1))
	{
	  *second_task = subbox.nbox_z_thisslice*(i1*subbox.nbox_y+i2)+(i3bis-subbox.nbox_z_thisslice*ThisSlice);
	  *second_z = find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,i3bis) + *box_z;
	}
      else
	{
	  *second_task=-1;
	  *second_z=-99;
	}
    }
  else if (*box_z>=TL-subbox.safe_z)   // upper border, second task is the next in z
    {
      i3bis=i3+1;
      if (i3bis>=subbox.nbox_z_allslices)
	i3bis-=subbox.nbox_z_allslices;
      if (i3bis>=subbox.nbox_z_thisslice*ThisSlice && i3bis<subbox.nbox_z_thisslice*(ThisSlice+1))
	{
	  *second_task = subbox.nbox_z_thisslice*(i1*subbox.nbox_y+i2)+(i3bis-subbox.nbox_z_thisslice*ThisSlice);
	  *second_z = *box_z - TL;
	}
      else
	{
	  *second_task=-1;
	  *second_z=-99;
	}
    }
  else
    {
      *second_task=-1;
      *second_z=-99;
    }

  /* z-coordinates are fixed to take into account the bounday layer */
  *box_z += subbox.safe_z;
  *second_z += subbox.safe_z;

#ifdef VERBOSE
  //  printf("find_task: Task %d, %d    %d %d %d    %d %d   %d %d\n", 
  //	 ThisTask, global_z, i1, i2, i3, main_task, *second_task, *box_z, *second_z);
  //fflush(stdout);
#endif

  return main_task;
}


int send_data(comm_struct *comm_data)
{
  /* communication -> sending */

  int j2,k2,j1,k1,k2start,k1start;
  int Lgrid_x, Lgrid_y, offset;
  size_t size;

#ifdef VERBOSE
  printf("SEND_DATA: Task %d sending to %d\n",ThisTask,comm_data->task);
#endif

  Lgrid_x = find_length(MyGrids[0].GSglobal_x,subbox.nbox_x,comm_data->box_x) + 2.*subbox.safe_x;
  Lgrid_y = find_length(MyGrids[0].GSglobal_y,subbox.nbox_y,comm_data->box_y) + 2.*subbox.safe_y;

  /* allocates the plane for communications */
  size = Lgrid_x * Lgrid_y;
  if (size>largest_size)   /* questo poi si toglie */
    {
      printf("ASSURDO: size>largest_size in SEND_DATA!!!\n");
      return 1;
    }

  /* starting coordinates of the sub-box (can be negative) */
  k2start=find_start(MyGrids[0].GSglobal_y,subbox.nbox_y,comm_data->box_y)-subbox.safe_y;
  k1start=find_start(MyGrids[0].GSglobal_x,subbox.nbox_x,comm_data->box_x)-subbox.safe_x;

  /* storing of information on sub_plane */
  /* NB: anche qui assume che la distribuzione della memoria sia in piani */
  offset = (comm_data->gz - MyGrids[0].GSstart_z) * MyGrids[0].GSlocal_x *  MyGrids[0].GSlocal_y;
  for (j2=0; j2<Lgrid_y; j2++)
    {
      k2=j2+k2start;
      if (k2<0) k2+=MyGrids[0].GSglobal_y;
      if (k2>=MyGrids[0].GSglobal_y) k2-=MyGrids[0].GSglobal_y;
      for (j1=0; j1<Lgrid_x; j1++)
	{
	  k1=j1+k1start;
	  if (k1<0) k1+=MyGrids[0].GSglobal_x;
	  if (k1>=MyGrids[0].GSglobal_x) k1-=MyGrids[0].GSglobal_x;

	  *(sub_plane + j1 + j2*Lgrid_x) = 
	    *(products + k1 + MyGrids[0].GSlocal_x * k2 + offset);
	}
    }

  /* sending sub_plane to the receiver */
  if (MPI_Send((void*)sub_plane, size * sizeof(product_data), 
	       MPI_BYTE, comm_data->task, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: failed communication in send_data\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  /* done */
  return 0;
}

int recv_data(comm_struct *comm_data)
{

  /* communication -> receiving */

  int i,offset;
  size_t size;
  MPI_Status status;

#ifdef VERBOSE
  printf("RECV_DATA: Task %d receiving from %d\n",ThisTask,comm_data->task);
#endif

  /* allocates the plane for communications */
  size = subbox.Lgwbl_x * subbox.Lgwbl_y;
  if (size>largest_size)   /* questo poi si toglie */
    {
      printf("ASSURDO: size>largest_size in RECV_DATA!!!\n");
      return 1;
    }

  /* receiving information from the sender */
  if (MPI_Recv((void*)sub_plane, size * sizeof(product_data), 
	       MPI_BYTE, comm_data->task, 0, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
    {
      printf("ERROR on task %d: failed communication in recv_data\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  /* copying information onto the frag structure */
  offset = comm_data->bz * size;
  for (i=0; i<size; i++)
    *(frag + i + offset) = *(sub_plane+i);
  
  /* done */
  return 0;
}

int keep_data(comm_struct *comm_data)
{
  int j2,k2,j1,k1,offmax,offrag;

  /* storing of information on the product structure */
  /* NB: anche qui assume che la distribuzione della memoria sia in piani */
  offmax = (comm_data->gz - MyGrids[0].GSstart_z) * MyGrids[0].GSlocal_x *  MyGrids[0].GSlocal_y;
  offrag = comm_data->bz * subbox.Lgwbl_x * subbox.Lgwbl_y;

  for (j2=0; j2<subbox.Lgwbl_y; j2++)
    {
      k2=j2+subbox.stabl_y;
      if (k2<0) k2+=MyGrids[0].GSlocal_y;
      if (k2>=MyGrids[0].GSlocal_y) k2-=MyGrids[0].GSlocal_y;
      for (j1=0; j1<subbox.Lgwbl_x; j1++)
	{
	  k1=j1+subbox.stabl_x;
	  if (k1<0) k1+=MyGrids[0].GSlocal_x;
	  if (k1>=MyGrids[0].GSlocal_x) k1-=MyGrids[0].GSlocal_x;

	  *(frag + j1 + j2*subbox.Lgwbl_x + offrag) =
	    *(products + k1 + MyGrids[0].GSlocal_x * k2 + offmax);
	}
    }

  return 0;
}
