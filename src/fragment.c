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
#include <sys/types.h>
#include <sys/stat.h>

#define R_THR 2
#define A_THR 9

//#define DEBUG
#define FTOZ(A) (A>0.0 ? A-1 : A)

int fragment(void);
int fragment_simpler(void);
int count_peaks(volume_data *);
int find_all_peaks();
void set_fragment_parameters(int);
void write_map(int);
int create_map(void);
int initialize_volume(int, int, int, int, int, int, volume_data *);
int merge_catalogs(volume_data *, int);

int check_consistency();

#ifdef RECOMPUTE_DISPLACEMENTS
int shift_all_displacements(void);
int recompute_group_velocities(void);
#endif

char *volume_space;
size_t volume_length;

static int index_compare_F(const void * a,const void * b)
{
  if ( groups[*((int *)a)].t_peak == groups[*((int *)b)].t_peak )
    return 0;
  else
    if ( groups[*((int *)a)].t_peak > groups[*((int *)b)].t_peak )
      return -1;
    else
      return 1;
}


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


int fragment_driver()
{
  /* This routine can be edited to call fragmentation for many combinations 
     of fragment parameters */

  double tmp;

  if (!ThisTask)
    printf("\n[%s] Second part: fragmentation of the collapsed medium\n",fdate());

  /* assigns fragmentation parameters */
  set_fragment_parameters(ORDER_FOR_GROUPS);
  
  /* reallocates the memory */
  if (reallocate_memory_for_fragmentation())
    return 1;

  /* timing */
  cputime.frag=MPI_Wtime();
  cputime.group=cputime.sort=cputime.distr=0.0;
#ifdef PLC
  cputime.plc=0.0;
#endif

  /* redistribution of products */
  if (!ThisTask)
    printf("[%s] Starting re-distribution of products\n",fdate());

  tmp=MPI_Wtime();

  map_to_be_used=0;
  if (distribute())
    return 1;

  tmp=MPI_Wtime()-tmp;
  cputime.distr += tmp;

  if (!ThisTask)
    printf("[%s] Re-distribution of products done, cputime = %14.6f\n", fdate(), tmp);

  /* sets arrays to zero */
  memset(group_ID,0,subbox.Nalloc*sizeof(int));
  memset(linking_list,0,subbox.Nalloc*sizeof(int)); // questa non serve piu`, ci metto i picchi ma poi la leveremo
  memset(groups,0,subbox.PredNpeaks*sizeof(group_data));
  memset(wheretoplace_mycat,0,subbox.PredNpeaks*sizeof(histories_data));
#ifdef PLC
  memset(plcgroups,0,plc.Nmax*sizeof(plcgroup_data));
#endif

  /* memory for the small volume */
  // QUESTA LA DOBBIAMO PIAZZARE DA QUALCHE PARTE NELLA MEMORIA GLOBALE
  int sizemax=MyGrids[0].GSglobal[0];
  volume_length = sizemax * sizemax * sizemax * (sizeof(product_data) + 2*sizeof(int) + sizeof(group_data));
  volume_space = (char*)malloc(volume_length);

  printf("riservato al volume: %lu\n",volume_length);  // LEVARE

  /* finds all the peaks and initializes the groups */
  int Npeaks = find_all_peaks();
  printf("Found a total of %d peaks\n",Npeaks);

  /* starts to fragment small volumes */
  volume_data my_volume;
  int Nsub=4;
  int size=MyGrids[0].GSglobal[0]/Nsub;  // questa deve essere scalata alla dimensione della subbox

  printf("size: %d, min N: %d\n",size,(int)pow((double)size/2.74,3.));

  check_consistency();  // verifica che i merger tree partano bene
  
  /* primo loop su una griglia */
  printf("primo loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  // le coordinate sono nel frame della subbox
	  if (initialize_volume(isub*size-1,jsub*size-1,ksub*size-1,
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("primo loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}

  if (check_consistency())
    return 1;


  /* secondo loop spostandosi di mezza griglia in x */
  printf("secondo loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1+(int)(size/2),jsub*size-1,ksub*size-1,
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("secondo loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}

  if (check_consistency())
    return 1;


  /* terzo loop spostandosi di mezza griglia in y */
  printf("terzo loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1,jsub*size-1+(int)(size/2),ksub*size-1,
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 
	  
	  //printf("terzo loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}

  if (check_consistency())
    return 1;

  /* quarto loop spostandosi di mezza griglia in z */
  printf("quarto loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1,jsub*size-1,ksub*size-1+(int)(size/2),
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("quarto loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}


  if (check_consistency())
    return 1;

  /* quinto loop spostandosi in diagonale xy */
  printf("quinto loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1+(int)(size/2),jsub*size-1+(int)(size/2),ksub*size-1,
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("quarto loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}


  if (check_consistency())
    return 1;



  /* sesto loop spostandosi in diagonale xz*/
  printf("sesto loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1+(int)(size/2),jsub*size-1,ksub*size-1+(int)(size/2),
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("quarto loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}


  if (check_consistency())
    return 1;

  /* settimo loop spostandosi in diagonale yz */
  printf("settimo loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1,jsub*size-1+(int)(size/2),ksub*size-1+(int)(size/2),
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("quarto loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}


  if (check_consistency())
    return 1;

  /* ottavo loop spostandosi in diagonale xyz */
  printf("ottavo loop\n");
  for (int isub=0; isub<Nsub; isub++)
    for (int jsub=0; jsub<Nsub; jsub++)
      for (int ksub=0; ksub<Nsub; ksub++)
	{
	  
	  if (initialize_volume(isub*size-1+(int)(size/2),jsub*size-1+(int)(size/2),ksub*size-1+(int)(size/2),
				size+2,size+2,size+2, &my_volume))
	    return 1;

	  if (build_groups_in_volume(&my_volume))
	    return 1; 

	  //printf("quarto loop: volume %d %d %d, ngroups=%d\n",isub,jsub,ksub,ngroups);

	  if (merge_catalogs(&my_volume,1))
	    return 1;

	}


  if (check_consistency())
    return 1;

  #define FACTOR 5

  /* ultimo loop centrato sugli aloni non ancora risolti */
  int ib,jb,kb;
  for (int g=FILAMENT+1; g<ngroups; g++)
    {
      if (!groups[g].good && groups[g].Mass>1)
	{
	  if (groups[g].halo_app!=g)
	    printf("Allora non ho capito un cazzo!\n");

	  /* the volume size is requested to be odd */
	  size=MyGrids[0].GSglobal[0]/Nsub+1;
	  if (!(size%2))
	    size++;
	  /* the center of the box is the peak */
	  INDEX_TO_COORD(groups[g].name,ib,jb,kb,MyGrids[0].GSglobal);  // DA TRASFORMARE IN COORDINATE DELLA SUBBOX!
	  if (!groups[g].good)
	    do 
	      {
		printf("Calling group %d of mass %d (%d %d %d) with size %d\n",g,groups[g].Mass,ib,jb,kb,size);
		if (initialize_volume(ib-size/2,jb-size/2,kb-size/2,size,size,size, &my_volume))
		  return 1;

		if (build_groups_in_volume(&my_volume))
		  return 1; 

		if (merge_catalogs(&my_volume,1))
		  return 1;

		if (!groups[g].good)
		  size=(int)(size*1.2); // THIS CAN BE IMPROVED!!!
	      }
	    while (!groups[g].good);
	}

    }

  if (check_consistency())
    return 1;


  /* halo catalog */
  if (write_catalog(outputs.n-1))
    return 1;

  /* halo mass function */
  if (compute_mf(outputs.n-1))
    return 1;

  /* merger histories, only at the last redshift */
  if (write_histories())
    return 1;

  free(volume_space);

  return 0;
}


int check_consistency()
{

  static int call=0;
  int err=0;

  // aux1 lo uso per capire se ho gia` visitato un alone
  for (int i=FILAMENT+1; i<ngroups; i++)
    groups[i].aux1=0;

  // controllo che tutte le merger histories siano consistenti
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (groups[i].halo_app==i)
      {
	int next=i;
	int count=0;
	int loop=0;
	do
	  {
	    if (groups[next].aux1)
	      printf("errore: branch %d del tree %d visto due volte: %d\n",next,i,groups[next].aux1);
	    groups[next].aux1=i;

	    if (count && !loop && groups[next].halo_app==next)
	      {
		printf("errore: da %d siamo atterrati sul tree %d \n",i,next);
		loop++;
		err++;
	      }

	    count++;
	    if (count>1000)
	      {
		printf("errore: alluppato sul tree %d\n",i);
		err++;
		break;
	      }

	    next=groups[next].ll;
	  }
	while (next != i);
      }

  // controllo se sono rimasti fuori aloni
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (!groups[i].aux1)
      {
	printf("errore: l`alone %d e` fuori da ogni tree\n",i);
	err++;
      }

  // conto il numero di aloni non assegnati
  int counter=0;
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (!groups[i].good)
      counter++;
  printf("Aloni non assegnati: %d\n",counter);
  counter=0;
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (groups[i].Mass>1 && !groups[i].good)
      counter++;
  printf("Aloni principali non assegnati: %d\n",counter);

  // controllo che gli aloni assegnati siano good
  counter=0;
  for (int i=FILAMENT+1; i<ngroups; i++)
    if ((!groups[i].good) && groups[i].halo_app!=i)
     counter++;

 if (counter)
   printf("ATTENZIONE, %d branch di qualche albero hanno good=0\n",counter);


  // consistenza dei group_ID
  // aux1 lo uso per contare le particelle in alone
  for (int i=FILAMENT+1; i<ngroups; i++)
    groups[i].aux1=0;

  for (int pos=0; pos<subbox.Npart; pos++)
    {
      if (group_ID[pos]>=ngroups)
	printf("QUESTO NO! %d %d\n",pos,group_ID[pos]);
      
      if (group_ID[pos]>FILAMENT)
	groups[group_ID[pos]].aux1++;
    }
 
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (groups[i].good && groups[i].halo_app==i && groups[i].aux1!=groups[i].Mass)
      printf("ERRORE: non mi torna la massa del gruppo %d: %d da groups, %d da group_ID\n",i,groups[i].Mass,groups[i].aux1);

  // consistenza delle linking list dei group_ID
  int cp=0,cR=0;
  for (int i=FILAMENT+1; i<ngroups; i++)
    {
      if (groups[i].good && groups[i].halo_app==i)
	{
	  int pos=groups[i].point;
	  counter=cR=cp=0;
	  do
	    {
	      if (pos<0)
		break;

	      if (group_ID[pos]!=i)
		cp++;
	      //  printf("ERRORE: la particella %d non risulta appartenere all'alone %d - %d\n",pos,i,group_ID[pos]);
	      if (groups[i].good && frag[pos].Rmax>=0)
		cR++;
	      // printf("ERRORE: la particella %d dell'alone %d ha un Rmax non negativo %d\n",pos,i,frag[pos].Rmax);
	      pos=linking_list[pos]; 
	      counter++;
	      if (counter>1000)
		break;
	    }
	  while (pos!=groups[i].point);

	  if ((counter!=groups[i].Mass) || cp || cR)
	   printf("ERRORE nella mappa delle particelle dell'alone %d: counter=%d, mass=%d, cp=%d, cR=%d\n",i,counter,groups[i].Mass,cp,cR);
	}
    }

  char filename[100];
  sprintf(filename,"check_%d",call);
  printf("writing %s\n",filename);
  FILE *fd=fopen(filename,"w");
  for (int i=FILAMENT+1; i<ngroups; i++)
    if (!groups[i].good)
      {
	int ib,jb,kb;
	INDEX_TO_COORD(groups[i].name,ib,jb,kb,MyGrids[0].GSglobal);
	fprintf(fd," %12d %12Ld %12d %12d %12.8f  %3d %3d %3d\n",i,groups[i].name,groups[i].Mass,groups[i].aux2,
		groups[i].aux2/pow(groups[i].Mass,1./3.),  ib,jb,kb);
      }
  fclose(fd);
  

  fflush(stdout);
  ++call;
  return err;

}



#define SET_PBC(I,F,L) ( (F?(I+L)%L:I) )

int initialize_volume(int ic, int jc, int kc, int sizex, int sizey, int sizez, volume_data *my_volume)
{

  /* writes on my_volume all the needed information */
  // IN QUESTO MOMENTO ASSUME DI AVERE TUTTE LE INFORMAZIONI NEL SUBVOLUME

  size_t count_memory;

  memset(volume_space, 0, volume_length);

  my_volume->Start[_x_]=ic;
  my_volume->Start[_y_]=jc;
  my_volume->Start[_z_]=kc;

  my_volume->GridSize[0]=sizex;
  my_volume->GridSize[1]=sizey;
  my_volume->GridSize[2]=sizez;
  my_volume->Npart=my_volume->GridSize[0]*my_volume->GridSize[1]*my_volume->GridSize[2];

  my_volume->Frag = (product_data *)volume_space;
  count_memory = my_volume->Npart*sizeof(product_data);
  my_volume->Group_ID = (int *)(volume_space + count_memory);
  count_memory += my_volume->Npart*sizeof(int);
  my_volume->Linking_list = (int *)(volume_space + count_memory);
  count_memory += my_volume->Npart*sizeof(int);

  //FILE *file=fopen("volumetto","w");

  for (int iv=0; iv<my_volume->GridSize[_x_]; iv++) 
    {
      int ib=SET_PBC(iv+ic,subbox.pbc[_x_],subbox.Lgwbl[_x_]);
      for (int jv=0; jv<my_volume->GridSize[_y_]; jv++)
	{
	  int jb=SET_PBC(jv+jc,subbox.pbc[_y_],subbox.Lgwbl[_y_]);
	  for (int kv=0; kv<my_volume->GridSize[_z_]; kv++)
	    {
	      int kb=SET_PBC(kv+kc,subbox.pbc[_z_],subbox.Lgwbl[_z_]);

	      int posv=COORD_TO_INDEX(iv,jv,kv,my_volume->GridSize);
	      int posb=COORD_TO_INDEX(ib,jb,kb,subbox.Lgwbl);

	      memcpy(my_volume->Frag+posv, frag+posb, sizeof(product_data));
	      //fprintf(file," %3d %3d %3d  %10d  %10d   %10.8f\n",iv,jv,kv,posv,posb,my_volume->Frag[posv].Fmax);
	    }
	}
    }

  //fclose(file);

  my_volume->Npeaks = count_peaks(my_volume);
  my_volume->Ngroups = my_volume->Npeaks + FILAMENT + 1;
  my_volume->Groups = (group_data *)(volume_space + count_memory);

  return 0;
}

int volume2subbox(int pv, volume_data *my_volume) // e` possibile evitare di passare my_volume?
{

  /* given the index of a particle in the small volume frame
     it returns the index of the same particle in the subbox frame */

  int iv,jv,kv;
  /* from small volume index to its coordinates */
  INDEX_TO_COORD(pv,iv,jv,kv,my_volume->GridSize);
  /* from coordinates in small volume frame to those in subbox frame */
  int is = SET_PBC(iv+my_volume->Start[_x_],subbox.pbc[_x_],subbox.Lgwbl[_x_]);
  int js = SET_PBC(jv+my_volume->Start[_y_],subbox.pbc[_y_],subbox.Lgwbl[_y_]);
  int ks = SET_PBC(kv+my_volume->Start[_z_],subbox.pbc[_z_],subbox.Lgwbl[_z_]);
  /* from coordinates to subbox index */
  return COORD_TO_INDEX(is,js,ks,subbox.Lgwbl);

}


int merge_catalogs(volume_data *my_volume, int update)
{
  int i,j,k,dist,dist2;


  /* computes the smallest distance from the border of each halo */
  for (int pos=0; pos<my_volume->Npart; pos++)
    if (my_volume->Group_ID[pos]>FILAMENT)
      {
	INDEX_TO_COORD(pos,i,j,k,my_volume->GridSize);
	dist=i;
	dist2 = my_volume->GridSize[_x_]-i-1;
	dist = (dist2<dist?dist2:dist);
	dist2 = j;
	dist = (dist2<dist?dist2:dist);
	dist2 = my_volume->GridSize[_y_]-j-1;
	dist = (dist2<dist?dist2:dist);
	dist2 = k;
	dist = (dist2<dist?dist2:dist);
	dist2 = my_volume->GridSize[_z_]-k-1;
	dist = (dist2<dist?dist2:dist);

	my_volume->Groups[my_volume->Group_ID[pos]].aux2=dist;

      }

  /* sorts the global group list in order of decreasing t_peak */
  for (int i=0; i<ngroups; i++)
    *(indices+i)=i;
  qsort((void *)indices, ngroups, sizeof(int), index_compare_F);

  FILE *f1=fopen("global","w");
  for (int i=FILAMENT+1; i<ngroups; i++)
    fprintf(f1,"%10d %10d %10Ld %20.15f\n",i,indices[i],groups[indices[i]].name,groups[indices[i]].t_peak);
  fclose(f1);

  f1=fopen("volume","w");
  for (int i=FILAMENT+1; i<my_volume->Ngroups; i++)
    fprintf(f1,"%10d %10Ld %20.15f\n",i,my_volume->Groups[i].name,my_volume->Groups[i].t_peak);
  fclose(f1);

  /* matches the two group lists and create new halos if necessary */
  int p=0;
  for (int gvol=FILAMENT+1; gvol<my_volume->Ngroups; gvol++)
    {

      /* the two lists are sorted, so it's possible to scroll them at the same time 
         p is the pointer to the next element in the global list */
      while (groups[indices[p]].name   != my_volume->Groups[gvol].name && 
	     groups[indices[p]].t_peak >= my_volume->Groups[gvol].t_peak && 
	     groups[indices[p]].t_peak >= outputs.Flast)
	p++;

      if (groups[indices[p]].name == my_volume->Groups[gvol].name)
	  /* this is the case of a matching halo, aux1 will contain its location in the global catalog (of the subbox) */	  
	my_volume->Groups[gvol].aux1=indices[p];

      else
	{
	  /* in this case the halo is not present in the list */
	  printf("Questo non dovrebbe succedere! %d %Ld %f\n",gvol,my_volume->Groups[gvol].name,my_volume->Groups[gvol].t_peak);  // LEVARE
	  /* the pointer is set to a new halo that will be initialized later */
	  // CHECK that you can add one more group
	  my_volume->Groups[gvol].aux1=ngroups;
	  groups[ngroups].aux2=groups[ngroups].good=0;
	  ngroups++;
	}
    }

  FILE *fup; 
  static int first=1;

  if (!update)
    {
      fup=fopen("match","w");
      for (int gvol=FILAMENT+1; gvol<my_volume->Ngroups; gvol++)
	fprintf(fup,"  %12d  %12d\n",gvol, my_volume->Groups[gvol].aux1);
      fclose(fup);

      if (first)
	{
	  fup=fopen("update","w");
	  first=0;
	  printf("writing to file update\n");
	}
      else
	{
	  fup=fopen("match","wa");
	  printf("appending to file update\n");
	}
    }

  /* updates the groups if necessary */
  for (int gvol=FILAMENT+1; gvol<my_volume->Ngroups; gvol++)
    {

      int global=my_volume->Groups[gvol].aux1; // per semplicita` di scrittura

      // check DA LEVARE
      if (global<=FILAMENT || global>ngroups-1)
	printf("OH CACCA!!! %d %d\n",gvol,global);

      // check DA LEVARE
      if (my_volume->Groups[gvol].name != groups[global].name)
	printf("OHIBO`!!!  %d %Ld %d %Ld\n",gvol,my_volume->Groups[gvol].name,
	       global, groups[global].name);

      /* merged halos are accounted for in the merger tree of their main halo */
      if (my_volume->Groups[gvol].halo_app != gvol)
	continue;


      if (!update)
	{
	  // qui dare output di cosa succede

	  fprintf(fup,"# halo name: %Ld - local halo %d matches with global halo %d\n",groups[global].name,gvol,global);
	  fprintf(fup,"# merger tree of local halo:\n");

	  if (groups[global].halo_app != global)
	    fprintf(fup,"# CASO STRANO!\n");

	  int nbvol=0;
	  int next=gvol;
	  fprintf(fup,"         name  branch   local  global      mw    mass     mam      zpeak    zmerge \n");
	  do
	    {
	      fprintf(fup," %12Ld  %6d  %6d  %6d  %6d  %6d  %6d  %9.4f %9.4f\n",
		      my_volume->Groups[next].name,
		      nbvol++,
		      next,
		      my_volume->Groups[next].aux1,
		      my_volume->Groups[next].merged_with,
		      my_volume->Groups[next].Mass,
		      my_volume->Groups[next].mass_at_merger,
		      FTOZ(my_volume->Groups[next].t_peak),
		      FTOZ(my_volume->Groups[next].t_merge));
	      next=my_volume->Groups[next].ll;  /* next branch in the tree */		     
	    }
	  while (next != gvol);
	  fprintf(fup,"# number of branches: %d\n",nbvol);
	  fprintf(fup,"# distance from the border: %d\n",my_volume->Groups[gvol].aux2);
	  fprintf(fup,"\n");

	  int start,nbglo;
	  if (groups[global].halo_app != global)
	    {
	      fprintf(fup,"# the matching halo %d is not a main halo! it belongs to %d, name: %Ld\n",
		     global,groups[global].halo_app,groups[groups[global].halo_app].name);
	      nbglo=0;
	      start=groups[global].halo_app;
	    }
	  else
	    {
	      fprintf(fup,"# the matching halo %d is a main halo\n",global);
	      nbglo=0;
	      start=global;
	    }

	  fprintf(fup,"# merger tree of global halo %d:\n", start);

	  int nb=0;
	  next=start;
	  fprintf(fup,"         name  branch          global      mw    mass     mam      zpeak    zmerge \n");
	  do
	    {
	      fprintf(fup," %12Ld  %6d          %6d  %6d  %6d  %6d  %9.4f %9.4f\n",
		      groups[next].name,
		      nb++,
		      next,
		      groups[next].merged_with,
		      groups[next].Mass,
		      groups[next].mass_at_merger,
		      FTOZ(groups[next].t_peak),
		      FTOZ(groups[next].t_merge));
	      next=groups[next].ll;  /* next branch in the tree */		     
	    }
	  while (next != start);
	  if (groups[global].halo_app == global)
	    nbglo=nb;

	  fprintf(fup,"# number of branches: %d - of the main tree: %d\n",nbglo,nb);
	  fprintf(fup,"# distance from the border: %d\n",groups[global].aux2);
	  fprintf(fup,"\n");

	  int same_hist = (nbglo == nbvol);
	  if (same_hist)
	    {
	      int nxv=gvol;
	      int nxb=global;
	      do
		{
		  same_hist &= my_volume->Groups[nxv].aux1 == nxb;
		  nxv=my_volume->Groups[nxv].ll;
		  nxb=groups[nxb].ll;
		}
	      while (nxv != gvol);

	    }
	  int same_tree = same_hist && (my_volume->Groups[gvol].Mass == groups[global].Mass);

	  int safe_halo = (my_volume->Groups[gvol].aux2 / pow(my_volume->Groups[gvol].Mass,1./3.) > R_THR || my_volume->Groups[gvol].aux2 > A_THR);

	  int better_halo = (my_volume->Groups[gvol].aux2 > groups[global].aux2);

	  // si potrebbe aggiungere una condizione sulla perfetta consistenza dei merger tree
	  if (same_tree)
	    {
	      fprintf(fup,"these two halos are considered the same, checking history... ");

	      int nxv=gvol;
	      int nxb=global;
	      int reallythesame=1;
	      do
		{
		  reallythesame &= my_volume->Groups[nxv].aux1 == nxb;
		  nxv=my_volume->Groups[nxv].ll;
		  nxb=groups[nxb].ll;
		}
	      while (nxv != gvol);

	      if (reallythesame)
		fprintf(fup,"confirmed!\n");
	      else
		fprintf(fup,"ALAS! they are not! (this is impossible)\n");
	    }
	  else
	    {
	      fprintf(fup,"these two halos are considered different\n");

	      if (nbvol==nbglo && !same_hist)
		fprintf(fup,"WARNING! they have the same number of branches but different histories\n");

	      /* check the difference in the merger trees */
	      if (!same_hist && better_halo && safe_halo)
		{
		  fprintf(fup,"this halo would be used in the catalog\n");

		  next=gvol;	      
		  int notsingle=0;
		  do
		    {
		      fprintf(fup,"   in the global tree, halo %d belongs to halo %d\n",my_volume->Groups[next].aux1,groups[my_volume->Groups[next].aux1].halo_app);
		      notsingle+=(my_volume->Groups[next].aux1 != groups[my_volume->Groups[next].aux1].halo_app);
		      next=my_volume->Groups[next].ll;
		    }
		  while (next != gvol);

		  if (notsingle)
		    fprintf(fup,"INTERESTING CASE, there is something to understand here\n");
		  else
		    fprintf(fup,"simple case: all progenitors are single peaks in the global catalog");
		}
	    }
	  fprintf(fup,"\n");
	  fprintf(fup,"\n");
	
	  fflush(fup);
	}


      /* conditions to store the halo in the global catalog: 
         - its matching peak has not been assigned
	 - the halo has been safely reconstructed */

      /* a halo is safely reconstructed if its distance from the border is > A_THR grid points 
         and its relative distance is > R_THR */
      int safe_halo = (my_volume->Groups[gvol].aux2 / pow(my_volume->Groups[gvol].Mass,1./3.) > R_THR &&
		       my_volume->Groups[gvol].aux2 > A_THR);

      if (global==1763)
	printf("pippo\n");

      if (safe_halo)
	{
	  /* here we store the whole merger tree */
	  int next=gvol;

	  do
	    {
	      if (groups[my_volume->Groups[next].aux1].good)
		printf("MERDA! %d %d\n",next, my_volume->Groups[next].aux1);

	      /* first copy the whole structure */
	      memcpy(groups+my_volume->Groups[next].aux1, my_volume->Groups+next, sizeof(group_data));
	      /* then fix the pointers */
	      groups[my_volume->Groups[next].aux1].ll          = my_volume->Groups[my_volume->Groups[next].ll].aux1;
	      groups[my_volume->Groups[next].aux1].halo_app    = my_volume->Groups[my_volume->Groups[next].halo_app].aux1;
	      groups[my_volume->Groups[next].aux1].merged_with = my_volume->Groups[my_volume->Groups[next].merged_with].aux1;
	      /* and mark the branch as good */
	      groups[my_volume->Groups[next].aux1].good=1;

	      next=my_volume->Groups[next].ll;
	    }
	  while (next != gvol);

	  /* here we reconstruct the pointers to particles, set group_ID and 
	     set to -1 the Rmax of all the particles included in the halo, so that they are not processed again */
	  int pv=my_volume->Groups[gvol].point;
	  int ps=volume2subbox(pv,my_volume);
	  groups[global].point=ps;
	  do
	   {
	     if (my_volume->Group_ID[pv]!=gvol)
	       printf("CAZZO! %d %d %d\n",gvol,pv,my_volume->Group_ID[pv]);

	     if (pv==my_volume->Groups[gvol].bottom)
	       groups[global].bottom = ps;

	     /* sets the Rmax of that particle to -1 */
	     group_ID[ps]     = global;
	     frag[ps].Rmax    = -1;
	     pv               = my_volume->Linking_list[pv];
	     linking_list[ps] = volume2subbox(pv,my_volume);
	     ps               = linking_list[ps];
	   }
	  while (pv != my_volume->Groups[gvol].point);

	}
      else
	{
	  /* here we update only the main halo */
	  groups[global].Mass = my_volume->Groups[gvol].Mass;
	  groups[global].aux2 = my_volume->Groups[gvol].aux2;
	}
    }

  if (!update)
    fclose(fup);

  return 0;

}


int count_peaks(volume_data *v)
{

  /* counts the number of Fmax peaks down to the final redshift;
     the boundary layer is assumed to be of one particle */

  int iz,i,j,k,nn,peak_cond,i1,j1,k1;
  int Npeaks=0;

  for (iz=0; iz<v->Npart; iz++)
    {
      /* position on the local box */
      INDEX_TO_COORD(iz,i,j,k,v->GridSize);

      peak_cond=0;
      /* avoid borders */
      if ( i==0 || i==v->GridSize[_x_]-1 ) continue;
      if ( j==0 || j==v->GridSize[_y_]-1 ) continue;
      if ( k==0 || k==v->GridSize[_z_]-1 ) continue;

      /* do not count peaks that collapse late */
      if (v->Frag[iz].Fmax < outputs.Flast) continue;

      /* peak condition */
      peak_cond=1;
      for (nn=0; nn<6; nn++)
	{
	  switch (nn)
	    {
	    case 0:
	      i1=i-1;
	      j1=j;
	      k1=k;
	      break;
	    case 1:
	      i1=i+1;
	      j1=j;
	      k1=k;
	      break;
	    case 2:
	      i1=i;
	      j1=j-1;
	      k1=k;
	      break;
	    case 3:
	      i1=i;
	      j1=j+1;
	      k1=k;
	      break;
	    case 4:
	      i1=i;
	      j1=j;
	      k1=k-1;
	      break;
	    case 5:
	      i1=i;
	      j1=j;
	      k1=k+1;
	      break;
	    }

	  peak_cond &= (v->Frag[iz].Fmax > v->Frag[COORD_TO_INDEX(i1,j1,k1,v->GridSize)].Fmax);

	  if (!peak_cond)
	    break;

	}

      if (peak_cond)
	{
	  Npeaks++;
	  /* long long particle_name =  */
	  /*   COORD_TO_INDEX((long long)((i + v->start[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]), */
	  /* 		   (long long)((j + v->start[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]), */
	  /* 		   (long long)((k + v->start[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]), */
	  /* 		   MyGrids[0].GSglobal); */
	  /* printf(" picco %d: [%d %d %d], %Ld, Fmax=%f\n",Npeaks,i,j,k,particle_name,v->Frag[iz].Fmax); */
	}
    }

  return Npeaks;
}

int find_all_peaks()
{

  /* finds the Fmax peaks down to the final redshift for the whole subbox
     this works within the classic fragmentation domain decomposition 
     it is assumed that the boundary layer is only one particle */

  int iz,i,j,k,nn,npeaks_tot,peak_cond,i1,j1,k1;

  /* initializes the global number of groups */
  groups[FILAMENT].Mass=0;
  ngroups=FILAMENT+1;

  npeaks_tot=0;     /* number of groups, group 1 is the filament group */

  /* loop on all particles in the subbox */
  for (iz=0; iz<subbox.Nstored; iz++)
    {

      /* skip particles that do not collapse */
      if (frag[iz].Fmax<outputs.Flast)
	continue;

      /* position on the local box */
      INDEX_TO_COORD(iz,i,j,k,subbox.Lgwbl);

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

	  peak_cond &= (frag[iz].Fmax > frag[COORD_TO_INDEX(i1,j1,k1,subbox.Lgwbl)].Fmax);

	  if (!peak_cond)
	    break;

	}

      if (peak_cond)
	{
	  /* initializes each group as the peak, but does not assign positions and velocities */

          groups[ngroups].t_peak=frag[iz].Fmax;
	  if (params.MinHaloMass<=1)
	    groups[ngroups].t_appear=frag[iz].Fmax;
	  else
	    groups[ngroups].t_appear=-1;

          groups[ngroups].t_merge=-1;
	  groups[ngroups].Mass=1;
	  groups[ngroups].name = 
	    COORD_TO_INDEX((long long)((i + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]),
			   (long long)((j + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]),
			   (long long)((k + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]),
			   MyGrids[0].GSglobal);
	  groups[ngroups].good=0;
	  groups[ngroups].point=iz;
	  groups[ngroups].bottom=iz;
	  groups[ngroups].ll=ngroups;
	  groups[ngroups].halo_app=ngroups;

	  /* aux2 is the distance to the border, 0 means that it is not reliably reconstructed */
	  groups[ngroups].aux2=0; 

	  linking_list[iz]=iz;
	  group_ID[iz]=ngroups;

	  ngroups++;
	  npeaks_tot++;
	}
    }

  for (int i=0; i<=FILAMENT; i++)
    groups[i].name=groups[i].aux2=groups[i].t_peak=0;

  return npeaks_tot;
}


/* int create_map() */
/* { */
/*   int i,j,k,i1,i2,j1,j2,k1,k2; */

/*   /\* sets to 1 all particles in the well-resolved region plus one row for each side (without PBCs) *\/ */
/*   memset(frag_map_update, 0, subbox.maplength*sizeof(unsigned int)); */
/*   if (!subbox.pbc[_x_]) */
/*     { */
/*       i1=subbox.safe[_x_]-1; */
/*       i2=subbox.Lgrid[_x_]+subbox.safe[_x_]+1; */
/*     } */
/*   else */
/*     { */
/*       i1=0; */
/*       i2=subbox.Lgrid[_x_]; */
/*     } */
/*   if (!subbox.pbc[_y_]) */
/*     { */
/*       j1=subbox.safe[_y_]-1; */
/*       j2=subbox.Lgrid[_y_]+subbox.safe[_y_]+1; */
/*     } */
/*   else */
/*     { */
/*       j1=0; */
/*       j2=subbox.Lgrid[_y_]; */
/*     } */
/*   if (!subbox.pbc[_z_]) */
/*     { */
/*       k1=subbox.safe[_z_]-1; */
/*       k2=subbox.Lgrid[_z_]+subbox.safe[_z_]+1; */
/*     } */
/*   else */
/*     { */
/*       k1=0; */
/*       k2=subbox.Lgrid[_z_]; */
/*     } */

/*   for (i=i1; i<i2; i++) */
/*     for (j=j1; j<j2; j++) */
/*       for (k=k1; k<k2; k++) */
/* 	set_mapup_bit(i,j,k); */

/*   return 0; */
/* } */


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


/* int fragment_simpler() */
/* { */
/*   int Npeaks, Ngood; */
/*   unsigned long long mynadd[2], nadd_all[2]; */
/*   double tmp; */

/*   /\* calls the routine that builds the groups *\/ */
/*   tmp=MPI_Wtime(); */

/*   if (build_groups(Npeaks,outputs.zlast,1)) */
/*     return 1; */

/*   tmp=MPI_Wtime()-tmp; */
/*   if (!ThisTask) */
/*     printf("[%s] Fragmentation done to redshift %6.4f, cputime = %14.6f\n", */
/* 	   fdate(),outputs.zlast,tmp); */
/*   cputime.group+=tmp; */

/*   /\* cpu time for fragmentation *\/ */
/*   cputime.frag = MPI_Wtime() - cputime.frag; */

/*   // save particle group name  ***SCALAR*** */
/*   char filename[LBLENGTH]; */
/*   FILE *file; */

/*   printf("Nalloc: %d, Nstored: %d, total local size: %d\n",subbox.Nstored,subbox.Nstored,MyGrids[0].total_local_size); */
/*   sprintf(filename,"groupIDs.%s.out",params.RunFlag); */
/*   file=fopen(filename,"w"); */
/*   for (int p=0; p<MyGrids[0].total_local_size; p++) */
/*     { */
/*       long long minus99=-99; */
/*       long long minus1=-1; */
/*       if (group_ID[p]>FILAMENT) */
/* 	fwrite(&groups[group_ID[p]].name,sizeof(long long), 1, file); */
/*       else if (group_ID[p]==FILAMENT) */
/* 	fwrite(&minus1,sizeof(long long), 1, file); */
/*       else */
/* 	fwrite(&minus99,sizeof(long long), 1, file); */
/*     } */
/*   fclose(file); */

/*   if (!ThisTask) */
/*     printf("[%s] Finishing fragment, total cputime = %14.6f\n",fdate(),cputime.frag); */

/*   return 0; */
/* } */


  /*

  // trova tutti i picchi - per ora li metto in linking_list ATTENZIONE CHE INT E NON LONG LONG
  int *peak_list = linking_list;
  int Npeaks=find_all_peaks(peak_list);
  volume_data my_volume;
  printf("Found a total of %d peaks\n",Npeaks);
  printf("\n");

  FILE *file=fopen("risultato","w");
  // loop sui picchi
  for (int ipeak=0; ipeak<Npeaks; ipeak++)
    {
      if (ipeak%100==0)
	printf("%d\n",ipeak);
      // alloca la memoria per un volumetto e trascrivi le informazioni
      // provare a spostare le allocazioni in initialize_volume, quando funziona
      int size=10;
      my_volume.GridSize[0]=2*size+1;
      my_volume.GridSize[1]=2*size+1;
      my_volume.GridSize[2]=2*size+1;
      my_volume.Npart=my_volume.GridSize[0]*my_volume.GridSize[1]*my_volume.GridSize[2];
      my_volume.Frag = (product_data *)calloc(my_volume.Npart,sizeof(product_data));
      my_volume.group_ID = (int *)calloc(my_volume.Npart,sizeof(int));
      my_volume.Linking_list = (int *)calloc(my_volume.Npart,sizeof(int));

      if (initialize_volume(peak_list[ipeak],size,&my_volume))
	return 1;

      my_volume.Npeaks = count_peaks(&my_volume);
      my_volume.Ngroups = my_volume.Npeaks + FILAMENT + 1;
      my_volume.groups = (group_data *)calloc(my_volume.Ngroups, sizeof(group_data));

      // lancia build_groups, che deve lavorare su altre posizioni di memoria
      if (build_groups_in_volume(&my_volume))
	return 1;

      for (int u=FILAMENT+1; u<my_volume.Ngroups; u++)
	if (my_volume.groups[u].name==peak_list[ipeak])
	  {
	    //printf("gruppo %d: peak: %f, name:%Ld, mass:%d, t_peak:%f, t_merge:%f, merged with:%d \n",ipeak, 
	    // 	   frag[peak_list[ipeak]].Fmax, 
	    // 	   my_volume.groups[u].name, my_volume.groups[u].Mass, my_volume.groups[u].t_peak, 
	    // 	   my_volume.groups[u].t_merge, my_volume.groups[u].merged_with); 
	    double tm = (my_volume.groups[u].t_merge==-1?-1:my_volume.groups[u].t_merge-1);
	    fprintf(file," %10d  %10Ld  %4d  %10.7f  %10.7f  %10d  %10d\n",peak_list[ipeak],
		    my_volume.groups[u].name, my_volume.groups[u].Mass, my_volume.groups[u].t_peak-1,
		    tm, my_volume.groups[u].merged_with, my_volume.groups[u].mass_at_merger);

	  }
    }

  fclose(file);
  */


  /* if (1) */
  /*   { */
  /*     int SIZE=12; */
  /*     int XC=70; */
  /*     int YC=28; */
  /*     int ZC=13; */
  /*     subbox.pbc[_x_] = 0; */
  /*     subbox.pbc[_y_] = 0; */
  /*     subbox.pbc[_z_] = 0; */
  /*     subbox.Lgrid[_x_] = 2*SIZE+1; */
  /*     subbox.Lgrid[_y_] = 2*SIZE+1; */
  /*     subbox.Lgrid[_z_] = 2*SIZE+1; */
  /*     subbox.safe[_x_] = 1; */
  /*     subbox.safe[_y_] = 1; */
  /*     subbox.safe[_z_] = 1; */
  /*     subbox.Lgwbl[_x_] = subbox.Lgrid[_x_] + 2*subbox.safe[_x_]; */
  /*     subbox.Lgwbl[_y_] = subbox.Lgrid[_y_] + 2*subbox.safe[_y_]; */
  /*     subbox.Lgwbl[_z_] = subbox.Lgrid[_z_] + 2*subbox.safe[_z_]; */
  /*     subbox.start[_x_] = XC-SIZE; */
  /*     subbox.start[_y_] = YC-SIZE; */
  /*     subbox.start[_z_] = ZC-SIZE; */
  /*     subbox.stabl[_x_] = subbox.start[_x_] - subbox.safe[_x_]; */
  /*     subbox.stabl[_y_] = subbox.start[_y_] - subbox.safe[_y_]; */
  /*     subbox.stabl[_z_] = subbox.start[_z_] - subbox.safe[_z_]; */
  /*   } */

  //if (fragment())
  //if (fragment_simpler())
  //return 1;







//      else
//	{
//
//	  /* counts the number of branches for the halo in the small volume */
//	  int nbvol=0;
//	  int nbglo=0;
//	  int next=gvol;
//	  do
//	    {
//	      ++nbvol;
//	      next=my_volume->Groups[next].ll;  
//	    }
//	  while (next != gvol);
//
//	  /* /\* counts the number of branches in the halo of the global catalog *\/ */
//	  /* if (groups[global].halo_app != global) */
//	  /*   { */
//	  /*     nbglo=0;     /\* merged halos have zero branches *\/ */
//	  /*   } */
//	  /* else */
//	  /*   { */
//
//	  next=global;
//	  do
//	    {
//	      ++nbglo;
//	      next=groups[next].ll;
//	    }
//	  while (next != global);
//
//	  int same_hist = (nbglo == nbvol);
//	  if (same_hist)
//	    {
//	      int nxv=gvol;
//	      int nxb=global;
//	      do
//		{
//		  same_hist &= my_volume->Groups[nxv].aux1 == nxb;
//		  nxv=my_volume->Groups[nxv].ll;
//		  nxb=groups[nxb].ll;
//		}
//	      while (nxv != gvol);
//
//	    }
//	  
//	  /* The tree is the same if it has the same mass and identical histories */
//	  int same_tree = same_hist && (my_volume->Groups[gvol].Mass == groups[global].Mass);
//	  int safe_halo = (my_volume->Groups[gvol].aux2 / pow(my_volume->Groups[gvol].Mass,1./3.) > R_THR || my_volume->Groups[gvol].aux2 > A_THR);
//
//	  int better_halo = (my_volume->Groups[gvol].aux2 > groups[global].aux2);
//
//	  /* if the halo has changed and its distance from the border has increased then update it */
//	  if (better_halo && !same_tree)
//	    {
//	      if (safe_halo)
//		{
//		  /* here we update the whole merger tree */
//		  int next=gvol;
//		  do
//		    {
//		      /* first copy the whole structure */
//		      memcpy(groups+my_volume->Groups[next].aux1, my_volume->Groups+next, sizeof(group_data));
//		      /* then fix what should be changed */
//		      groups[my_volume->Groups[next].aux1].ll          = my_volume->Groups[my_volume->Groups[next].ll].aux1;
//		      groups[my_volume->Groups[next].aux1].halo_app    = my_volume->Groups[my_volume->Groups[next].halo_app].aux1;
//		      groups[my_volume->Groups[next].aux1].merged_with = my_volume->Groups[my_volume->Groups[next].merged_with].aux1;
//
//		      next=my_volume->Groups[next].ll;
//		    }
//		  while (next != gvol);
//		}
//	      else
//		{
//		  /* here we update only the main halo */
//		  groups[global].Mass = my_volume->Groups[gvol].Mass;
//		  groups[global].aux2 = my_volume->Groups[gvol].aux2;
//		}
//	    }
//	}
//    }
//


//  if (update)
//    {
//      /* removes ghost halos */
//      for (int i=FILAMENT+1; i<ngroups; i++)
//	groups[i].aux1=0;
//
//      for (int i=FILAMENT+1; i<ngroups; i++)
//	if (groups[i].halo_app==i)
//	  {
//	    int next=i;
//	    do
//	      {
//		groups[next].aux1=i;
//		next=groups[next].ll;
//	      }
//	    while (next != i);
//	  }
//
//      for (int i=FILAMENT+1; i<ngroups; i++)
//	if (!groups[i].aux1)
//	  {
//	    groups[i].t_merge=-1;
//	    groups[i].Mass=1;
//	    groups[i].ll=i;
//	    groups[i].halo_app=i;
//	    groups[i].aux1=0; 
//	    groups[i].aux2=0; 
//	  }
//    }

  /* /\* merges the newly produced catalog with the main one *\/ */
  /* groups[FILAMENT].Mass += my_volume->Groups[FILAMENT].Mass; */
  /* for (int ig=FILAMENT+1; ig<my_volume->Ngroups; ig++) */
  /*   { */
  /*     memcpy(groups+ngroups, my_volume->Groups+ig, sizeof(group_data)); */
  /*     ngroups++; */
  /*   } */
