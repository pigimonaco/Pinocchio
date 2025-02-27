/* ######HEADER###### */

#include "pinocchio.h"
#include <gsl/gsl_sf.h>

#define NCOUNTERS 15

unsigned long long int particle_name;
int good_particle;

#ifdef PLC
const gsl_root_fsolver_type *brent;
gsl_root_fsolver *solver;
gsl_function cPLC;
int thisgroup;
int replicate[3];
double brent_err;
#define SAFEPLC 3

double condition_PLC(PRODFLOAT);
double condition_F(double, void *);
int store_PLC(PRODFLOAT);
double find_brent(double, double);
#endif

int build_groups(int Npeaks, double zstop, int first_call)
{

  int merge[NV][NV], neigh[NV], fil_list[NV][4];
  static int iout,nstep,nstep_p;
  int nn,ifil,this_z,neigrp,nf,pos;
  int iz,i1,j1,k1,skip;
  int ig3,small,large,to_group,accgrp,ig1,ig2;
  int accrflag,nmerge,peak_cond;
  double ratio,best_ratio,d2,r2,cputmp;
  int merge_flag;
  int ibox,jbox,kbox;

  static int last_z=0;

    /*
      counters:
      0 number of peaks
      i number of particles with i neighbours (i<=6)
      7 number of accretion events
      8 number of accretion events before checking a merger
      9 number of accretion events after checking a merger
      10 number of merging events
      11 number of major mergers
      12 number of filament particles
      13 number of accreted filament particles
      14 number of good halos
    */
  static unsigned long long counters[NCOUNTERS],all_counters[NCOUNTERS];

#ifdef PLC
  static int plc_started=0, last_check_done=0;
  int irep, save, mysave, storex, storey, storez;
  static double NextF_PLC, DeltaF_PLC;
  double aa, bb, Fplc;
#endif

  if (first_call)
    {
      /* Initializations */
      ngroups=FILAMENT;           /* group list starts from FILAMENTS */
      /* filaments are not grouped! */
      for (i1=0; i1<=FILAMENT; i1++)
	{
	  groups[i1].point=-1;
	  groups[i1].bottom=-1;
	  groups[i1].good=0;
	}
      iout=0;

      for (i1=0; i1<NCOUNTERS; i1++)
	{
	  counters[i1]=0; 
	  all_counters[i1]=0;
	}


      if (!ThisTask)
	printf("[%s] Starting the fragmentation process to redshift %7.4f\n",fdate(),zstop);

      /* Calculates the number of steps required */

#ifdef CLASSIC_FRAGMENTATION
      nstep=0;
      while (frag[indices[nstep]].Fmax >= outputs.Flast)
	nstep++;
#else
      nstep=subbox.Nstored;
#endif
      nstep_p=nstep/20;

#ifdef PLC
      /* initialization of the root finding routine */
      cPLC.function = &condition_F;
      brent  = gsl_root_fsolver_brent;
      solver = gsl_root_fsolver_alloc (brent);
      DeltaF_PLC  = 0.9; // CAPIRE COME FISSARLO
      NextF_PLC   = plc.Fstart * DeltaF_PLC;
      brent_err   = 1.e-2 * params.InterPartDist;
      plc.Nstored = plc.Nstored_last = 0;
#endif

      first_call=0;
    }
  else
    {
      if (!ThisTask)
	printf("[%s] Restarting the fragmentation process to redshift %7.4f\n",fdate(),zstop);
    }


  /************************************************************************
                        START OF THE CYCLE ON POINTS
   ************************************************************************/
  for (this_z=last_z; this_z<nstep; this_z++)
    {
#ifdef CLASSIC_FRAGMENTATION
      iz=indices[this_z];
#else
      iz=this_z;
#endif

      if (frag[iz].Fmax < zstop+1.0)
	{
	  if (!ThisTask)
	    printf("[%s] Pausing fragmentation process\n",fdate());
	  last_z=this_z;
	  return 0;
	}

#ifdef PLC
      /* checks if it is time to write the stored PLC halos */
      if (plc_started && frag[iz].Fmax < NextF_PLC && frag[iz].Fmax >= plc.Fstop)
	{
	  // LEVARE le informazioni in piu`
	  if (!ThisTask)
	    printf("[%s] Syncing tasks for PLC...  Nmax=%d, Nstored=%d, Nstored_last=%d, Fmax=%f\n",fdate(),plc.Nmax, plc.Nstored, plc.Nstored_last,frag[iz].Fmax);

	  /* check whether the plc buffer is full so that it is necessary to write the plc catalog 
	     the criterion is that you need at least space for 
	     SAFEPLC times the number of halos that have been updated after last check */
	  mysave = (plc.Nmax - plc.Nstored < SAFEPLC * (plc.Nstored - plc.Nstored_last));

	  MPI_Reduce(&mysave, &save, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Bcast(&save, 1, MPI_INT, 0, MPI_COMM_WORLD);

	  // questo lo mettiamo dentro l'if dopo i test
	  if (!ThisTask)
	    printf("[%s] %d tasks require to store PLC halos at F=%6.3F...\n",fdate(),save,frag[iz].Fmax);

	  if (save)
	    {
	      /* Save PLC halos if at least one task asks to */
	      if (write_PLC(0))
		return 1;
	      plc.Nstored=0;
	    }

	  NextF_PLC *= DeltaF_PLC;
	  plc.Nstored_last = plc.Nstored;
	    
	}
#endif

      /* More initializations */

      neigrp=0;               /* number of neighbouring groups */
      nf=0;                   /* number of neighbouring filament points */
      accrflag=0;             /* if =1 all the neighbouring filaments are accreted */
      for (i1=0; i1<NV; i1++)
	neigh[i1]=0;          /* number of neighbours */

      /* grid coordinates from the indices (sub-box coordinates) */
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);
#endif

      /* skips if the point is at the border (and PBCs are not active) */
      skip=0;
      if ( !subbox.pbc[_x_] && (ibox==0 || ibox==subbox.Lgwbl[_x_]-1) ) ++skip;
      if ( !subbox.pbc[_y_] && (jbox==0 || jbox==subbox.Lgwbl[_y_]-1) ) ++skip;
      if ( !subbox.pbc[_z_] && (kbox==0 || kbox==subbox.Lgwbl[_z_]-1) ) ++skip;

      particle_name = 
	COORD_TO_INDEX((long long)((ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]),
		       (long long)((jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]),
		       (long long)((kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]),
		       MyGrids[0].GSglobal);

      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      if (!skip)
	{
	  peak_cond=1;	 
	  /* checks whether the neighbouring particles collapse later */
	  for (nn=0; nn<NV; nn++)
	    {
	      switch (nn)
		{
		case 0:
		  i1=( subbox.pbc[_x_] && ibox==0 ? subbox.Lgwbl[_x_]-1 : ibox-1 );
		  j1=jbox;
		  k1=kbox;
		  break;
		case 1:
		  i1=( subbox.pbc[_x_] && ibox==subbox.Lgwbl[_x_]-1 ? 0 : ibox+1 );
		  j1=jbox;
		  k1=kbox;
		  break;
		case 2:
		  i1=ibox;
		  j1=( subbox.pbc[_y_] && jbox==0 ? subbox.Lgwbl[_y_]-1 : jbox-1 );
		  k1=kbox;
		  break;
		case 3:
		  i1=ibox;
		  j1=( subbox.pbc[_y_] && jbox==subbox.Lgwbl[_y_]-1 ? 0 : jbox+1 );
		  k1=kbox;
		  break;
		case 4:
		  i1=ibox;
		  j1=jbox;
		  k1=( subbox.pbc[_z_] && kbox==0 ? subbox.Lgwbl[_z_]-1 : kbox-1 );
		  break;
		case 5:
		  i1=ibox;
		  j1=jbox;
		  k1=( subbox.pbc[_z_] && kbox==subbox.Lgwbl[_z_]-1 ? 0 : kbox+1 );
		  break;
		}

#ifdef CLASSIC_FRAGMENTATION
	      pos = COORD_TO_INDEX(i1,j1,k1,subbox.Lgwbl);
	      neigh[nn] = group_ID[pos];
	      peak_cond &= (frag[iz].Fmax > frag[pos].Fmax);
#else
	      pos = find_location(i1,j1,k1);
	      if (pos>=0)
		{
		  neigh[nn] = group_ID[indices[pos]];
		  peak_cond &= (frag[iz].Fmax > frag[indices[pos]].Fmax);
		}
	      else
		neigh[nn] = 0;
#endif
	      if (neigh[nn]==FILAMENT)
		{
		  neigh[nn]=0;
		  fil_list[nf][0]=i1;
		  fil_list[nf][1]=j1;
		  fil_list[nf][2]=k1;
#ifdef CLASSIC_FRAGMENTATION
		  fil_list[nf][3]=pos;
#else
		  fil_list[nf][3]=indices[pos];
#endif
		  nf++;
		}

	    }

	  /* Cleans the list of neighbouring groups */
	  clean_list(neigh);

	  /* Number of neighbouring groups */
	  for (nn=neigrp=0; nn<NV; nn++)
	    if (neigh[nn]>FILAMENT)
	      neigrp++;

	  if (neigrp>0 && good_particle) 
	    counters[neigrp]++;

#ifdef PLC
	  /* Past light cone on-the-fly reconstruction: */

	  /* is it time to reconstruct the PLC? */
	  if (frag[iz].Fmax<plc.Fstart && frag[iz].Fmax>=plc.Fstop)
	    {

	      if (!plc_started)
		{
		  plc_started=1;
		  if (!ThisTask)
		    printf("[%s] Starting PLC reconstruction, Task 0 will store at most %d halos\n",fdate(),plc.Nmax);
		  cputmp=MPI_Wtime();
		}

	      /* switching PBCs off for the check */
	      storex=subbox.pbc[_x_];
	      storey=subbox.pbc[_y_];
	      storez=subbox.pbc[_z_];
	      subbox.pbc[_x_]=subbox.pbc[_y_]=subbox.pbc[_z_]=0;

	      for (ig1=0; ig1<neigrp; ig1++)
		{
		  /* is the group good and massive enough? */
		  if (neigh[ig1] > FILAMENT &&
		      groups[neigh[ig1]].good && 
		      groups[neigh[ig1]].Mass >= params.MinHaloMass)
		    {
		      thisgroup=neigh[ig1];

		      /* loop on replications */
		      for (irep=0; irep<plc.Nreplications; irep++)
			if (!(frag[iz].Fmax > plc.repls[irep].F1 || 
			      groups[thisgroup].Flast < plc.repls[irep].F2))
			  {
			    replicate[0]=plc.repls[irep].i;
			    replicate[1]=plc.repls[irep].j;
			    replicate[2]=plc.repls[irep].k;
			    bb=condition_PLC(frag[iz].Fmax);

			    /* in this unlikely case we catch the group just on the PLC */
			    if (bb==0.0)
			      {
				if (store_PLC(frag[iz].Fmax))
				  return 1;
			      }
			    else if (bb>0.)
			      {
				/* if it is outside the PLC then check whether it has just passed */
				aa=condition_PLC(groups[thisgroup].Flast);
				if (aa<0.0)
				  {
				    /* in this case the group has passed through the PLC since last time */
				    if ( (Fplc=find_brent(groups[thisgroup].Flast,frag[iz].Fmax)) == -99.00)
				      return 1;
				    if (store_PLC(Fplc))
				      return 1;
				  }
			      }
			  }
		    }
		  /* updates the Flast field */
		  groups[neigh[ig1]].Flast=frag[iz].Fmax;
	      
		}
	      /* restoring PBCs */
	      subbox.pbc[_x_]=storex;
	      subbox.pbc[_y_]=storey;
	      subbox.pbc[_z_]=storez;

	    }
	  else if (plc.Fstart>0. && frag[iz].Fmax<plc.Fstart)
	    {
	      for (ig1=0; ig1<neigrp; ig1++)
		groups[neigh[ig1]].Flast=frag[iz].Fmax;
	    }

#endif
	}
      else
	{
	  peak_cond=0;
	  neigrp=0;
	}

      /* Is the point a peak? */
      if (peak_cond)
       {
	 /**********************************************************************
                                   FIRST CASE: PEAK
	 **********************************************************************/

	 /* New group */
        if (good_particle) 
	  counters[0]++;

        ngroups++;

	// LEVARE
        if (ngroups > Npeaks+2)
	  {
	    printf("OH MY DEAR, TOO MANY GROUPS, THIS SHOULD NOT HAPPEN!\n");
	    return 1;
	  }

	groups[ngroups].t_peak=frag[iz].Fmax;
	groups[ngroups].t_appear=-1;
	groups[ngroups].t_merge=-1;
	groups[ngroups].Pos[0]=ibox+SHIFT;
	groups[ngroups].Pos[1]=jbox+SHIFT;
	groups[ngroups].Pos[2]=kbox+SHIFT;
	groups[ngroups].Vel[0]=frag[iz].Vel[0];
	groups[ngroups].Vel[1]=frag[iz].Vel[1];
	groups[ngroups].Vel[2]=frag[iz].Vel[2];
#ifdef TWO_LPT
        groups[ngroups].Vel_2LPT[0]=frag[iz].Vel_2LPT[0];
        groups[ngroups].Vel_2LPT[1]=frag[iz].Vel_2LPT[1];
        groups[ngroups].Vel_2LPT[2]=frag[iz].Vel_2LPT[2];
#ifdef THREE_LPT
        groups[ngroups].Vel_3LPT_1[0]=frag[iz].Vel_3LPT_1[0];
        groups[ngroups].Vel_3LPT_1[1]=frag[iz].Vel_3LPT_1[1];
        groups[ngroups].Vel_3LPT_1[2]=frag[iz].Vel_3LPT_1[2];
        groups[ngroups].Vel_3LPT_2[0]=frag[iz].Vel_3LPT_2[0];
        groups[ngroups].Vel_3LPT_2[1]=frag[iz].Vel_3LPT_2[1];
        groups[ngroups].Vel_3LPT_2[2]=frag[iz].Vel_3LPT_2[2];
#endif
#endif
#ifdef RECOMPUTE_DISPLACEMENTS
	groups[ngroups].Vel_after[0]=frag[iz].Vel_after[0];
	groups[ngroups].Vel_after[1]=frag[iz].Vel_after[1];
	groups[ngroups].Vel_after[2]=frag[iz].Vel_after[2];
#ifdef TWO_LPT
        groups[ngroups].Vel_2LPT_after[0]=frag[iz].Vel_2LPT_after[0];
        groups[ngroups].Vel_2LPT_after[1]=frag[iz].Vel_2LPT_after[1];
        groups[ngroups].Vel_2LPT_after[2]=frag[iz].Vel_2LPT_after[2];
#ifdef THREE_LPT
        groups[ngroups].Vel_3LPT_1_after[0]=frag[iz].Vel_3LPT_1_after[0];
        groups[ngroups].Vel_3LPT_1_after[1]=frag[iz].Vel_3LPT_1_after[1];
        groups[ngroups].Vel_3LPT_1_after[2]=frag[iz].Vel_3LPT_1_after[2];
        groups[ngroups].Vel_3LPT_2_after[0]=frag[iz].Vel_3LPT_2_after[0];
        groups[ngroups].Vel_3LPT_2_after[1]=frag[iz].Vel_3LPT_2_after[1];
        groups[ngroups].Vel_3LPT_2_after[2]=frag[iz].Vel_3LPT_2_after[2];
#endif
#endif
#endif
	groups[ngroups].Mass=1;
	groups[ngroups].name=particle_name;
	groups[ngroups].good = good_particle;
        groups[ngroups].point = iz;
        groups[ngroups].bottom = iz;
        groups[ngroups].ll=ngroups;
        groups[ngroups].halo_app=ngroups;
#ifdef PLC
	if (frag[iz].Fmax > plc.Fstart)
	  groups[ngroups].Flast=plc.Fstart;
	else
	  groups[ngroups].Flast=frag[iz].Fmax;
#endif

        group_ID[iz]=ngroups;
        linking_list[iz]=iz;
	if (params.MinHaloMass==1)
	  {
	    groups[ngroups].t_appear=frag[iz].Fmax;
#ifdef SNAPSHOT
	    frag[iz].zacc=frag[iz].Fmax-1;
#endif
	  }

       }
     else if (neigrp==1)
       {

	 /**********************************************************************
                                  SECOND CASE: 1 GROUP
	  **********************************************************************/

	 /* if the points touches only one group, check whether to accrete the point on it */

	 condition_for_accretion(1,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[0],&d2,&r2);

	   if (d2<r2)
	     {

	       /*********************
                accretion on a group!
	        *********************/

	       if (good_particle) 
		 counters[7]++;
	       accrflag=1;
	       to_group=neigh[0];
	       accretion(to_group,ibox,jbox,kbox,iz,frag[iz].Fmax);
	     }
	   else
	     {
	       /*********
		filament!
	        *********/

	       if (good_particle) 
		 counters[12]++;
	       groups[FILAMENT].Mass++;
	       group_ID[iz]=FILAMENT;
	       linking_list[iz]=iz;
	     }

       }
     else if (neigrp>1)
       {
	 /**********************************************************************
                                   THIRD CASE: >1 GROUP
	  **********************************************************************/

	 /* In this case the point touches more than one group */

	 /*********************
          accretion on a group?
 	  *********************/

	 best_ratio=pow(10.*subbox.Lgwbl[_x_],2.0);
	 accgrp=-1;
	 for (ig1=0; ig1<neigrp; ig1++)
	   {
	     condition_for_accretion(2,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[ig1],&d2,&r2);
	       ratio=d2/r2;
	     if (ratio<1.0 && ratio<best_ratio)
	     {
	       best_ratio=ratio;
	       accgrp=ig1;
	     }
	   }

	 if (accgrp>=0)
	   {
	     if (good_particle) 
	       {
		 counters[7]++;
		 counters[8]++;
	       }
	     accrflag=1;
	     to_group=neigh[accgrp];
	     accretion(neigh[accgrp],ibox,jbox,kbox,iz,frag[iz].Fmax);
	   }

	 /* Then checks whether the groups must be merged together */

	 nmerge=0;
	 for (ig1=0; ig1<neigrp; ig1++)
	   for (ig2=0; ig2<ig1; ig2++)
	     {
	       merge[ig1][ig2]=0;
	       condition_for_merging(frag[iz].Fmax,neigh[ig1],neigh[ig2],&merge_flag);
	       if (merge_flag)
		 {
		   merge[ig1][ig2]=1;
		   nmerge++;
		 }
	     }

	 /******************
          merging of groups!
	  ******************/

	 /* The group number of the largest group is preserved */

	 if (nmerge>0)
	   {
	     for (ig1=0; ig1<neigrp; ig1++)
	       for (ig2=0; ig2<ig1; ig2++)
		 if (merge[ig1][ig2]==1 && neigh[ig1]!=neigh[ig2])
		   {
		     if (good_particle) 
		       counters[10]++;
		     if (groups[neigh[ig1]].Mass > groups[neigh[ig2]].Mass)
		       {
			 merge_groups(neigh[ig1],neigh[ig2],frag[iz].Fmax);
			 large=neigh[ig1];
			 small=neigh[ig2];
		       }
                    else
		      {
			merge_groups(neigh[ig2],neigh[ig1],frag[iz].Fmax);
			small=neigh[ig1];
			large=neigh[ig2];
		      }
		     if (to_group==small) 
		       to_group=large;
		     for (ig3=0; ig3<neigrp; ig3++)
		       if (neigh[ig3]==small) 
			 neigh[ig3]=large;

		     if (groups[large].Mass < 5*groups[small].Mass && good_particle)
		       counters[11]++;
		   }
	   }

	 /* If relevant, it tries again to accrete the particle */

	 if (accgrp==-1)
	   {
	     clean_list(neigh);

	     /* Number of neighbouring groups */
	     for (nn=neigrp=0; nn<NV; nn++)
	       if (neigh[nn]>FILAMENT)
		 neigrp++;

	     best_ratio=pow(10.*subbox.Lgwbl[_x_],2.0);
	     accgrp=-1;
	     for (ig1=0; ig1<neigrp; ig1++)
	       {
		 condition_for_accretion(3,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[ig1],&d2,&r2);
		   ratio=d2/r2;
		 if (ratio<best_ratio)
		   {
		     best_ratio=ratio;
		     accgrp=ig1;
		   }
	       }

	     if (best_ratio<1)
	       {
		 if (good_particle) 
		   {
		     counters[7]++;
		     counters[9]++;
		   }
		 accrflag=1;
		 to_group=neigh[accgrp];
		 accretion(neigh[accgrp],ibox,jbox,kbox,iz,frag[iz].Fmax);
	       }
	     else
	       {
		 /* If the point has not been accreted at all: */

		 /*********
		  filament!
 		  *********/

		 if (good_particle)
		   counters[12]++;
		 groups[FILAMENT].Mass++;
		 group_ID[iz]=FILAMENT;
		 linking_list[iz]=iz;
	       }
	   }
       }
     else
       {
	 /**********************************************************************
                                 FOURTH CASE: FILAMENTS
	  **********************************************************************/
	 if (good_particle)
	   counters[12]++;
	 groups[FILAMENT].Mass++;
	 group_ID[iz]=FILAMENT;
	 linking_list[iz]=iz;

	 /* end of cases */

       }

      /* Checks whether to accrete all the neighbouring filaments;
	 first it checks conditions for all filament particles,  
	 then it accretes those that should */

      if (accrflag && nf && !skip)
	{
	  for (ifil=0; ifil<nf; ifil++)
	    {
	      condition_for_accretion(4,fil_list[ifil][0], fil_list[ifil][1],fil_list[ifil][2],
				      fil_list[ifil][3],frag[iz].Fmax,to_group, &d2,&r2);
	      if (d2<r2)
		fil_list[ifil][3]*=-1;
	    }
	  for (ifil=0; ifil<nf; ifil++)
	    if (fil_list[ifil][3]<0)
	      {
		fil_list[ifil][3]*=-1;
		accretion(to_group, fil_list[ifil][0], fil_list[ifil][1],
			  fil_list[ifil][2],fil_list[ifil][3],frag[iz].Fmax);

		groups[FILAMENT].Mass--;

		if ( fil_list[ifil][0] >= subbox.safe[_x_] && 
		     fil_list[ifil][0] <  subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
		     fil_list[ifil][1] >= subbox.safe[_y_] && 
		     fil_list[ifil][1] <  subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
		     fil_list[ifil][2] >= subbox.safe[_z_] && 
		     fil_list[ifil][2] <  subbox.Lgwbl[_z_]-subbox.safe[_z_] )
		  {
		    counters[7]++;
		    counters[13]++;
		    counters[12]--;
		  }
	      }
      	  }

#ifdef PLC
      /* if this is the end of the cycle for PLC, it performs a last check on all halos */
      if (plc.Fstart>0 && !last_check_done &&
	  (this_z==nstep-1 || frag[iz].Fmax<plc.Fstop))
	{

	  /* first write to disc what you have */
	  if (write_PLC(0))
	    return 1;
	  plc.Nstored=0;

	  /* switching PBCs off for the check */
	  storex=subbox.pbc[_x_];
	  storey=subbox.pbc[_y_];
	  storez=subbox.pbc[_z_];
	  subbox.pbc[_x_]=subbox.pbc[_y_]=subbox.pbc[_z_]=0;

	  last_check_done=1;
	  for (ig1=FILAMENT+1; ig1<=ngroups; ig1++)
	    {
	      /* is the group alive, good and massive enough? */
	      if (groups[ig1].point >= 0 && groups[ig1].good &&
		  groups[ig1].Mass >= params.MinHaloMass)
		{
		  thisgroup=ig1;

		  /* loop on replications */
		  for (irep=0; irep<plc.Nreplications; irep++)
		    if (groups[ig1].Flast > plc.repls[irep].F2)  // CAPIRE SE VA BENE
		      {
			replicate[0]=plc.repls[irep].i;
			replicate[1]=plc.repls[irep].j;
			replicate[2]=plc.repls[irep].k;

			bb=condition_PLC(plc.Fstop);

			/* in this unlikely case we catch the group just on the PLC */
			if (bb==0.0)
			  {
			    if (store_PLC(plc.Fstop))
			      return 1;
			  }
			else if (bb>0.)
			  {

			    /* if it is outside the PLC then check whether it has just passed */
			    aa=condition_PLC(groups[ig1].Flast);
			    if (aa<0.0)
			      {
				/* in this case the group has passed through the PLC since last time */
				if ( (Fplc=find_brent(groups[ig1].Flast,plc.Fstop)) == -99.0)
				  return 1;
				if (store_PLC(Fplc))
				  return 1;
			      }
			  }
		      }
		}
	    }

	  /* restoring PBCs */
	  subbox.pbc[_x_]=storex;
	  subbox.pbc[_y_]=storey;
	  subbox.pbc[_z_]=storez;

	  if (!ThisTask)
	    printf("[%s] PLC: Last check on groups done, Task 0 stored %d halos (max:%d)\n",
		   fdate(),plc.Nstored,plc.Nmax);
	  cputime.plc += MPI_Wtime()-cputmp;

	  if (write_PLC(1))
	    return 1;
	}
#endif


      /************************************************************************
                          END OF DO-CYCLE ON COLLAPSED POINTS
       ************************************************************************/

      /* Fraction of steps done */     
      if (!ThisTask && !(this_z%nstep_p))
        printf("[%s] *** %3d%% done, F = %6.2f,  z = %6.2f\n",fdate(),
	       this_z/(nstep_p)*5,frag[iz].Fmax,
	       frag[iz].Fmax-1.0
	       );

      /* Give output if relevant */
      while (this_z==nstep-1 || frag[iz].Fmax < outputs.F[iout])
	{
	  cputmp=MPI_Wtime();

	  if (!ThisTask)
	    printf("[%s] Writing output at z=%f\n",fdate(), 
		   outputs.z[iout]);

	  fflush(stdout);
	  MPI_Barrier(MPI_COMM_WORLD);

	  if (write_catalog(iout))
	    return 1;

	  if (compute_mf(iout))
	    return 1;

#ifdef ONLYFORFIRSTBHS
	  if (iout==outputs.n-1 || iout==10)  /* Here 10 corresponds to z=8 */
	    {	      
	      if (write_histories(iout))
		return 1;
	    }
#else
	  if (iout==outputs.n-1)
	    {	      
	      if (write_histories())
		return 1;
	    }
#endif

	  cputime.io += MPI_Wtime() - cputmp;

	  iout++;
	  if (this_z==nstep-1)
	    break;
	}

    }

  /* Counters */
  for (ig1=FILAMENT+1; ig1<=ngroups; ig1++)
    if (groups[ig1].point >= 0 && groups[ig1].good) 
      counters[14]++;

  MPI_Reduce(counters, all_counters, NCOUNTERS, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      printf("Total number of peaks:                 %Lu\n",all_counters[0]);
      printf("Total number of good halos:            %Lu\n",all_counters[14]);
      printf("Particles with N neighbouring groups:  %Lu %Lu %Lu %Lu %Lu %Lu\n",
	     all_counters[1],all_counters[2],all_counters[3],
	     all_counters[4],all_counters[5],all_counters[6]);
      printf("Total number of accretion events:      %Lu\n",all_counters[ 7]);
      printf("Accretion before evaluating merger:    %Lu\n",all_counters[ 8]);
      printf("Accretion after evaluating merger:     %Lu\n",all_counters[ 9]);
      printf("Accretion of filament particles:       %Lu\n",all_counters[13]);
      printf("\n");
      printf("Global stats at the final redshift:\n");      
      printf("Total number of merger events:         %Lu\n",all_counters[10]);
      printf("Total number of major merger events:   %Lu\n",all_counters[11]);
      printf("Final number of filament particles:    %Lu\n",all_counters[12]);
      printf("Final number of particles in halos:    %Lu\n",all_counters[0]+all_counters[ 7]);
      printf("Total number of collapsed particles:   %Lu\n",all_counters[0] + all_counters[ 7] + all_counters[12]);
      printf("Total number of uncollapsed particles: %Lu\n",MyGrids[0].Ntotal - all_counters[0] - all_counters[ 7] - all_counters[12]);
      printf("\n");
    }

/* Saving group_ID for the snapshot*/

#ifdef SNAPSHOT
	groups[0].name=0;
   	groups[FILAMENT].name=FILAMENT;
   	for (iz=0; iz<subbox.Nstored; iz++){
     	frag[iz].group_ID = groups[group_ID[iz]].name;
	}	
#endif

#ifdef PLC

  gsl_root_fsolver_free (solver);

#endif

  /* Bye! */

  return 0;
}

void clean_list(int *arr)
{
  int i,j,a,done;

  for (j=1; j<NV; j++)
    {
      for (i=j-1; i>=0; i--)
	if (arr[i]==arr[j])
	  arr[j]=0;

      a=arr[j];
      done=0;

      for (i=j-1; i>=0; i--)
	{
	  if (arr[i]>0)
	    {
	      done=1;
	      break;
	    }
	  arr[i+1]=arr[i];
	}

      if (!done) 
	i=-1;
      arr[i+1]=a;
    }
}


#define RMS0 1.0

PRODFLOAT virial(int grp,PRODFLOAT F,int flag)
{
  /* Gives the SQUARE OF the "virial radius" of a group */

  PRODFLOAT r2,rlag,sigmaD;

  rlag = pow((double)grp,0.333333333333333);
  int S=Smoothing.Nsmooth-1;
  sigmaD = sqrt(Smoothing.TrueVariance[S]) * 
#ifdef SCALE_DEPENDENT
    GrowingMode((double)F-1.,Smoothing.k_GM_dens[S]); //  ATTENZIONE, QUESTO NON E` DEL TUTTO CORRETTO IN MG
#else
    GrowingMode((double)F-1.,params.k_for_GM);
#endif
  if (!flag)
    /*  merging */
    r2 = pow( f_m * pow(rlag, espo) * (sigmaD>sigmaD0 ? 1.0+(sigmaD-sigmaD0)*f_rm : 1.0) , 2.0 ) + pow( f_200 * rlag, 2.0 );
  else
    /* accretion */
    r2 = pow( f_a * pow(rlag, espo) * (sigmaD>sigmaD0 ? 1.0+(sigmaD-sigmaD0)*f_ra : 1.0) , 2.0 ) + pow( f_200 * rlag, 2.0 );

  return r2;
}

/* ACCRETION AND MERGING */

void merge_groups(int grp1,int grp2,PRODFLOAT time)
{
  /* Merges grp2 into grp1 */

  int i1;

#ifdef SNAPSHOT
  /* updates zacc of particles in case one of the halos (or both)
     was below the MinHaloMass threshold but the merged halo is above
     the threshold */
  if ( (groups[grp1].t_appear==-1 || groups[grp2].t_appear==-1)
       && (groups[grp1].Mass + groups[grp2].Mass >= params.MinHaloMass) )
    {
      if ( groups[grp1].t_appear==-1 )
	{
	  i1=groups[grp1].point;
	  while (linking_list[i1]!=groups[grp1].point)
	    {
	      frag[i1].zacc=time-1.0;
	      i1=linking_list[i1];
	    }
	  frag[i1].zacc=time-1.0;	  
	}

      if ( groups[grp2].t_appear==-1 )
	{
	  i1=groups[grp2].point;
	  while (linking_list[i1]!=groups[grp2].point)
	    {
	      frag[i1].zacc=time-1.0;
	      i1=linking_list[i1];
	    }
	  frag[i1].zacc=time-1.0;
	}      
    }
#endif

  /* updates the linking list */
  i1=groups[grp2].point;
  while (linking_list[i1] != groups[grp2].point)
    {
      group_ID[i1]=grp1;
      i1=linking_list[i1];
    }
  group_ID[i1]=grp1;

  linking_list[groups[grp1].bottom]=groups[grp2].point;
  linking_list[groups[grp2].bottom]=groups[grp1].point;
  groups[grp1].bottom=groups[grp2].bottom;
  groups[grp2].point=-1;
  groups[grp2].bottom=-1;

  /* updates substructure in the new group */
  if (groups[grp1].Mass >= params.MinHaloMass && groups[grp2].Mass >= params.MinHaloMass)
    update_history(grp1,grp2,time);

  /* updates the centre position and masses, and the existence flag */

  set_obj(grp1,time,&obj1);
  set_obj(grp2,time,&obj2);

  update(&obj1,&obj2);
  set_group(grp1,&obj1);

  if (groups[grp1].Mass >= params.MinHaloMass && groups[grp1].t_appear==-1)
    groups[grp1].t_appear=time;

  /* Bye! */
}


void update_history(int g1,int g2,PRODFLOAT time)
{
  /* UPGRADE THE SUBSTRUCTURE OF A GROUP: M(G1)>M(G2) IE G2 --> G1 */

  int old_i;

  /* if groups have no substructure */
  if (groups[g1].ll == g1 && groups[g2].ll == g2)
    {
      groups[g1].ll=g2;
      groups[g2].ll=g1;
    }
  /* if g1 has substructure but g2 is a single halo */
  else if (groups[g1].ll != g1 && groups[g2].ll == g2)
    {
      groups[g2].ll=g1;
      old_i=g1;
      while (groups[old_i].ll != g1)
        old_i=groups[old_i].ll;
      groups[old_i].ll=g2;
    }
  /* if g1 is a single halo and g2 has substructure (rare!) */
  else if (groups[g1].ll == g1 && groups[g2].ll != g2)
    {
      old_i=g2;
      while (groups[old_i].ll != g2)
	{
	  old_i=groups[old_i].ll;
	  groups[old_i].halo_app=g1;
	}
      groups[g2].halo_app=g1;
      groups[g1].ll=groups[g2].ll;
      groups[g2].ll=g1;
    }
  else
    /* both the groups have substructure */
    {
      old_i=g2;
      while (groups[old_i].ll != g2)
	{
	  old_i=groups[old_i].ll;
	  groups[old_i].halo_app=g1;
	}
      old_i=g1;
      while (groups[old_i].ll != g1)
        old_i=groups[old_i].ll;
      groups[old_i].ll=groups[g2].ll;
      groups[g2].ll=g1;
    }

  groups[g2].halo_app=g1;
  groups[g2].t_merge=time;
  groups[g2].mass_at_merger=groups[g1].Mass;
  groups[g2].merged_with=g1;
}


void accretion(int group,int i,int j,int k,int indx,PRODFLOAT F)
{
  /* Accretes the point (i,j,k) on the group */

  set_obj(group,F,&obj1);
  set_point(i,j,k,indx,F,&obj2);
  update(&obj1,&obj2);
  set_group(group,&obj1);

  /* if the group goes above the MinHaloMass threshold, set its appearing time */
  if (groups[group].Mass >= params.MinHaloMass && groups[group].t_appear==-1)
    {
      groups[group].t_appear=F;

#ifdef SNAPSHOT
      /* updates zacc for all group particles */
      int i1=groups[group].point;
      while (linking_list[i1]!=groups[group].point)
	{
	   frag[i1].zacc=F-1.0;
	   i1=linking_list[i1];
	}
      frag[i1].zacc=F-1.0;
#endif
    }

  group_ID[indx]=group;

  /* Updates the linking list */

  linking_list[groups[group].bottom]=indx;
  groups[group].bottom=indx;
  linking_list[indx]=groups[group].point;

#ifdef SNAPSHOT
  /* accretion redshift */
  if (groups[group].Mass >= params.MinHaloMass)
    frag[indx].zacc=F-1.0;
#endif
}

/* CONDITIONS */


void condition_for_accretion(int call, int i,int j,int k,int ind,PRODFLOAT Fmax, int grp,double *dd,double *rr)
{
  /* Condition for accretion */

  double dx=0,dy=0,dz=0,d2;

  *rr=virial(groups[grp].Mass,Fmax,1);
  *dd=100.* *rr;

  set_point(i,j,k,ind,Fmax,&obj1);
  set_obj(grp,Fmax,&obj2);

  dy=dz=0.0;

  dx=distance(0,&obj1,&obj2);
  d2=dx*dx;
  if (d2<*rr)
    {
      dy=distance(1,&obj1,&obj2);
      d2+=dy*dy;
      if (d2<*rr)
  	{
  	  dz=distance(2,&obj1,&obj2);
  	  d2+=dz*dz;
  	  if (d2<=*rr)
  	    *dd=d2;
  	}
    }

}


void condition_for_merging(PRODFLOAT Fmax,int grp1,int grp2,int *merge_flag)
{
  /* Checks whether two groups must be merged */

  double rvir1,rvir2,dd,rr,dx,dy,dz;

  *merge_flag=0;
  rvir1=virial(groups[grp1].Mass,Fmax,0);
  rvir2=virial(groups[grp2].Mass,Fmax,0);
  rr=(rvir1>rvir2 ? rvir1 : rvir2);

  set_obj(grp1,Fmax,&obj1);
  set_obj(grp2,Fmax,&obj2);

  dx=distance(0,&obj1,&obj2);
  dd=dx*dx;
  if (dd<rr)
    {
      dy=distance(1,&obj1,&obj2);
      dd+=dy*dy;
      if (dd<rr)
	{
	  dz=distance(2,&obj1,&obj2);
	  dd+=dz*dz;
	  if (dd<=rr)
	    *merge_flag=1;
	}
    }
}



/* DISPLACEMENTS */

void set_obj(int grp,PRODFLOAT F,pos_data *myobj)
{
  int i;

  myobj->M=groups[grp].Mass;
  myobj->z=F-1.0;

#ifdef SCALE_DEPENDENT

  /* Group velocities are averaged over group extent, so their growth  
     should be computed with a different scale */

  /* Lagrangian radius of the object */
  myobj->R=pow((double)(myobj->M)*3./4./PI,1./3.)*params.InterPartDist;  // COSTANTE GIUSTA?
  double interp=(1.-myobj->R/Smoothing.Rad_GM[0])*(double)(Smoothing.Nsmooth-1);
  interp=(interp<0.?0.:interp);
  int indx=(int)interp;
  double w=interp-(double)indx;

  /* linear interpolation of log k */
  double myk= pow(10., log10(Smoothing.k_GM_displ[indx])*(1.-w) + log10(Smoothing.k_GM_displ[indx+1])*w);

  /* growing modes for displacements */
  myobj->D=GrowingMode(myobj->z,myk)/GrowingMode(0.0,myk); 
#ifdef TWO_LPT
  myobj->D2=GrowingMode_2LPT(myobj->z,myk)*GrowingMode_2LPT(0.0,0.0)/GrowingMode_2LPT(0.0,myk);
#ifdef THREE_LPT
  myobj->D31=GrowingMode_3LPT_1(myobj->z,myk)*GrowingMode_3LPT_1(0.0,0.0)/GrowingMode_3LPT_1(0.0,myk);
  myobj->D32=GrowingMode_3LPT_2(myobj->z,myk)*GrowingMode_3LPT_2(0.0,0.0)/GrowingMode_3LPT_2(0.0,myk);
#endif
#endif

#else

  /* growing modes for displacements */
  myobj->D=GrowingMode(myobj->z,params.k_for_GM);
#ifdef TWO_LPT
  myobj->D2=GrowingMode_2LPT(myobj->z,params.k_for_GM);
#ifdef THREE_LPT
  myobj->D31=GrowingMode_3LPT_1(myobj->z,params.k_for_GM);
  myobj->D32=GrowingMode_3LPT_2(myobj->z,params.k_for_GM);
#endif
#endif

#endif

#ifdef RECOMPUTE_DISPLACEMENTS
  if (!Segment.mine)
    myobj->w = 1.0;
  else
    myobj->w = (myobj->z - Segment.z[Segment.mine-1])/
      (Segment.z[Segment.mine] - Segment.z[Segment.mine-1]);
#endif

  myobj->M=groups[grp].Mass;
  for (i=0;i<3;i++)
    {
      myobj->q[i]=groups[grp].Pos[i];
      myobj->v[i]=groups[grp].Vel[i];
#ifdef TWO_LPT
      myobj->v2[i]=groups[grp].Vel_2LPT[i];
#ifdef THREE_LPT
      myobj->v31[i]=groups[grp].Vel_3LPT_1[i];
      myobj->v32[i]=groups[grp].Vel_3LPT_2[i];
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
      myobj->v_aft[i]=groups[grp].Vel_after[i];
#ifdef TWO_LPT
      myobj->v2_aft[i]=groups[grp].Vel_2LPT_after[i];
#ifdef THREE_LPT
      myobj->v31_aft[i]=groups[grp].Vel_3LPT_1_after[i];
      myobj->v32_aft[i]=groups[grp].Vel_3LPT_2_after[i];
#endif
#endif
#endif

    }

}

void set_obj_vel(int grp,PRODFLOAT F,pos_data *myobj)
{

#ifdef SCALE_DEPENDENT

  /* this routine sets the growth rates for peculiar velocities for a group */

  double interp=(1.-myobj->R/Smoothing.Rad_GM[0])*(double)(Smoothing.Nsmooth-1);
  interp=(interp<0.?0.:interp);
  int indx=(int)interp;
  double w=interp-(double)indx;

  /* linear interpolation of log k */
  double myk=pow(10., log10(Smoothing.k_GM_vel[indx])*(1.-w) + log10(Smoothing.k_GM_vel[indx+1])*w);
  myobj->Dv = fomega(myobj->z,myk) * GrowingMode(myobj->z,myk) * GrowingMode(0.0,0.0) / GrowingMode(0.0,myk);
    //* (Smoothing.norm_GM1_vel[indx+1]*(1.-w) + Smoothing.norm_GM1_vel[indx+1]*w);
#ifdef TWO_LPT
  myobj->D2v = fomega_2LPT(myobj->z,myk) * GrowingMode_2LPT(myobj->z,myk) * GrowingMode_2LPT(0.0,0.0) / GrowingMode_2LPT(0.0,myk);
    //* (Smoothing.norm_GM2_vel[indx+1]*(1.-w) + Smoothing.norm_GM2_vel[indx+1]*w);
#ifdef THREE_LPT
  myobj->D31v = fomega_3LPT_1(myobj->z,myk) * GrowingMode_3LPT_1(myobj->z,myk) * GrowingMode_3LPT_1(0.0,0.0) / GrowingMode_3LPT_1(0.0,myk); 
    //* (Smoothing.norm_GM31_vel[indx+1]*(1.-w) + Smoothing.norm_GM31_vel[indx+1]*w);
  myobj->D32v = fomega_3LPT_2(myobj->z,myk) * GrowingMode_3LPT_2(myobj->z,myk) * GrowingMode_3LPT_2(0.0,0.0) /  GrowingMode_3LPT_2(0.0,myk);
    //* (Smoothing.norm_GM32_vel[indx+1]*(1.-w) + Smoothing.norm_GM32_vel[indx+1]*w);
#endif
#endif

#else

  double myk=params.k_for_GM;
  myobj->Dv = fomega(myobj->z,myk) * GrowingMode(myobj->z,myk);
#ifdef TWO_LPT
  myobj->D2v = fomega_2LPT(myobj->z,myk) * GrowingMode_2LPT(myobj->z,myk);
#ifdef THREE_LPT
  myobj->D31v = fomega_3LPT_1(myobj->z,myk) * GrowingMode_3LPT_1(myobj->z,myk);
  myobj->D32v = fomega_3LPT_2(myobj->z,myk) * GrowingMode_3LPT_2(myobj->z,myk);
#endif
#endif

#endif
}


void set_point(int i,int j,int k,int ind,PRODFLOAT F,pos_data *myobj)
{

  myobj->z=F-1.0;

#ifdef SCALE_DEPENDENT

  int S=Smoothing.Nsmooth-1;
  double myk=Smoothing.k_GM_displ[S];
  myobj->D=GrowingMode(myobj->z,myk) / GrowingMode(0.0,myk);
#ifdef TWO_LPT
  myobj->D2=GrowingMode_2LPT(myobj->z,myk) * GrowingMode_2LPT(0.0,0.0) / GrowingMode_2LPT(0.0,myk);
#ifdef THREE_LPT
  myobj->D31=GrowingMode_3LPT_1(myobj->z,myk) * GrowingMode_3LPT_1(0.0,0.0) / GrowingMode_3LPT_1(0.0,myk);
  myobj->D32=GrowingMode_3LPT_2(myobj->z,myk) * GrowingMode_3LPT_2(0.0,0.0) / GrowingMode_3LPT_2(0.0,myk);
#endif
#endif

#else

  double myk=params.k_for_GM;
  myobj->D=GrowingMode(myobj->z,myk);
#ifdef TWO_LPT
  myobj->D2=GrowingMode_2LPT(myobj->z,myk);
#ifdef THREE_LPT
  myobj->D31=GrowingMode_3LPT_1(myobj->z,myk);
  myobj->D32=GrowingMode_3LPT_2(myobj->z,myk);
#endif
#endif

#endif

  myobj->M=1;
  myobj->q[0]=i+SHIFT;
  myobj->q[1]=j+SHIFT;
  myobj->q[2]=k+SHIFT;
  myobj->v[0]=frag[ind].Vel[0];
  myobj->v[1]=frag[ind].Vel[1];
  myobj->v[2]=frag[ind].Vel[2];
#ifdef TWO_LPT
  myobj->v2[0]=frag[ind].Vel_2LPT[0];
  myobj->v2[1]=frag[ind].Vel_2LPT[1];
  myobj->v2[2]=frag[ind].Vel_2LPT[2];
#ifdef THREE_LPT
  myobj->v31[0]=frag[ind].Vel_3LPT_1[0];
  myobj->v31[1]=frag[ind].Vel_3LPT_1[1];
  myobj->v31[2]=frag[ind].Vel_3LPT_1[2];
  myobj->v32[0]=frag[ind].Vel_3LPT_2[0];
  myobj->v32[1]=frag[ind].Vel_3LPT_2[1];
  myobj->v32[2]=frag[ind].Vel_3LPT_2[2];
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS

  if (!Segment.mine)
    myobj->w = 1.0;
  else
    myobj->w = (myobj->z - Segment.z[Segment.mine-1])/
      (Segment.z[Segment.mine] - Segment.z[Segment.mine-1]);

  myobj->v_aft[0]=frag[ind].Vel_after[0];
  myobj->v_aft[1]=frag[ind].Vel_after[1];
  myobj->v_aft[2]=frag[ind].Vel_after[2];
#ifdef TWO_LPT
  myobj->v2_aft[0]=frag[ind].Vel_2LPT_after[0];
  myobj->v2_aft[1]=frag[ind].Vel_2LPT_after[1];
  myobj->v2_aft[2]=frag[ind].Vel_2LPT_after[2];
#ifdef THREE_LPT
  myobj->v31_aft[0]=frag[ind].Vel_3LPT_1_after[0];
  myobj->v31_aft[1]=frag[ind].Vel_3LPT_1_after[1];
  myobj->v31_aft[2]=frag[ind].Vel_3LPT_1_after[2];
  myobj->v32_aft[0]=frag[ind].Vel_3LPT_2_after[0];
  myobj->v32_aft[1]=frag[ind].Vel_3LPT_2_after[1];
  myobj->v32_aft[2]=frag[ind].Vel_3LPT_2_after[2];
#endif
#endif

#endif

}


void set_group(int grp,pos_data *myobj)
{
  int i;

  groups[grp].Mass=myobj->M;
  for (i=0;i<3;i++)
    {
      groups[grp].Pos[i]=myobj->q[i];
      groups[grp].Vel[i]=myobj->v[i];
#ifdef TWO_LPT
      groups[grp].Vel_2LPT[i]=myobj->v2[i];
#ifdef THREE_LPT
      groups[grp].Vel_3LPT_1[i]=myobj->v31[i];
      groups[grp].Vel_3LPT_2[i]=myobj->v32[i];
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS
      groups[grp].Vel_after[i]=myobj->v_aft[i];
#ifdef TWO_LPT
      groups[grp].Vel_2LPT_after[i]=myobj->v2_aft[i];
#ifdef THREE_LPT
      groups[grp].Vel_3LPT_1_after[i]=myobj->v31_aft[i];
      groups[grp].Vel_3LPT_2_after[i]=myobj->v32_aft[i];
#endif
#endif
#endif
    }
}


PRODFLOAT q2x(int i, pos_data *myobj, int pbc, double Box, int order)
{
  /* it moves an object from the Lagrangian to the Eulerian space,
   in sub-box coordinates */

  PRODFLOAT pos;

#ifndef RECOMPUTE_DISPLACEMENTS

  pos = myobj->q[i] + myobj->v[i]*myobj->D;
#ifdef TWO_LPT
  if (order>1)
    pos += myobj->D2*myobj->v2[i];
#ifdef THREE_LPT
  if (order>2)
    pos += myobj->D31*myobj->v31[i] + myobj->D32*myobj->v32[i];
#endif
#endif

#else

// CAMBIARE

  pos = myobj->q[i] + (myobj->v[i]*(1.-myobj->w) + myobj->v_aft[i]*myobj->w); // * myobj->D;
#ifdef TWO_LPT
  if (order>1)
    pos += (myobj->v2[i]*(1.-myobj->w) + myobj->v2_aft[i]*myobj->w); // * myobj->D2;
#ifdef THREE_LPT
  if (order>2)
    pos += (myobj->v31[i]*(1.-myobj->w) + myobj->v31_aft[i]*myobj->w); // * myobj->D31
        +  (myobj->v32[i]*(1.-myobj->w) + myobj->v32_aft[i]*myobj->w); // * myobj->D32;
#endif
#endif

#endif

  if (pbc)
    {
      if (pos>=Box) pos-=Box;
      if (pos< 0.0) pos+=Box;
    }

  return pos;
}



PRODFLOAT vel(int i,pos_data *myobj)
{
  /* it gives the velocity of a group in km/s */
  PRODFLOAT vv,fac;

  fac=Hubble(myobj->z)/(1.+myobj->z)*params.InterPartDist;

#ifndef RECOMPUTE_DISPLACEMENTS

  vv = myobj->v[i] * fac * myobj->Dv;
#ifdef TWO_LPT
  vv += myobj->v2[i] * fac * myobj->D2v;
#ifdef THREE_LPT
  vv += myobj->v31[i] * fac * myobj->D31v
    +   myobj->v32[i] * fac * myobj->D32v;
#endif
#endif

#else

  vv = (myobj->v[i]*(1.-myobj->w) + myobj->v_aft[i]*myobj->w) * fac * myobj->Dv;
#ifdef TWO_LPT
  vv += (myobj->v2[i]*(1.-myobj->w) + myobj->v2_aft[i]*myobj->w) * fac * myobj->D2v;
#ifdef THREE_LPT
  vv += (myobj->v31[i]*(1.-myobj->w) + myobj->v31_aft[i]*myobj->w) * fac * myobj->D31v
    +   (myobj->v32[i]*(1.-myobj->w) + myobj->v32_aft[i]*myobj->w) * fac * myobj->D32v;
#endif
#endif

#endif

  return vv;
}


PRODFLOAT distance(int i,pos_data *obj1,pos_data *obj2)
{
  PRODFLOAT d;

  /* here displacements are computed at the ORDER_FOR_GROUPS order */
  d=q2x(i,obj2,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_GROUPS)-q2x(i,obj1,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_GROUPS);

  if (subbox.pbc[i])
    {
      PRODFLOAT halfL=(PRODFLOAT)subbox.Lgwbl[i]/2.;
      if (d >  halfL) d -= subbox.Lgwbl[i];
      if (d < -halfL) d += subbox.Lgwbl[i];
    }

  return d;
}


/* ACCRETION AND MERGING */


void update(pos_data *obj1, pos_data *obj2)
{
  /* Updates the centre position of group */

  double d;
  int i;

  for (i=0;i<3;i++)
    {
      d = fabs(obj1->q[i] - obj2->q[i]);
      if (!subbox.pbc[i])
        obj1->q[i] = (obj1->q[i]*obj1->M + obj2->q[i]*obj2->M)/(double)(obj1->M + obj2->M);
      else
	{
	  PRODFLOAT halfL=subbox.Lgwbl[i]/2.;
	  /* PBC must be considered */
	  if (d <= halfL)
	    /* in this case the distance without PBC is correct */
	    obj1->q[i] = (obj1->q[i]*obj1->M + obj2->q[i]*obj2->M)/(double)(obj1->M + obj2->M);
	  else if (obj1->q[i] > halfL)
	    /* in this case xc is in the second half, and obj2->q[i] -> obj2->q[i]+Lgrid */
	    obj1->q[i] = (obj1->q[i]*obj1->M+(obj2->q[i]+subbox.Lgwbl[i])*obj2->M)/(double)(obj1->M+obj2->M);
	  else
	    /* in this case xc is in the first half, and obj2->q[i] -> obj2->q[i]-Lgrid */
	    obj1->q[i] = (obj1->q[i]*obj1->M+(obj2->q[i]-subbox.Lgwbl[i])*obj2->M)/(double)(obj1->M+obj2->M);

	  /* checks PBC again */
	  if (subbox.pbc[i] && obj1->q[i]>subbox.Lgwbl[i]) obj1->q[i]-=subbox.Lgwbl[i];
	  if (subbox.pbc[i] && obj1->q[i]<0) obj1->q[i]+=subbox.Lgwbl[i];
	}

      /* velocity */
      obj1->v[i] = (obj1->v[i]*obj1->M + obj2->v[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef TWO_LPT
      obj1->v2[i] = (obj1->v2[i]*obj1->M + obj2->v2[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef THREE_LPT
      obj1->v31[i] = (obj1->v31[i]*obj1->M + obj2->v31[i]*obj2->M)/(double)(obj1->M + obj2->M);
      obj1->v32[i] = (obj1->v32[i]*obj1->M + obj2->v32[i]*obj2->M)/(double)(obj1->M + obj2->M);
#endif
#endif

#ifdef RECOMPUTE_DISPLACEMENTS

      obj1->v_aft[i] = (obj1->v_aft[i]*obj1->M + obj2->v_aft[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef TWO_LPT
      obj1->v2_aft[i] = (obj1->v2_aft[i]*obj1->M + obj2->v2_aft[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef THREE_LPT
      obj1->v31_aft[i] = (obj1->v31_aft[i]*obj1->M + obj2->v31_aft[i]*obj2->M)/(double)(obj1->M + obj2->M);
      obj1->v32_aft[i] = (obj1->v32_aft[i]*obj1->M + obj2->v32_aft[i]*obj2->M)/(double)(obj1->M + obj2->M);
#endif
#endif

#endif
    }

  obj1->M+=obj2->M;

  /* bye! */
}

#ifdef PLC
#define MAX_ITER 100

double condition_F(double F, void *p)
{
  return condition_PLC((PRODFLOAT)F);
}

double condition_PLC(PRODFLOAT F)
{
  /* Condition for the identification of groups in past light cone */

  int i;
  double diff1,condition;

  set_obj(thisgroup,F,&obj1);

  for (i=0, condition=0.0; i<3; i++)
    {
      /* displacement is done up to ORDER_FOR_CATALOG */
      diff1 = q2x(i,&obj1,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_CATALOG) + subbox.stabl[i] - ( plc.center[i] -
					  MyGrids[0].GSglobal[i]*replicate[i] );
      condition+=diff1*diff1;
    }

  condition = sqrt(condition) - ComovingDistance((double)F-1.0)/params.InterPartDist;

  return condition;
}


int store_PLC(PRODFLOAT F)
{
  int i,ii;
  PRODFLOAT x[3];
  double rhor,theta,phi;
  static int give_message=0;

  if (plc.Nstored==plc.Nmax)
    {
      if (!give_message)
	{
	  printf("ERROR on task %d: PLC storage overshooted\n",ThisTask);
	  printf("The PLC output will be incomplete\n");
	  fflush(stdout);
	  give_message=1;
	}
      return 0;
    }

  set_obj(thisgroup,F,&obj1);
  set_obj_vel(thisgroup,F,&obj1);
  for (i=0; i<3; i++)
    {
      /* displacement is done up to ORDER_FOR_CATALOG */
      x[i] = params.InterPartDist * 
	( q2x(i,&obj1,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_CATALOG) + subbox.stabl[i] - ( plc.center[i] - MyGrids[0].GSglobal[i]*replicate[i] ) );
    }

  coord_transformation_cartesian_polar(x,&rhor,&theta,&phi);
  if (90.-theta<params.PLCAperture)
    {

      plcgroups[plc.Nstored].z    = F-1.0;
      plcgroups[plc.Nstored].Mass = groups[thisgroup].Mass;
      plcgroups[plc.Nstored].name = groups[thisgroup].name;

      for (i=0; i<3; i++)
	{
#ifdef ROTATE_BOX
	  ii=i-1;
	  if (ii==-1)
	    ii=2;
#else
	  ii=i;
#endif
	  /* displacement is done up to ORDER_FOR_CATALOG */
	  plcgroups[plc.Nstored].x[ii] = x[i];
	  plcgroups[plc.Nstored].v[ii] = vel(i,&obj1);
	}
      /* plcgroups[plc.Nstored].rhor=rhor; */
      /* plcgroups[plc.Nstored].theta=theta; */
      /* plcgroups[plc.Nstored].phi=phi; */


      int iz = (int)((plcgroups[plc.Nstored].z - params.LastzForPLC)/plc.delta_z);
      if (iz==plc.nzbins)
	iz--;
      plc.nz[iz]+=1.0;

      plc.Nstored++;
    }

  return 0;
}

void coord_transformation_cartesian_polar(PRODFLOAT *x, double *rho, double *theta, double *phi)
{
  /* transformation from cartesian coordinates to polar */

  *rho   = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  if (*rho>0)
    {
      *theta = -acos((x[0]*plc.zvers[0]+x[1]*plc.zvers[1]+x[2]*plc.zvers[2])/ *rho) * 180./PI + 90.;
      *phi   = atan2(x[0]*plc.yvers[0]+x[1]*plc.yvers[1]+x[2]*plc.yvers[2],
		     x[0]*plc.xvers[0]+x[1]*plc.xvers[1]+x[2]*plc.xvers[2]) * 180./PI;  
      if (*phi<0) *phi+=360.;
    }
  else
    {
      *theta=90.0;
      *phi=0.0;
    }
}


double find_brent(double x_hi, double x_lo)
{
  int iter, status;
  double r;

  gsl_root_fsolver_set (solver, &cPLC, x_lo, x_hi);

  iter=0;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (solver);
      r = gsl_root_fsolver_root (solver);
      x_lo = gsl_root_fsolver_x_lower (solver);
      x_hi = gsl_root_fsolver_x_upper (solver);
      /* status = gsl_root_test_interval (x_lo, x_hi, brent_err, 0.001); */
      /* if (status == GSL_SUCCESS) */
      /* 	return r; */
      
      if (fabs(condition_PLC(r)) < brent_err)
	return r;

      else
	status=GSL_CONTINUE;

    }
  while (status == GSL_CONTINUE && iter < MAX_ITER);

  printf("ERROR on task %d: find_brent could not converge - %f %f\n",ThisTask,x_hi,x_lo);
  return -99.0;

}

#endif

/* Quick construction of groups */
int quick_build_groups(int Npeaks)
{
  /* limited and quick version of build_groups:
     no PLC, no counters, no outputs */
  int merge[NV][NV], neigh[NV], fil_list[NV][4];
  int nn,ifil,neigrp,pos;
  int iz,i1,j1,k1,skip;
  int ig3,small,large,nf,to_group,accgrp,ig1,ig2;
  int accrflag,nstep,nmerge,nstep_p,peak_cond;
  double ratio,best_ratio,d2,r2;
  int merge_flag;
  int ibox,jbox,kbox;

  /* Initializations */

  ngroups=FILAMENT;           // number of groups + filaments
  groups[FILAMENT].point = groups[FILAMENT].bottom = subbox.Nstored; // filaments are not grouped!
  for (i1=0; i1<FILAMENT; i1++)  // this is probably unnecessary, but better add it
    {
      groups[i1].point=-1;
      groups[i1].good=0;
    }

  if (!ThisTask)
    printf("[%s] Starting the quick fragmentation process\n",fdate());

  /* Calculates the number of steps required */

  nstep=subbox.Nstored;
  nstep_p=nstep/5;

  /************************************************************************
                        START OF THE CYCLE ON POINTS
   ************************************************************************/
  for (iz=0; iz<nstep; iz++)
    {
      /* More initializations */

      neigrp=0;               // number of neighbouring groups
      nf=0;                   // number of neighbouring filament points
      accrflag=0;             // if =1 all the neighbouring filaments are accreted
      for (i1=0; i1<NV; i1++)
	neigh[i1]=0;          // number of neighbours

      /* grid coordinates from the indices (sub-box coordinates) */
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);

      /* skips if the point is at the border (and PBCs are not active) */
      skip=0;
      if ( !subbox.pbc[_x_] && (ibox==0 || ibox==subbox.Lgwbl[_x_]-1) ) ++skip;
      if ( !subbox.pbc[_y_] && (jbox==0 || jbox==subbox.Lgwbl[_y_]-1) ) ++skip;
      if ( !subbox.pbc[_z_] && (kbox==0 || kbox==subbox.Lgwbl[_z_]-1) ) ++skip;

      particle_name = 
	COORD_TO_INDEX((long long)((ibox + subbox.stabl[_x_] + MyGrids[0].GSglobal[_x_])%MyGrids[0].GSglobal[_x_]),
		       (long long)((jbox + subbox.stabl[_y_] + MyGrids[0].GSglobal[_y_])%MyGrids[0].GSglobal[_y_]),
		       (long long)((kbox + subbox.stabl[_z_] + MyGrids[0].GSglobal[_z_])%MyGrids[0].GSglobal[_z_]),
		       MyGrids[0].GSglobal);

      good_particle = ( ibox>=subbox.safe[_x_] && ibox<subbox.Lgwbl[_x_]-subbox.safe[_x_] && 
			jbox>=subbox.safe[_y_] && jbox<subbox.Lgwbl[_y_]-subbox.safe[_y_] && 
			kbox>=subbox.safe[_z_] && kbox<subbox.Lgwbl[_z_]-subbox.safe[_z_] );

      if (!skip)
	{
	  peak_cond=1;	 
	  /* checks whether the neighbouring particles collapse later */
	  for (nn=0; nn<NV; nn++)
	    {
	      switch (nn)
		{
		case 0:
		  i1=( subbox.pbc[_x_] && ibox==0 ? subbox.Lgwbl[_x_]-1 : ibox-1 );
		  j1=jbox;
		  k1=kbox;
		  break;
		case 1:
		  i1=( subbox.pbc[_x_] && ibox==subbox.Lgwbl[_x_]-1 ? 0 : ibox+1 );
		  j1=jbox;
		  k1=kbox;
		  break;
		case 2:
		  i1=ibox;
		  j1=( subbox.pbc[_y_] && jbox==0 ? subbox.Lgwbl[_y_]-1 : jbox-1 );
		  k1=kbox;
		  break;
		case 3:
		  i1=ibox;
		  j1=( subbox.pbc[_y_] && jbox==subbox.Lgwbl[_y_]-1 ? 0 : jbox+1 );
		  k1=kbox;
		  break;
		case 4:
		  i1=ibox;
		  j1=jbox;
		  k1=( subbox.pbc[_z_] && kbox==0 ? subbox.Lgwbl[_z_]-1 : kbox-1 );
		  break;
		case 5:
		  i1=ibox;
		  j1=jbox;
		  k1=( subbox.pbc[_z_] && kbox==subbox.Lgwbl[_z_]-1 ? 0 : kbox+1 );
		  break;
		}

	      pos = find_location(i1,j1,k1);
	      if (pos>=0)
		{
		  neigh[nn] = group_ID[indices[pos]];
		  peak_cond &= (frag[iz].Fmax > frag[indices[pos]].Fmax);
		}
	      else
		neigh[nn] = 0;

	      if (neigh[nn]==FILAMENT)
		{
		  neigh[nn]=0;
		  fil_list[nf][0]=i1;
		  fil_list[nf][1]=j1;
		  fil_list[nf][2]=k1;
		  fil_list[nf][3]=indices[pos];
		  nf++;
		}

	    }

	  /* Cleans the list of neighbouring groups */
	  clean_list(neigh);

	  /* Number of neighbouring groups */
	  for (nn=neigrp=0; nn<NV; nn++)
	    if (neigh[nn]>FILAMENT)
	      neigrp++;
	}
      else
	{
	  peak_cond=0;
	  neigrp=0;
	}


      /* Is the point a peak? */
      if (peak_cond)
       {
	 /**********************************************************************
                                   FIRST CASE: PEAK
	 **********************************************************************/
	 /* New group */
        ngroups++;
	groups[ngroups].t_peak=frag[iz].Fmax;
	groups[ngroups].t_appear=-1;
	groups[ngroups].t_merge=-1;
	groups[ngroups].Pos[0]=ibox+SHIFT;
	groups[ngroups].Pos[1]=jbox+SHIFT;
	groups[ngroups].Pos[2]=kbox+SHIFT;
	groups[ngroups].Vel[0]=frag[iz].Vel[0];
	groups[ngroups].Vel[1]=frag[iz].Vel[1];
	groups[ngroups].Vel[2]=frag[iz].Vel[2];
#ifdef TWO_LPT
        groups[ngroups].Vel_2LPT[0]=frag[iz].Vel_2LPT[0];
        groups[ngroups].Vel_2LPT[1]=frag[iz].Vel_2LPT[1];
        groups[ngroups].Vel_2LPT[2]=frag[iz].Vel_2LPT[2];
#ifdef THREE_LPT
        groups[ngroups].Vel_3LPT_1[0]=frag[iz].Vel_3LPT_1[0];
        groups[ngroups].Vel_3LPT_1[1]=frag[iz].Vel_3LPT_1[1];
        groups[ngroups].Vel_3LPT_1[2]=frag[iz].Vel_3LPT_1[2];
        groups[ngroups].Vel_3LPT_2[0]=frag[iz].Vel_3LPT_2[0];
        groups[ngroups].Vel_3LPT_2[1]=frag[iz].Vel_3LPT_2[1];
        groups[ngroups].Vel_3LPT_2[2]=frag[iz].Vel_3LPT_2[2];
#endif
#endif
	groups[ngroups].Mass=1;
	groups[ngroups].name=particle_name;
	groups[ngroups].good = good_particle;
        groups[ngroups].point = iz;
        groups[ngroups].bottom = iz;
        groups[ngroups].ll=ngroups;
        groups[ngroups].halo_app=ngroups;
        group_ID[iz]=ngroups;
        linking_list[iz]=iz;

       }
     else if (neigrp==1)
       {

	 /**********************************************************************
                                  SECOND CASE: 1 GROUP
	  **********************************************************************/
	 /* if the points touches only one group, check whether to accrete the point on it */

	 condition_for_accretion(1,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[0],&d2,&r2);
	 if (d2<r2)
	     {
	       accrflag=1;
	       to_group=neigh[0];
	       accretion(to_group,ibox,jbox,kbox,iz,frag[iz].Fmax);
	     }
	   else
	     {
	       groups[FILAMENT].Mass++;
	       group_ID[iz]=FILAMENT;
	       linking_list[iz]=iz;
	     }
       }
     else if (neigrp>1)
       {
	 /**********************************************************************
                                   THIRD CASE: >1 GROUP
	  **********************************************************************/
	 /* In this case the point touches more than one group */

	 best_ratio=pow(10.*subbox.Lgwbl[_x_],2.0);
	 accgrp=-1;
	 for (ig1=0; ig1<neigrp; ig1++)
	   {
	     condition_for_accretion(2,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[ig1],&d2,&r2);
	       ratio=d2/r2;
	     if (ratio<1.0 && ratio<best_ratio)
	     {
	       best_ratio=ratio;
	       accgrp=ig1;
	     }
	   }
	 if (accgrp>=0)
	   {
	     accrflag=1;
	     to_group=neigh[accgrp];
	     accretion(neigh[accgrp],ibox,jbox,kbox,iz,frag[iz].Fmax);
	   }
	 /* Then checks whether the groups must be merged together */
	 nmerge=0;
	 for (ig1=0; ig1<neigrp; ig1++)
	   for (ig2=0; ig2<ig1; ig2++)
	     {
	       merge[ig1][ig2]=0;
	       condition_for_merging(frag[iz].Fmax,neigh[ig1],neigh[ig2],&merge_flag);
	       if (merge_flag)
		 {
		   merge[ig1][ig2]=1;
		   nmerge++;
		 }
	     }

	 /******************
          merging of groups!
	  ******************/

	 /* The group number of the largest group is preserved */

	 if (nmerge>0)
	   {
	     for (ig1=0; ig1<neigrp; ig1++)
	       for (ig2=0; ig2<ig1; ig2++)
		 if (merge[ig1][ig2]==1 && neigh[ig1]!=neigh[ig2])
		   {
		     if (groups[neigh[ig1]].Mass > groups[neigh[ig2]].Mass)
		       {
			 merge_groups(neigh[ig1],neigh[ig2],frag[iz].Fmax);
			 large=neigh[ig1];
			 small=neigh[ig2];
		       }
                    else
		      {
			merge_groups(neigh[ig2],neigh[ig1],frag[iz].Fmax);
			small=neigh[ig1];
			large=neigh[ig2];
		      }
		     if (to_group==small) 
		       to_group=large;
		     for (ig3=0; ig3<neigrp; ig3++)
		       if (neigh[ig3]==small) 
			 neigh[ig3]=large;
		   }
	   }

	 /* If relevant, it tries again to accrete the particle */

	 if (accgrp==-1)
	   {
	     clean_list(neigh);

	     /* Number of neighbouring groups */
	     for (nn=neigrp=0; nn<NV; nn++)
	       if (neigh[nn]>FILAMENT)
		 neigrp++;

	     best_ratio=pow(10.*subbox.Lgwbl[_x_],2.0);
	     accgrp=-1;
	     for (ig1=0; ig1<neigrp; ig1++)
	       {
		 condition_for_accretion(3,ibox,jbox,kbox,iz,frag[iz].Fmax,neigh[ig1],&d2,&r2);
		   ratio=d2/r2;
		 if (ratio<best_ratio)
		   {
		     best_ratio=ratio;
		     accgrp=ig1;
		   }
	       }

	     if (best_ratio<1)
	       {
		 accrflag=1;
		 to_group=neigh[accgrp];
		 accretion(neigh[accgrp],ibox,jbox,kbox,iz,frag[iz].Fmax);
	       }
	     else
	       {
		 /* If the point has not been accreted at all: */
		 groups[FILAMENT].Mass++;
		 group_ID[iz]=FILAMENT;
		 linking_list[iz]=iz;
	       }
	   }
       }
     else
       {
	 /**********************************************************************
                                 FOURTH CASE: FILAMENTS
	  **********************************************************************/
	 groups[FILAMENT].Mass++;
	 group_ID[iz]=FILAMENT;
	 linking_list[iz]=iz;

	 /* end of cases */
       }

      /* Checks whether to accrete all the neighbouring filaments;
	 first it checks conditions for all filament particles,  
	 then it accretes those that should */
      if (accrflag && nf && !skip)
	{
	  for (ifil=0; ifil<nf; ifil++)
	    {
	      condition_for_accretion(4,fil_list[ifil][0], fil_list[ifil][1],fil_list[ifil][2],
				      fil_list[ifil][3],frag[iz].Fmax,to_group, &d2,&r2);
	      if (d2<r2)
		fil_list[ifil][3]*=-1;
	    }
	  
	  for (ifil=0; ifil<nf; ifil++)
	    if (fil_list[ifil][3]<0)
	      {
		fil_list[ifil][3]*=-1;
		accretion(to_group, fil_list[ifil][0], fil_list[ifil][1],
			  fil_list[ifil][2],fil_list[ifil][3],frag[iz].Fmax);
		groups[FILAMENT].Mass--;
	      }
      	  }

      /************************************************************************
                          END OF DO-CYCLE ON COLLAPSED POINTS
       ************************************************************************/

      /* Fraction of steps done */
      if (!ThisTask && !(iz%nstep_p))
        printf("[%s] *** %3d%% done, F = %6.2f,  z = %6.2f\n",fdate(),
	       iz/(nstep_p)*20,frag[iz].Fmax,
	       frag[iz].Fmax-1.0
	       );

    }

  return 0;
}


int update_map(unsigned int *nadd)
{

  int group,ig,jg,kg,size,size2,rr,i,j,k,i1,j1,k1;
  memset(frag_map_update, 0, subbox.maplength*sizeof(unsigned int));
  nadd[0]=nadd[1]=0;

  for (group=FILAMENT+1; group<ngroups; group++)
    {
      ig=(int)(groups[group].Pos[0]+0.5);
      jg=(int)(groups[group].Pos[1]+0.5);
      kg=(int)(groups[group].Pos[2]+0.5);
      size=(int)(params.BoundaryLayerFactor*pow((double)groups[group].Mass/4.188790205,0.333333333333333)+0.5);
      size2=size*size;
      
      for (i1=ig-size; i1<ig+size; i1++)
	{
	  if (i1<0 || i1>=subbox.Lgwbl[_x_])
		{
		  if (subbox.pbc[_x_])
		    i = ( i1<0 ? i1+subbox.Lgwbl[_x_] : i1-subbox.Lgwbl[_x_] );
		  else
		    i = -1;
		}
	      else
		i=i1;

	  for (j1=jg-size; j1<jg+size; j1++)
	    {
	      if (j1<0 || j1>=subbox.Lgwbl[_y_])
		{
		  if (subbox.pbc[_y_])
		    j = ( j1<0 ? j1+subbox.Lgwbl[_y_] : j1-subbox.Lgwbl[_y_] );
		  else
		    j = -1;
		}
	      else
		j=j1;
	      
	      for (k1=kg-size; k1<kg+size; k1++)
		{
		  if (k1<0 || k1>=subbox.Lgwbl[_z_])
		    {
		      if (subbox.pbc[_z_])
			k = ( k1<0 ? k1+subbox.Lgwbl[_z_] : k1-subbox.Lgwbl[_z_] );
		      else
			k = -1;
		    }
		  else
		    k=k1;

		  if (i<0 || j<0 || k<0)
		    {
		      nadd[1]++;
		      continue;
		    }

		  if (!get_map_bit(i,j,k))
		    {
		      rr=(i1-ig)*(i1-ig)+(j1-jg)*(j1-jg)+(k1-kg)*(k1-kg);
		      if (rr<=size2)
			{
			  set_mapup_bit(i,j,k);
			  nadd[0]++;
			}
		    }
		}
	    }
	}
    }

  return 0;
}

