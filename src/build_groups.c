/*****************************************************************
 *                        PINOCCHIO  V5.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco, Luca Tornatore, Marius Lepinzan, Chiara Moretti, Giuliano Taffoni
 Copyright (C) 2024
 
 github repository: https://github.com/pigimonaco/Pinocchio
 
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




/* 
   This file contains code to construct groups (dark matter halos)
   from the list of collapsed particles.
   The main function is build_groups, that performs the group
   construction down to redshift zstop.
   quick_build_groups implements a quick version of build_groups, used
   to determine which part of the boundary region should be
   communicated to perform group construction.
*/

#include "pinocchio.h"
#include <gsl/gsl_sf.h>

#define NCOUNTERS 15

unsigned long long int particle_name;
int good_particle;
pos_data obj1,obj2;

void set_weight(pos_data *);

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

  /* The algorithm for group construction runs as follows:
     loop on all collapsed particle
     + for each particle check its neighbours
     + if it is a peak of Fmax, make it a one-particle halo
     + if it touches only one halo:
       -> check if it should be accreted
       -> if it should, accrete it
       -> otherwise, tag it as filament
     + if it touches more than one halo:
       -> check if it should be accreted to one halo
       -> if it should, choose the one it gets nearest to
       -> then check if all the halo pairs should merge
       -> if they should, merge them
       -> if the particle was not accreted before, re-check it
       -> otherwise, tag it as filament
     + if it touches only filaments, tag it as filament
   */


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
    counters of various events:
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
      /* this part contains initializations that should be done at
         first call */
      ngroups=FILAMENT;           /* group list starts from FILAMENTS */
      /* filaments are not grouped */
      for (i1=0; i1<=FILAMENT; i1++)
        {
          groups[i1].point=-1;
          groups[i1].bottom=-1;
          groups[i1].good=0;
        }
      iout=0;

      /* sets the counter to zero */
      for (i1=0; i1<NCOUNTERS; i1++)
        {
          counters[i1]=0; 
          all_counters[i1]=0;
        }


      if (!ThisTask)
        printf("[%s] Starting the fragmentation process to redshift %7.4f\n",fdate(),zstop);

      /* nstep is the number of collapsed particles that will be
         checked by the code. In classic fragmentation this is the
         number of particles that collapse by the end of the run. In
         default fragmentation this selection is performed at
         distribution time, so this is the number of stored particles.
      */
#ifdef CLASSIC_FRAGMENTATION
      nstep=0;
      while (frag[indices[nstep]].Fmax >= outputs.Flast)
        nstep++;
#else
      nstep=subbox.Nstored;
#endif
      nstep_p=nstep/20;

#ifdef PLC
      /* initialization of the root finding routine used by the PLC code */
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
                    START OF THE CYCLE ON COLLAPSED PARTICLES
   ************************************************************************/
  for (this_z=last_z; this_z<nstep; this_z++)
    {
      /* In classic fragmentation the particles must be addressed in
         order of collapse time, while in default fragmentation they
         are already in that order */
#ifdef CLASSIC_FRAGMENTATION
      iz=indices[this_z];
#else
      iz=this_z;
#endif


#ifdef PLC
      /* here it checks if it is time to write the stored PLC halos.
         If it is the case, it waits for all the tasks to get to the
         same point and writes the catalog up to this point
       */
      if (plc_started && frag[iz].Fmax < NextF_PLC && frag[iz].Fmax >= plc.Fstop)
        {
          if (!ThisTask)
            printf("[%s] Syncing tasks for PLC...  Nmax=%d, Nstored=%d, Nstored_last=%d, Fmax=%f\n",fdate(),plc.Nmax, plc.Nstored, plc.Nstored_last,frag[iz].Fmax);

          /* each task checks if the plc buffer is full and it is
             necessary to write the plc catalog. The criterion is that
             you need at least space for SAFEPLC times the number of
             halos that have been updated after last check */
          mysave = (plc.Nmax - plc.Nstored < SAFEPLC * (plc.Nstored - plc.Nstored_last));

          MPI_Reduce(&mysave, &save, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
          MPI_Bcast(&save, 1, MPI_INT, 0, MPI_COMM_WORLD);

          if (!ThisTask)
            printf("[%s] %d tasks require to store PLC halos at F=%6.3F...\n",fdate(),save,frag[iz].Fmax);

          if (save)
            {
              /* Save PLC halos if at least one task asks to */
              if (write_PLC(0))
                return 1;
              plc.Nstored=0;
            }

          /* next update will be at this time */
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

      /*******************************************/
      /* PARTICLE COORDINATES AND PEAK CONDITION */
      /*******************************************/
#ifdef CLASSIC_FRAGMENTATION
      INDEX_TO_COORD(iz,ibox,jbox,kbox,subbox.Lgwbl);
#else
      INDEX_TO_COORD(frag_pos[iz],ibox,jbox,kbox,subbox.Lgwbl);
#endif

      /* skips the peak condition if the point is at the border (and PBCs are not active) */
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
              /* coordinates of the neighbouring particle */
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

              /* accessing the neighbouring particle differs in
                 classic and default fragmentation: in classic
                 fragmentation the information on the neighbour is
                 immediately obtained, in standard fragmentation the
                 neighbour must be seeked in the particle list */
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

              /* neighbouring filaments are stored separately */
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

          /* Cleans the list of neighbouring groups removing duplicates */
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

              /* check if this is the first call */
              if (!plc_started)
                {
                  plc_started=1;
                  if (!ThisTask)
                    printf("[%s] Starting PLC reconstruction, Task 0 will store at most %d halos\n",fdate(),plc.Nmax);
                  cputmp=MPI_Wtime();
                }

              /* PBCs are switched off for this check */
              storex=subbox.pbc[_x_];
              storey=subbox.pbc[_y_];
              storez=subbox.pbc[_z_];
              subbox.pbc[_x_]=subbox.pbc[_y_]=subbox.pbc[_z_]=0;

              /* the check is performed on all neighbouring groups */
              for (ig1=0; ig1<neigrp; ig1++)
                {
                  /* is the group good and massive enough? */
                  if (neigh[ig1] > FILAMENT &&
                      groups[neigh[ig1]].good && 
                      groups[neigh[ig1]].Mass >= params.MinHaloMass)
                    {
                      thisgroup=neigh[ig1];

                      /* loop on replications */
                      /* this loop may be threaded or ported to GPUs,
                         though it's relatively fast */
                      for (irep=0; irep<plc.Nreplications; irep++)
                        /* checks that the redshift falls in the
                           range of the replication */
                        if (!(frag[iz].Fmax > plc.repls[irep].F1 || 
                              groups[thisgroup].Flast < plc.repls[irep].F2))
                          {
                            replicate[0]=plc.repls[irep].i;
                            replicate[1]=plc.repls[irep].j;
                            replicate[2]=plc.repls[irep].k;
                            /* this computes the difference between
                               the group distance from the observer
                               and its comoving distance at that z */
                            bb=condition_PLC(frag[iz].Fmax);

                            /* in this unlikely case we catch the
                               group just on the PLC */
                            if (bb==0.0)
                              {
                                if (store_PLC(frag[iz].Fmax))
                                  return 1;
                              }
                            else if (bb>0.)
                              {
                                /* if it is outside the PLC then check
                                   whether it has just passed since
                                   last check */
                                aa=condition_PLC(groups[thisgroup].Flast);
                                if (aa<0.0)
                                  {
                                    /* in this case the group has
                                       passed through the PLC since
                                       last check, then compute its
                                       time of PLC cross */
                                    if ( (Fplc=find_brent(groups[thisgroup].Flast,frag[iz].Fmax)) == -99.00)
                                      return 1;
                                    /* and store it */
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
              /* if nothing must be done, only update the Flast field */
              for (ig1=0; ig1<neigrp; ig1++)
                groups[neigh[ig1]].Flast=frag[iz].Fmax;
            }

#endif
        }
      else
        {
          /* this closes the if (!skip) condition above: if the
             particle is at the border of the domain it will not be a
             peak and will have no neighbours */
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

          /* this paranoid check will raise an error in case of major bugs */
          if (ngroups > Npeaks+2)
            {
              printf("OH MY DEAR, TASK %d FOUND TOO MANY GROUPS AT z=%f, STEP %d OF %d, THIS SHOULD NOT HAPPEN!\n",ThisTask,frag[iz].Fmax-1.0,this_z,nstep);
              return 1;
            }

          /* sets the properties of a one-particle group */
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
          groups[ngroups].Vel_prev[0]=frag[iz].Vel_prev[0];
          groups[ngroups].Vel_prev[1]=frag[iz].Vel_prev[1];
          groups[ngroups].Vel_prev[2]=frag[iz].Vel_prev[2];
#ifdef TWO_LPT
          groups[ngroups].Vel_2LPT_prev[0]=frag[iz].Vel_2LPT_prev[0];
          groups[ngroups].Vel_2LPT_prev[1]=frag[iz].Vel_2LPT_prev[1];
          groups[ngroups].Vel_2LPT_prev[2]=frag[iz].Vel_2LPT_prev[2];
#ifdef THREE_LPT
          groups[ngroups].Vel_3LPT_1_prev[0]=frag[iz].Vel_3LPT_1_prev[0];
          groups[ngroups].Vel_3LPT_1_prev[1]=frag[iz].Vel_3LPT_1_prev[1];
          groups[ngroups].Vel_3LPT_1_prev[2]=frag[iz].Vel_3LPT_1_prev[2];
          groups[ngroups].Vel_3LPT_2_prev[0]=frag[iz].Vel_3LPT_2_prev[0];
          groups[ngroups].Vel_3LPT_2_prev[1]=frag[iz].Vel_3LPT_2_prev[1];
          groups[ngroups].Vel_3LPT_2_prev[2]=frag[iz].Vel_3LPT_2_prev[2];
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

          /* Then checks by pairs whether the groups must be merged together */

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

	  /****************/
          /* end of cases */
	  /****************/

        }

      /* If the particle was accreted somewhere, checks whether to
         accrete to the same halo all the neighbouring filaments;
         first it checks conditions for all filament particles, then
         it accretes those that should */

      if (accrflag && nf && !skip)
        {
	  /* loop on filament particles */
          for (ifil=0; ifil<nf; ifil++)
            {
              condition_for_accretion(4,fil_list[ifil][0], fil_list[ifil][1],fil_list[ifil][2],
                                      fil_list[ifil][3],frag[iz].Fmax,to_group, &d2,&r2);
	      /* tags filaments that should be accreted (the condition
		 should be checked without changing the halo) */
              if (d2<r2)
                fil_list[ifil][3]*=-1;
            }
	  /* accrete all filament particles */
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
      /* if this is the end of the cycle for PLC, perform a last check
	 on all halos to see if some has passed through the PLC in
	 the meantime */
      /* Note: this part is largely a replication of the PLC code above */
      if (plc.Fstart>0 && !last_check_done &&
          (this_z==nstep-1 || frag[iz].Fmax<plc.Fstop))
        {

          /* first write to disc what you have */
          if (write_PLC(0))
            return 1;
          plc.Nstored=0;

	  /* PBCs are switched off for this check */
          storex=subbox.pbc[_x_];
          storey=subbox.pbc[_y_];
          storez=subbox.pbc[_z_];
          subbox.pbc[_x_]=subbox.pbc[_y_]=subbox.pbc[_z_]=0;

          last_check_done=1;
	  /* the check is performed on ALL groups */
          for (ig1=FILAMENT+1; ig1<=ngroups; ig1++)
            {
              /* is the group alive, good and massive enough? */
              if (groups[ig1].point >= 0 && groups[ig1].good &&
                  groups[ig1].Mass >= params.MinHaloMass)
                {
                  thisgroup=ig1;

                  /* loop on replications */
		  /* this loop may be threaded or ported to GPUs,
		     though it's relatively fast */
                  for (irep=0; irep<plc.Nreplications; irep++)
		    /* checks that the redshift falls in the
		       range of the replication */
                    if (groups[ig1].Flast > plc.repls[irep].F2)
                      {
                        replicate[0]=plc.repls[irep].i;
                        replicate[1]=plc.repls[irep].j;
                        replicate[2]=plc.repls[irep].k;
			/* this computes the difference between the
			   group distance from the observer and its
			   comoving distance at that z */
                        bb=condition_PLC(plc.Fstop);

                        /* in this unlikely case we catch the group
			   just on the PLC */
                        if (bb==0.0)
                          {
                            if (store_PLC(plc.Fstop))
                              return 1;
                          }
                        else if (bb>0.)
                          {
			    /* if it is outside the PLC then check
			       whether it has just passed since
			       last check */
                            aa=condition_PLC(groups[ig1].Flast);
                            if (aa<0.0)
                              {
                                /* in this case the group has passed
				   through the PLC since last time */
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
                                  CLOSING THE CYCLE
       ************************************************************************/

      /* prints the fraction of steps done */
      if (!ThisTask && !(this_z%nstep_p))
        printf("[%s] *** %3d%% done, F = %6.2f,  z = %6.2f\n",fdate(),
               this_z/(nstep_p)*5,frag[iz].Fmax,
               frag[iz].Fmax-1.0);

      /* Write output if relevant. 
         The while cycle is here because it may happen that some
         outputs are requested when there are very few collapsed
         particles, and if no particle collapses between two catalogs
         one catalog may not be written; better to have an empty
         catalog, or two identical ones */
      while (this_z==nstep-1 || frag[iz].Fmax < outputs.F[iout])
        {
          cputmp=MPI_Wtime();

          if (!ThisTask)
            printf("[%s] Writing output at z=%f\n",fdate(), 
                   outputs.z[iout]);

          fflush(stdout);
          MPI_Barrier(MPI_COMM_WORLD);

	  /* halo catalog */
          if (write_catalog(iout))
            return 1;

	  /* halo mass function */
          if (compute_mf(iout))
            return 1;

	  /* merger histories, only at the last redshift */
          if (iout==outputs.n-1)
            {         
              if (write_histories())
                return 1;
            }

          cputime.io += MPI_Wtime() - cputmp;

	  /* if we are at the end, just break out of the while cycle */
          iout++;
          if (this_z==nstep-1)
            break;
        }


      if (this_z!=nstep-1 && frag[iz].Fmax < zstop+1.0)
        {
          if (!ThisTask)
            printf("[%s] Pausing fragmentation process\n",fdate());
          last_z=this_z+1;
          return 0;
        }

      /************************************************************************
                          END OF DO-CYCLE ON COLLAPSED POINTS
       ************************************************************************/
    }

  /* Counters */
  for (ig1=FILAMENT+1; ig1<=ngroups; ig1++)
    if (groups[ig1].point >= 0 && groups[ig1].good) 
      counters[14]++;

  MPI_Reduce(counters, all_counters, NCOUNTERS, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

  /* statistics */
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


#ifdef SNAPSHOT
  /* Saving group_ID in the frag structure for the snapshot*/
  groups[0].name=0;
  groups[FILAMENT].name=FILAMENT;
  for (iz=0; iz<subbox.Nstored; iz++)
    frag[iz].group_ID = groups[group_ID[iz]].name;
#endif

#ifdef PLC

  /* close the root solver */
  gsl_root_fsolver_free (solver);

#endif

  /* Bye! */

  return 0;
}

void clean_list(int *arr)
{
  /* this function removes replications in the neighbour list and
     moves all non-empty values at the beginning */
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



PRODFLOAT virial(int mass,PRODFLOAT F,int flag)
{
  /* This function is at the root of the decision of whether particles
     should be accreted and halos should be merged. For a group of a
     given mass (in number of particles) it returns the SQUARE of its
     "virial radius", that is the capture radius for accretion and
     merging events. 
     flag=1 for accretion, flag=1 for merging.
  */

  PRODFLOAT r2,rlag,sigmaD;

  /* the Lagrangian radius in grid units is just the cubic root of the
  halo mass */
  rlag = pow((double)mass,0.333333333333333);
  int S=Smoothing.Nsmooth-1;
  /* this is the linear density rms on the grid */
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



   /*-----------------MODIFIED GRAVITY WITH f(R)-----------------*/

  // DEVO CAPIRE COS'E` QUESTO CODICE QUI SOTTO!!!

#ifdef SALTALA
  // FILE *output_file_rlag = fopen("output_rlag_values.txt", "w");
  int interpolation_Done = 0;

  for (int iradius = 0; iradius <= S  ; iradius++) {
    if (rlag * params.InterPartDist >= Smoothing.Radius[iradius]) {

      // Check if interpolation has already been done for this rlag
      if (interpolation_Done == 0) { 

	// Perform reverse linear interpolation to find the k value
	//printf("your rlag is smaller than this smoothing radius: %.6f\n", Smoothing.Radius[iradius - 1]);
	double r_lower = Smoothing.Radius[iradius];
	double r_upper = Smoothing.Radius[iradius - 1];
	double k_lower = Smoothing.k_GM_dens[iradius];
	double k_upper = Smoothing.k_GM_dens[iradius - 1];

	// Reverse linear interpolation formula
	// Risultato identico all'interpolazione lineare di GSL. Si potrebbe usare per l'interpolazione nei tempi di collasso? Così si può portare su GPU senza GSL !
	// double k_interpolated = k_lower + (k_upper - k_lower) * ((rlag * params.InterPartDist - r_lower) / (r_upper - r_lower));
                                
	// Define arrays for r and k values
	double r_values[2] = {r_lower, r_upper};
	double k_values[2] = {k_lower, k_upper};

	// Perform linear interpolation using GSL
	gsl_interp *interp = gsl_interp_alloc(gsl_interp_linear, 2);
	gsl_interp_init(interp, r_values, k_values, 2);
                
	double r_value = rlag * params.InterPartDist;
	double k_interpolated = gsl_interp_eval(interp, r_values, k_values, r_value, NULL);

	//Print data to the terminal for debugging in an organized layout
	//printf("rlag: %.6f\t r_lower: %.6f\t r_upper: %.6f\t k_lower: %.6f\t k_upper: %.6f\t k_interpolated: %.6f\n",
	//rlag * params.InterPartDist, r_lower, r_upper, k_lower, k_upper, k_interpolated);
	//printf("\n");
	// Set the flag to indicate that interpolation has been done for this rlag
	interpolation_Done = 1;
	// gsl_interp_free(interp);
      }
      
    }
    
  }
        
  // fclose(output_file_rlag);

#endif

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

  /* joins the linking lists */
  linking_list[groups[grp1].bottom]=groups[grp2].point;
  linking_list[groups[grp2].bottom]=groups[grp1].point;
  groups[grp1].bottom=groups[grp2].bottom;
  groups[grp2].point=-1;
  groups[grp2].bottom=-1;

  /* updates merger history in the new group */
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
  /* Upgrade the history of a group: g1 flows into g2 */

  int old_i;

  /* if groups have no branches */
  if (groups[g1].ll == g1 && groups[g2].ll == g2)
    {
      groups[g1].ll=g2;
      groups[g2].ll=g1;
    }
  /* if g1 has branches but g2 is a single halo */
  else if (groups[g1].ll != g1 && groups[g2].ll == g2)
    {
      groups[g2].ll=g1;
      old_i=g1;
      while (groups[old_i].ll != g1)
        old_i=groups[old_i].ll;
      groups[old_i].ll=g2;
    }
  /* if g1 is a single halo and g2 has branches */
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
    /* both the groups have branches */
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
  /* Checks whether a particle should be accreted on a group */

  double dx=0,dy=0,dz=0,d2;

  /* compute the capture radius of the group */
  *rr=virial(groups[grp].Mass,Fmax,1);
  *dd=100.* *rr;

  set_point(i,j,k,ind,Fmax,&obj1);
  set_obj(grp,Fmax,&obj2);

  dy=dz=0.0;

  /* to shorten the calculation, distances are checked dim by dim */
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
  /* sets all the quantities needed to handle a group */

  myobj->M=groups[grp].Mass;
  myobj->z=F-1.0;

#ifdef SCALE_DEPENDENT
  /* Group velocities are averaged over group particles, so their growth  
     should be computed on a different scale */

  /* Lagrangian radius of the object */
  // QUESTA INTERPOLAZIONE SAREBBE DA CONTROLLARE
  myobj->R=pow((double)(myobj->M)*3./4./PI,1./3.)*params.InterPartDist;  // COSTANTE GIUSTA?
  double interp=(1.-myobj->R/Smoothing.Rad_GM[0])*(double)(Smoothing.Nsmooth-1);
  interp=(interp<0.?0.:interp);
  int indx=(int)interp;
  double w=interp-(double)indx;

  /* the scale to compute the growth rate is obtaine by linear
     interpolation of log k in time */
  myobj->myk= pow(10., log10(Smoothing.k_GM_displ[indx])*(1.-w) + log10(Smoothing.k_GM_displ[indx+1])*w);
#else
  myobj->myk=params.k_for_GM;
#endif

  set_weight(myobj);

  myobj->M=groups[grp].Mass;
  for (int i=0;i<3;i++)
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
      myobj->v_prev[i]=groups[grp].Vel_prev[i];
#ifdef TWO_LPT
      myobj->v2_prev[i]=groups[grp].Vel_2LPT_prev[i];
#ifdef THREE_LPT
      myobj->v31_prev[i]=groups[grp].Vel_3LPT_1_prev[i];
      myobj->v32_prev[i]=groups[grp].Vel_3LPT_2_prev[i];
#endif
#endif
#endif

    }

}


void set_weight(pos_data *myobj)
{
  /* sets the weight in a pos_data object for interpolating displacements */
  if (!ScaleDep.myseg)
    {
      /* at the first fragmentation segment, the interpolation weight
	 is just the growing mode at redshift F-1 divided by the "final" one */
      myobj->w   = GrowingMode(       myobj->z, myobj->myk) / GrowingMode(       ScaleDep.z[ScaleDep.myseg], myobj->myk);
#ifdef TWO_LPT
      myobj->w2  = GrowingMode_2LPT(  myobj->z, myobj->myk) / GrowingMode_2LPT(  ScaleDep.z[ScaleDep.myseg], myobj->myk);
#ifdef THREE_LPT
      myobj->w31 = GrowingMode_3LPT_1(myobj->z, myobj->myk) / GrowingMode_3LPT_1(ScaleDep.z[ScaleDep.myseg], myobj->myk);
      myobj->w32 = GrowingMode_3LPT_2(myobj->z, myobj->myk) / GrowingMode_3LPT_2(ScaleDep.z[ScaleDep.myseg], myobj->myk);
#endif
#endif
    }
  else
    {
      /* after, the weight linearly interpolates between the two redshifts */
      myobj->w   = (GrowingMode(                         myobj->z, myobj->myk) - GrowingMode(       ScaleDep.z[ScaleDep.myseg-1], myobj->myk)) /
                   (GrowingMode(       ScaleDep.z[ScaleDep.myseg], myobj->myk) - GrowingMode(       ScaleDep.z[ScaleDep.myseg-1], myobj->myk) );
#ifdef TWO_LPT
      myobj->w2  = (GrowingMode_2LPT(                    myobj->z, myobj->myk) - GrowingMode_2LPT(  ScaleDep.z[ScaleDep.myseg-1], myobj->myk)) /
                   (GrowingMode_2LPT(  ScaleDep.z[ScaleDep.myseg], myobj->myk) - GrowingMode_2LPT(  ScaleDep.z[ScaleDep.myseg-1], myobj->myk) );
#ifdef THREE_LPT
      myobj->w31 = (GrowingMode_3LPT_1(                  myobj->z, myobj->myk) - GrowingMode_3LPT_1(ScaleDep.z[ScaleDep.myseg-1], myobj->myk)) /
                   (GrowingMode_3LPT_1(ScaleDep.z[ScaleDep.myseg], myobj->myk) - GrowingMode_3LPT_1(ScaleDep.z[ScaleDep.myseg-1], myobj->myk));
      myobj->w32 = (GrowingMode_3LPT_2(                  myobj->z, myobj->myk) - GrowingMode_3LPT_2(ScaleDep.z[ScaleDep.myseg-1], myobj->myk)) /
                   (GrowingMode_3LPT_2(ScaleDep.z[ScaleDep.myseg], myobj->myk) - GrowingMode_3LPT_2(ScaleDep.z[ScaleDep.myseg-1], myobj->myk) );   
#endif
#endif
    }

}

void set_obj_vel(int grp,PRODFLOAT F,pos_data *myobj)
{
  /* this sets the growth rates for peculiar velocities for a group */
  PRODFLOAT fac=Hubble(myobj->z)/(1.+myobj->z)*params.InterPartDist;

  myobj->Dv  = fac * fomega(myobj->z, myobj->myk);
#ifdef TWO_LPT
  myobj->D2v = fac * fomega_2LPT(myobj->z, myobj->myk);
#ifdef THREE_LPT
  myobj->D31v = fac * fomega_3LPT_1(myobj->z, myobj->myk);
  myobj->D32v = fac * fomega_3LPT_2(myobj->z, myobj->myk);
#endif
#endif

}



void set_point(int i,int j,int k,int ind,PRODFLOAT F,pos_data *myobj)
{
  /* sets all the quantities needed to handle a particle */

  myobj->z=F-1.0;

#ifdef SCALE_DEPENDENT
  int S=Smoothing.Nsmooth-1;
  myobj->myk=Smoothing.k_GM_displ[S];
#else
  myobj->myk=params.k_for_GM;
#endif

  set_weight(myobj);

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
  myobj->v_prev[0]=frag[ind].Vel_prev[0];
  myobj->v_prev[1]=frag[ind].Vel_prev[1];
  myobj->v_prev[2]=frag[ind].Vel_prev[2];
#ifdef TWO_LPT
  myobj->v2_prev[0]=frag[ind].Vel_2LPT_prev[0];
  myobj->v2_prev[1]=frag[ind].Vel_2LPT_prev[1];
  myobj->v2_prev[2]=frag[ind].Vel_2LPT_prev[2];
#ifdef THREE_LPT
  myobj->v31_prev[0]=frag[ind].Vel_3LPT_1_prev[0];
  myobj->v31_prev[1]=frag[ind].Vel_3LPT_1_prev[1];
  myobj->v31_prev[2]=frag[ind].Vel_3LPT_1_prev[2];
  myobj->v32_prev[0]=frag[ind].Vel_3LPT_2_prev[0];
  myobj->v32_prev[1]=frag[ind].Vel_3LPT_2_prev[1];
  myobj->v32_prev[2]=frag[ind].Vel_3LPT_2_prev[2];
#endif
#endif

#endif

}


void set_group(int grp,pos_data *myobj)
{
  /* copies relevant information from an object to a group data */

  groups[grp].Mass=myobj->M;
  for (int i=0;i<3;i++)
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
      groups[grp].Vel_prev[i]=myobj->v_prev[i];
#ifdef TWO_LPT
      groups[grp].Vel_2LPT_prev[i]=myobj->v2_prev[i];
#ifdef THREE_LPT
      groups[grp].Vel_3LPT_1_prev[i]=myobj->v31_prev[i];
      groups[grp].Vel_3LPT_2_prev[i]=myobj->v32_prev[i];
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

  if (!ScaleDep.myseg)
    {
      /* if the redshift falls before the end of the first segment
         (including the case when fragmentation is not segmented)
         then there is no need to interpolate */
      pos = myobj->q[i] + myobj->w * myobj->v[i];
#ifdef TWO_LPT
      if (order>1)
        pos += myobj->w2 * myobj->v2[i];
#ifdef THREE_LPT
      if (order>2)
        pos += myobj->w31 * myobj->v31[i]
            +  myobj->w32 * myobj->v32[i];
#endif
#endif

    }
#ifdef RECOMPUTE_DISPLACEMENTS
  else
    {
      /* else interpolate the velocities among two segments */
      pos = myobj->q[i] + (1.-myobj->w)*myobj->v_prev[i] + myobj->w*myobj->v[i];
#ifdef TWO_LPT
      if (order>1)
        pos += (1.-myobj->w2)*myobj->v2_prev[i] + myobj->w2*myobj->v2[i];
#ifdef THREE_LPT
      if (order>2) 
        pos += (1.-myobj->w31)*myobj->v31_prev[i] + myobj->w31*myobj->v32[i]
	    +  (1.-myobj->w32)*myobj->v32_prev[i] + myobj->w32*myobj->v32[i];
#endif
#endif
    }
#endif

  if (pbc)
    {
      /* impose periodic boundary conditions if required */
      if (pos>=Box) pos-=Box;
      if (pos< 0.0) pos+=Box;
    }

  return pos;
}



PRODFLOAT vel(int i, pos_data *myobj)
{
  /* this gives the velocity of a group in km/s */
  PRODFLOAT vv;

  if (!ScaleDep.myseg)
    {
      /* if the redshift falls before the end of the first segment
         (including the case when fragmentation is not segmented)
         then there is no need to interpolate */
      vv = myobj->v[i] * myobj->Dv * myobj->w;
#ifdef TWO_LPT
      vv += myobj->v2[i] * myobj->D2v  * myobj->w2;
#ifdef THREE_LPT
      vv += myobj->v31[i] * myobj->D31v * myobj->w31
         +  myobj->v32[i] * myobj->D32v * myobj->w32;
#endif
#endif

    }
#ifdef RECOMPUTE_DISPLACEMENTS
  else
    {
      /* else interpolate the velocities among two segments */
      vv = (myobj->v_prev[i]*(1.-myobj->w) + myobj->v[i]*myobj->w) * myobj->Dv;
#ifdef TWO_LPT
      vv += (myobj->v2_prev[i]*(1.-myobj->w2) + myobj->v2[i]*myobj->w2) * myobj->D2v;
#ifdef THREE_LPT
      vv += (myobj->v31_prev[i]*(1.-myobj->w31) + myobj->v31[i]*myobj->w31) * myobj->D31v
         +  (myobj->v32_prev[i]*(1.-myobj->w32) + myobj->v32[i]*myobj->w32) * myobj->D32v;
#endif
#endif
    }
#endif

  return vv;
}


PRODFLOAT distance(int i,pos_data *obj1,pos_data *obj2)
{
  /* return the 1D distance of two objects, respecting PBCs */

  PRODFLOAT d;

  /* here displacements are computed at the ORDER_FOR_GROUPS order */
  d = q2x(i,obj2,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_GROUPS)
    - q2x(i,obj1,subbox.pbc[i],(double)subbox.Lgwbl[i],ORDER_FOR_GROUPS);

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
  /* Updates a group after a merging or accretion */

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

      obj1->v_prev[i] = (obj1->v_prev[i]*obj1->M + obj2->v_prev[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef TWO_LPT
      obj1->v2_prev[i] = (obj1->v2_prev[i]*obj1->M + obj2->v2_prev[i]*obj2->M)/(double)(obj1->M + obj2->M);
#ifdef THREE_LPT
      obj1->v31_prev[i] = (obj1->v31_prev[i]*obj1->M + obj2->v31_prev[i]*obj2->M)/(double)(obj1->M + obj2->M);
      obj1->v32_prev[i] = (obj1->v32_prev[i]*obj1->M + obj2->v32_prev[i]*obj2->M)/(double)(obj1->M + obj2->M);
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
  /* Condition for a group being inside or outside the PLC 
     It computes the difference between the comoving distance of the
     group from the observer (given the replication) and the comoving
     distance at redshift z=F-1 */

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
  /* store halos in the PLC structure, ready to be written on file */

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
  /* brent root finder for a numerical equation */

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
      /*        return r; */
      
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

                  if (!get_map_bit_coord(i,j,k))
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

