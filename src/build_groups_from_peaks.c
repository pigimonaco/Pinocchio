#include "pinocchio.h"
#include <gsl/gsl_sf.h>

// TENTATIVO ABORTITO


PRODFLOAT nearest_encounter(pos_data, pos_data, double);
double LagDist(int, int, int);


int build_groups_from_peaks(int Npeaks, double zstop)
{

  pos_data obj1,obj2;

  // fuffa preliminare
  printf("Picchi (%d)\n",Npeaks);
  for (int i=FILAMENT+1; i<15; i++)
    printf(" i=%d, Fmax: %7.4f, q=[%7.4f,%7.4f,%7.4f], v=[%7.4f,%7.4f,%7.4f]\n",i,
	   groups[i].t_peak,groups[i].Pos[0],groups[i].Pos[1],groups[i].Pos[2],groups[i].Vel[0],groups[i].Vel[1],groups[i].Vel[2]);

  set_obj(2,1,&obj1);
  set_obj(3,1,&obj2);

  PRODFLOAT tmerge=nearest_encounter(obj1,obj2,(double)MyGrids[0].GSglobal[0]);
  printf("tmerge: %7.4f\n",tmerge);

  // da scalare con la dimensione Lagrangiana del piu` grande alone nella box
  double dmsq=pow(20./params.InterPartDist,2.);

  double *auxtimes=(double*)calloc(subbox.PredNpeaks,sizeof(double));  // dovremo capire dove mettere queste informazioni

  /* list of potential mergers */
  int nmerge=0;
  /* loop on peak pairs */
  for (int g1=FILAMENT+1; g1<FILAMENT+1+Npeaks; g1++)
    for (int g2=g1+1; g2<FILAMENT+1+Npeaks; g2++)
      {
	double dx=LagDist(0,g1,g2);
	double d2=dx*dx;
	if (d2<dmsq)
	  {
	    double dy=LagDist(1,g1,g2);
	    d2+=dy*dy;
	    if (d2<dmsq)
	      {
		double dz=LagDist(2,g1,g2);
		d2+=dz*dz;
		if (d2<=dmsq)
		  {
		    /* pairs beyond a threshold Lagrangian radius are excluded */
		    
		    set_obj(g1,1.0,&obj1);
		    set_obj(g2,1.0,&obj2);
		    
		    double zmerge=nearest_encounter(obj1,obj2,(double)MyGrids[0].GSglobal[0]);
		    
		    /* printf(" g1=%d, g2=%d, dq=[%7.4f,%7.4f,%7.4f], dv=[%7.4f,%7.4f,%7.4f], zmerge=%7.4f\n",g1,g2, */
		    /* 	   groups[g2].Pos[0]-groups[g1].Pos[0], */
		    /* 	   groups[g2].Pos[1]-groups[g1].Pos[1], */
		    /* 	   groups[g2].Pos[2]-groups[g1].Pos[2], */
		    /* 	   groups[g2].Vel[0]-groups[g1].Vel[0], */
		    /* 	   groups[g2].Vel[1]-groups[g1].Vel[1], */
		    /* 	   groups[g2].Vel[2]-groups[g1].Vel[2],zmerge); */

		    if (zmerge > 0)  // da controllare se non c'e` da andare un pelino piu` in la`
		      {
			group_ID[nmerge]=g1;      // uso due array che al momento non servono, poi sistemiamo
			linking_list[nmerge]=g2;
			auxtimes[nmerge]=zmerge;
			//printf(" %d  %d  %d  %f\n",nmerge, g1, g2, zmerge);
			nmerge++;
		      }
		  }
	      }
	  }

      }

  printf("number of potential mergings: %d\n",nmerge);

  return 0;
}



PRODFLOAT nearest_encounter(pos_data obj1, pos_data obj2, double Box)
{

  PRODFLOAT nearest;

  /* fix periodic boundary conditions for the pair */
  for (int ia=0; ia<3; ia++)
    if (fabs(obj1.q[ia]-obj2.q[ia]) > Box/2.)
      {
	if (obj1.q[ia] < obj2.q[ia])
	  obj1.q[ia] += Box;
	else
	  obj2.q[ia] += Box;
      }

  double dqdv = 
    (obj2.q[0]-obj1.q[0]) * (obj2.v[0]-obj1.v[0]) +
    (obj2.q[1]-obj1.q[1]) * (obj2.v[1]-obj1.v[1]) +
    (obj2.q[2]-obj1.q[2]) * (obj2.v[2]-obj1.v[2]);

  if (dqdv>0)
    return -1;

  double dvdv = 
    (obj2.v[0]-obj1.v[0]) * (obj2.v[0]-obj1.v[0]) +
    (obj2.v[1]-obj1.v[1]) * (obj2.v[1]-obj1.v[1]) +
    (obj2.v[2]-obj1.v[2]) * (obj2.v[2]-obj1.v[2]);

  double guess = - dqdv / dvdv;

  //printf("dqdv=%10.4f,  dvdv=%10.4f\n",dqdv,dvdv);

  nearest = guess;
  return nearest;
}


double LagDist(int ia, int g1, int g2) // da passare a PRODFLOAT
{

  /* computes the Lagrangian distance of two groups on dimension ia */
  double d = groups[g2].Pos[ia]-groups[g1].Pos[ia];
  if (d>MyGrids[0].GSglobal[ia]/2)
    d-=MyGrids[0].GSglobal[ia];
  else if (d<-MyGrids[0].GSglobal[ia]/2)
    d+=MyGrids[0].GSglobal[ia];

  return d;
}
