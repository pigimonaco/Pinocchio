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

double inverse_collapse_time(double *, double *, double *, double *, int *);
double ell(double,double,double,double,double);
void ord(double *,double *,double *);

int compute_collapse_times(int ismooth)
{

  int local_x,local_y,local_z,index,i,fail;
  double Fnew,diff_ten[6],delta,local_average,local_variance,global_average,global_variance,Np;
  double lambda1,lambda2,lambda3;

  /* initialize the variance of the slice */
  local_variance = 0.0;
  local_average  = 0.0;

  /* initialize Fmax if it is the first smoothing */
  if (!ismooth)
    for (i=0; i<MyGrids[0].total_local_size; i++)
      {
	products[i].Fmax   = -10.0;
	products[i].Rmax   = -1;
	products[i].Vmax[0] = 0.0;
	products[i].Vmax[1] = 0.0;
	products[i].Vmax[2] = 0.0;
#ifdef TWO_LPT
	products[i].Vmax_2LPT[0] = 0.0;
	products[i].Vmax_2LPT[1] = 0.0;
	products[i].Vmax_2LPT[2] = 0.0;
#ifdef THREE_LPT
	products[i].Vmax_3LPT_1[0] = 0.0;
	products[i].Vmax_3LPT_1[1] = 0.0;
	products[i].Vmax_3LPT_1[2] = 0.0;
	products[i].Vmax_3LPT_2[0] = 0.0;
	products[i].Vmax_3LPT_2[1] = 0.0;
	products[i].Vmax_3LPT_2[2] = 0.0;
#endif
#endif
      }

  /* loop on all particles in the slice */
  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
    {
      for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
	{
	  for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
	    {

	      index = local_x + (MyGrids[0].GSlocal_x) * 
		(local_y + local_z * MyGrids[0].GSlocal_y);

	      /* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */

	      /* sum all contributions */
	      for (i=0; i<6; i++)
		diff_ten[i] = second_derivatives[0][i][index];

	      /* Computation of the variance of the linear density field */
	      delta = diff_ten[0]+diff_ten[1]+diff_ten[2];
	      local_average  = local_average  + delta;
	      local_variance = local_variance + delta*delta;

	      /* Computation of ellipsoidal collapse */
	      Fnew = inverse_collapse_time(diff_ten, &lambda1, &lambda2, &lambda3, &fail);

#ifdef SCALE_DEPENDENT_GROWTH
	      SDGM.flag=0;
	      SDGM.ismooth=ismooth;
#endif
	      /* comment these two lines out if F is the inverse growing mode */
	      if (Fnew>0)
		Fnew = 1.+InverseGrowingMode(1./Fnew);

	      if (fail)
		{
		  printf("ERROR on task %d: failure in inverse_collapse_time\n",ThisTask);
		  fflush(stdout);
		  return 1;
		}

	      if (products[index].Fmax < Fnew)
		{
		  products[index].Fmax = Fnew;
		  products[index].Rmax = ismooth;
		}

	    }
	}
    }

  /* Computes the obtained variance */
  global_variance = 0.0;
  global_average  = 0.0;

  MPI_Reduce(&local_variance, &global_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&local_average , &global_average , 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!ThisTask)
    {
      Np = (double)MyGrids[0].GSglobal_x * (double)MyGrids[0].GSglobal_y *(double)MyGrids[0].GSglobal_z;
      global_variance /= Np;
      global_average  /= Np;
    }

  MPI_Bcast(&global_variance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  Smoothing.TrueVariance[ismooth]=global_variance;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=0;
  SDGM.ismooth=ismooth;
  Smoothing.TrueVariance[ismooth] *= GrowingMode(params.camb.ReferenceRedshift);
#endif


  return 0;
}



double inverse_collapse_time(double *deformation_tensor, double *x1, double *x2, double *x3, int *fail)
{

  int m;
  double mu1,mu2,mu3,c,q,r,sq,t;

  /* mu1, mu2 and mu3 are the principal invariants of the 3x3 tensor 
     of second derivatives. */

  *fail=0;

  /* diagonalization of the tensor */
  mu1=deformation_tensor[0]+deformation_tensor[1]+deformation_tensor[2];
  mu2=.5*mu1*mu1;
  for (m=0; m<6; m++)
    {
      if (m<3)
	c=0.5;
      else
	c=1.0;
      mu2 -= c*deformation_tensor[m]*deformation_tensor[m];
    }

  mu3 =  deformation_tensor[0]*deformation_tensor[1]*deformation_tensor[2]
    + 2.*deformation_tensor[3]*deformation_tensor[4]*deformation_tensor[5]
       - deformation_tensor[0]*deformation_tensor[5]*deformation_tensor[5]
       - deformation_tensor[1]*deformation_tensor[4]*deformation_tensor[4]
       - deformation_tensor[2]*deformation_tensor[3]*deformation_tensor[3];

  q=(mu1*mu1-3.0*mu2)/9.0;

  if (q==0.)
    {
      /* in this case the tensor is already diagonal */
      *x1=deformation_tensor[0];
      *x2=deformation_tensor[1];
      *x3=deformation_tensor[2];
    }
  else
    {
      r=-(2.*mu1*mu1*mu1 - 9.0*mu1*mu2 + 27.0*mu3)/54.;

      if (q*q*q < r*r || q<0.0)
	{
	  *fail=1;
	  return -10.0;
	}

      sq=sqrt(q);
      t=acos(r/q/sq);
      *x1=-2.*sq*cos(t/3.)+mu1/3.;
      *x2=-2.*sq*cos((t+2.*PI)/3.)+mu1/3.;
      *x3=-2.*sq*cos((t+4.*PI)/3.)+mu1/3.;
    }

  /* ordering and inverse collapse time */
  /* ord(&x1,&x2,&x3); */
  /* t=ell(x1,x2,x3,mu1,mu3); */
  ord(x1,x2,x3);
  t=ell(*x1,*x2,*x3,mu1,mu3);

  /* return value */
  if (t>0.)
    return 1./t;
  else
    return -10.;
}



double ell(double l1,double l2,double l3,double del,double det)
{

  /*
    This rouine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
  */


  double a1,a2,a3,den,dis,r,q,sq,t,s1,s2,s3,ell;

  den=det/126.+5.*l1*del*(del-l1)/84.;
  if (l1==0.)
    ell=-0.1;

  else if (den==0.)
    {
     if (del==l1)
       {
        if (l1>0.0)
	  ell=1./l1;
        else
	  ell=-.1;
       }
     else
       {
	 dis=7.*l1*(l1+6.*del);
	 if (dis<0.0)
           ell=-.1;
	 else
	   {
	     ell=(7.*l1-sqrt(dis))/(3.*l1*(l1-del));
	     if (ell<0.) 
	       ell=-.1;
	   }
       }
    }
  else
    {
      a1=3.*l1*(del-l1)/14./den;
      a2=l1/den;
      a3=-1./den;
      q=(a1*a1-3.*a2)/9.;
      r=(2.*a1*a1*a1-9.*a1*a2+27.*a3)/54.;
      if (r*r>q*q*q)
	{
	  sq=pow(sqrt(r*r-q*q*q)+fabs(r),0.333333333333333);
	  ell=-fabs(r)/r*(sq+q/sq)-a1/3.;
	  if (ell<0.) 
	    ell=-.1;
	}
      else
	{
	  sq=sqrt(q);
	  t=acos(r/q/sq);
	  s1=-2.*sq*cos(t/3.)-a1/3.;
	  s2=-2.*sq*cos((t+2.*PI)/3.)-a1/3.;
	  s3=-2.*sq*cos((t+4.*PI)/3.)-a1/3.;
	  if (s1<0.) 
	    s1=1.e10;
	  if (s2<0.) 
	    s2=1.e10;
	  if (s3<0.) 
	    s3=1.e10;
	  ell=(s1<s2?s1:s2);
	  ell=(s3<ell?s3:ell);
	  if (ell==1.e10) 
	    ell=-.1;
	}
    }
  if (del > 0. && ell > 0.)
    ell+=-.364/del*exp(-6.5*(l1-l2)/del-2.8*(l2-l3)/del);

  return ell;
}


#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

void ord(double *a, double *b, double *c)
{
  /* orders a,b,c in decreasing order a>b>c */
  double lo,hi;

  hi=max( max(*a,*b), *c );
  lo=min( min(*a,*b), *c );

  *b=*a+*b+*c-lo-hi;
  *a=hi;
  *c=lo;
}


int compute_velocities(int ismooth)
{

  int local_x, local_y, local_z, index, i;

#ifdef SMOOTH_VELOCITIES
  int final=(ismooth==Smoothing.Nsmooth-1);
#endif
  /* loop on all particles in the slice */
  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
	{

	  index = local_x + (MyGrids[0].GSlocal_x) * 
	    (local_y + local_z * MyGrids[0].GSlocal_y);

	  /* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */

#ifdef SMOOTH_VELOCITIES
	  if (products[index].Rmax==ismooth || (products[index].Fmax<0.0 && final) )
#endif
	    for (i=0;i<3;i++)
	      products[index].Vmax[i]=first_derivatives[0][i][index];

#ifdef TIMELESS_SNAPSHOT
	  /* sets accretion redshift to -1 default value */
	  products[index].zacc=-1;
#endif
	}

  return 0;
}
