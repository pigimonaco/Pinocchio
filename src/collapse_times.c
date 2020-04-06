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

/* TODO
   - write table in binary format
   - let the code read the table if required
*/

#include "pinocchio.h"
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

//#define ELL_SNG
//#define ELL_CLASSIC

double inverse_collapse_time(int, double *, double *, double *, double *, int *);
double ell(int,double,double,double);
double ell_classic(int,double,double,double);
double ell_sng(int,double,double,double);
void ord(double *,double *,double *);
#ifdef TABULATED_CT
double interpolate_collapse_time(int ismooth, double l1, double l2, double l3);

// quanto segue alla fine potra` essere spostato giu` dove inizia 

//#define DEBUG
//#define PRINTJUNK

//#define ALL_SPLINE
//#define TRILINEAR
//#define HISTO
#define BILINEAR_SPLINE

#define CT_NBINS_XY (50)  // 50
#define CT_NBINS_D  (100) // 100
#define CT_SQUEEZE (1.2)  // 1.2
#define CT_EXPO   (1.75)  // 1.75
#define CT_RANGE_D (7.0)  // 7.0
#define CT_RANGE_X (3.5)
#define CT_DELTA0 (-1.0)
double *CT_table=0x0;
double bin_x;
static gsl_interp_accel *accel=0x0;
gsl_spline ***CT_Spline;
double *delta_vector;
int Ncomputations, start, length;

#ifdef TRILINEAR
#ifdef HISTO
int *CT_histo=0x0;
#endif
#endif

#ifdef DEBUG
double ave_diff=0.0,var_diff=0.0;
int counter=0;
#ifdef PRINTJUNK
FILE *JUNK;
#endif
#endif

#endif

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
	products[i].Vel[0] = 0.0;
	products[i].Vel[1] = 0.0;
	products[i].Vel[2] = 0.0;
#ifdef TWO_LPT
	products[i].Vel_2LPT[0] = 0.0;
	products[i].Vel_2LPT[1] = 0.0;
	products[i].Vel_2LPT[2] = 0.0;
#ifdef THREE_LPT
	products[i].Vel_3LPT_1[0] = 0.0;
	products[i].Vel_3LPT_1[1] = 0.0;
	products[i].Vel_3LPT_1[2] = 0.0;
	products[i].Vel_3LPT_2[0] = 0.0;
	products[i].Vel_3LPT_2[1] = 0.0;
	products[i].Vel_3LPT_2[2] = 0.0;
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
	      Fnew = inverse_collapse_time(ismooth, diff_ten, &lambda1, &lambda2, &lambda3, &fail);
		
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

  return 0;
}



double inverse_collapse_time(int ismooth, double *deformation_tensor, double *x1, double *x2, double *x3, int *fail)
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
  ord(x1,x2,x3);

#ifdef TABULATED_CT
  t=interpolate_collapse_time(ismooth,*x1,*x2,*x3);

#ifdef DEBUG
  double t2=ell(ismooth,*x1,*x2,*x3);
  if (t2>0.9)
    {
      ave_diff+=t-t2;
      var_diff+=pow(t-t2,2.0);
      counter+=1;
#ifdef PRINTJUNK
      double ampl = sqrt(Smoothing.Variance[ismooth]);
      fprintf(JUNK,"%d  %f %f %f   %f %f %f   %20g %20g %20g\n",ismooth, (*x1+*x2+*x3)/ampl, (*x1-*x2)/ampl, (*x2-*x3)/ampl, *x1, *x2, *x3, t, t2, t-t2);
#endif
    }

#endif

#else
  t=ell(ismooth,*x1,*x2,*x3);
#endif  

    return t;
}


#define SMALL ((double)1.e-20)

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

  /* loop on all particles in the slice */
  for (local_z=0; local_z < MyGrids[0].GSlocal_z; local_z++)
    for (local_y=0; local_y < MyGrids[0].GSlocal_y; local_y++)
      for (local_x=0; local_x < MyGrids[0].GSlocal_x; local_x++)
	{

	  index = local_x + (MyGrids[0].GSlocal_x) * 
	    (local_y + local_z * MyGrids[0].GSlocal_y);

	  /* NB: THIS IS THE POINT WHERE DIFFERENT BOXES SHOULD BE CONSIDERED */

	  for (i=0;i<3;i++)
	    products[index].Vel[i]=first_derivatives[0][i][index];
	}

  return 0;
}

#ifdef TABULATED_CT

int initialize_collapse_times(int ismooth)
{
  double l1, l2, l3, del, ampl, x, y, interval, ref_interval, deltaf;
  int id,ix,iy,i,j;

  /* the bin size for delta varies like (delta - CT_DELTA0)^CT_EXPO, 
     so it is smallest around CT_DELTA0 if CT_EXPO>1. The bin size 
     is never smaller than CT_MIN_INTERVAL.
     The delta vector is computed at the first smoothing radius. */

  if (!ismooth)
    {
      delta_vector=(double*)malloc(CT_NBINS_D * sizeof(double));

      /* sets the sampling in density delta */
      if (CT_EXPO==1)
	{
	  /* in this case the spacing is even */
	  interval=2.*CT_RANGE_D/(double)CT_NBINS_D;
	  for (id=0; id<CT_NBINS_D; id++)
	    delta_vector[id]=id*interval -CT_RANGE_D;
	}
      else
	{
	  /* sampling is finer around CT_DELTA0 */
	  deltaf = pow(CT_SQUEEZE/CT_EXPO,1./(CT_EXPO-1.));
	  if (CT_EXPO==2)
	    ref_interval=( ( log((CT_RANGE_D-CT_DELTA0)/deltaf) + log((CT_RANGE_D+CT_DELTA0)/deltaf) )/ 
			   CT_EXPO + 2.*deltaf/CT_SQUEEZE ) / (CT_NBINS_D-2.0);
	  else
	    ref_interval=( ( pow(CT_RANGE_D-CT_DELTA0,2.-CT_EXPO) + pow(CT_RANGE_D+CT_DELTA0,2.-CT_EXPO) 
			     - 2.*pow(deltaf,2.-CT_EXPO) )/ 
			   CT_EXPO/(2.-CT_EXPO) + 2.*deltaf/CT_SQUEEZE ) / (CT_NBINS_D-2.0);
	  id=0;
	  del=-CT_RANGE_D;
	  do 
	    {
	      delta_vector[id]=del;
	      interval = CT_EXPO * ref_interval * pow(fabs(del-CT_DELTA0),CT_EXPO-1.0);
	      interval = (interval/ref_interval<CT_SQUEEZE? ref_interval*CT_SQUEEZE : interval);
	      del+=interval;
	      id++;
	    }  
	  while (id<CT_NBINS_D);
	}

      if (!ThisTask)
	printf("[%s] Grid for interpolating collapse times: CT_NBINS_D=%d, CT_NBINS_XY=%d\n",fdate(),CT_NBINS_D,CT_NBINS_XY);

      /* computations are divided among tasks */
      Ncomputations = CT_NBINS_D * CT_NBINS_XY * CT_NBINS_XY;
      bin_x = CT_RANGE_X /(double)(CT_NBINS_XY);

      /* allocates spline for interpolations */
      accel = gsl_interp_accel_alloc();
      CT_Spline=(gsl_spline***)calloc(CT_NBINS_XY,sizeof(gsl_spline**));
      for (i=0; i<CT_NBINS_XY; i++)
	{
	  CT_Spline[i]=(gsl_spline**)calloc(CT_NBINS_XY,sizeof(gsl_spline*));
	  for (j=0; j<CT_NBINS_XY; j++)
	    CT_Spline[i][j]=gsl_spline_alloc(gsl_interp_cspline,CT_NBINS_D);
	}
      CT_table=(double *)calloc(Ncomputations, sizeof(double));
#ifdef TRILINEAR
#ifdef HISTO
      CT_histo=(int *)calloc(Ncomputations, sizeof(int));
#endif
#endif

#ifdef DEBUG
#ifdef PRINTJUNK
      char filename[50];
      sprintf(filename,"debug.task%d.txt",ThisTask);
      JUNK=fopen(filename,"w");
#endif
#endif

    }

  /* The collapse time table is computed for each smoothing radius */
  //  if (!ismooth)
    {

      ampl = sqrt(Smoothing.Variance[ismooth]);

      /* map of computations of all tasks */
      int *counts=(int*)malloc(NTasks*sizeof(int));
      int *displs=(int*)malloc(NTasks*sizeof(int));
      int bunch = Ncomputations / NTasks;
      int remainder = Ncomputations % NTasks;

      for (i=0, ix=0; i<NTasks; i++)
	{
	  counts[i]=bunch + (i < remainder);
	  displs[i]=bunch*i + ((remainder > i)? i : remainder );
	}

      /* computations for this task */
      for (i=displs[ThisTask]; i<displs[ThisTask]+counts[ThisTask]; i++)
	{
	  id=i%CT_NBINS_D;
	  ix=(i/CT_NBINS_D)%CT_NBINS_XY;
	  iy=i/CT_NBINS_D/CT_NBINS_XY;

	  x   = ix * bin_x;
	  y   = iy * bin_x;
	  l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
	  l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
	  l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

	  CT_table[i]= ell(ismooth,l1,l2,l3);
	}

      MPI_Allgatherv( MPI_IN_PLACE, counts[ThisTask], MPI_DOUBLE, CT_table, counts, displs, MPI_DOUBLE, MPI_COMM_WORLD );

      free(displs);
      free(counts);

      /* now initialize the splines */
      for (i=0; i<CT_NBINS_XY; i++)
	for (j=0; j<CT_NBINS_XY; j++)
	  {
	    gsl_spline_init(CT_Spline[i][j],delta_vector,&CT_table[i*CT_NBINS_D+j*CT_NBINS_D*CT_NBINS_XY],CT_NBINS_D);
	  }


      /* Writes the table on a file */
      if (!ThisTask)
	{
	  char fname[LBLENGTH];
	  FILE *fd;
	  sprintf(fname,"pinocchio.%s.CTtable.out",params.RunFlag);
	  if (!ismooth)
	    fd=fopen(fname,"w");
	  else
	    fd=fopen(fname,"a");
	  fprintf(fd,"################\n# smoothing %2d #\n################\n",ismooth);
	  for (i=0; i<Ncomputations; i++)
	    {
	      id=i%CT_NBINS_D;
	      ix=(i/CT_NBINS_D)%CT_NBINS_XY;
	      iy=i/CT_NBINS_D/CT_NBINS_XY;

	      x   = ix * bin_x;
	      y   = iy * bin_x;
	      l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
	      l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
	      l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

	      fprintf(fd," %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20g\n",
		      ismooth,id,ix,iy,delta_vector[id],x,y,l1,l2,l3,CT_table[i]);
	    }
	  fclose(fd);
	}
    }

  return 0;
}


int reset_collapse_times(int ismooth)
{
  

#ifdef DEBUG
  if (counter)
    {
#ifdef TRILINEAR
      int flag=1;
#endif
#ifdef BILINEAR_SPLINE
      int flag=2;
#endif
#ifdef ALL_SPLINE
      int flag=3;
#endif

      double aa=ave_diff/(double)counter;
      printf("%d   %2d  %10d   %3d %3d %5.2f %5.2f   %20g  %20g DEBUG\n",
	     flag, ismooth, counter, CT_NBINS_D, CT_NBINS_XY, CT_EXPO, CT_SQUEEZE,
	     aa,sqrt(var_diff/(double)counter)-aa*aa);
    }
  /* else */
  /*   printf("N=0, no average or std dev\n"); */
  ave_diff=var_diff=0.0;
  counter=0;

#endif

#ifdef TRILINEAR
#ifdef HISTO
  // scrive l'istogramma della tabella
  char fname[50];
  FILE *fd;
  int i,id,ix,iy;
  double x,y,l1,l2,l3;
  sprintf(fname,"popolazione.txt");
  if (!ismooth)
    fd=fopen(fname,"w");
  else
    fd=fopen(fname,"a");
  fprintf(fd,"################\n# smoothing %2d #\n################\n",ismooth);
  double ampl = sqrt(Smoothing.Variance[ismooth]);
  for (i=0; i<Ncomputations; i++)
    {
      id=i%CT_NBINS_D;
      ix=(i/CT_NBINS_D)%CT_NBINS_XY;
      iy=i/CT_NBINS_D/CT_NBINS_XY;

      x   = ix * bin_x;
      y   = iy * bin_x;
      l1  = (delta_vector[id] + 2.*x +    y)/3.0*ampl;
      l2  = (delta_vector[id] -    x +    y)/3.0*ampl;
      l3  = (delta_vector[id] -    x - 2.*y)/3.0*ampl;

      fprintf(fd," %d   %3d %3d %3d   %8f %8f %8f   %8f %8f %8f   %20d\n",
	      ismooth,id,ix,iy,delta_vector[id],x,y,l1,l2,l3,CT_histo[i]);
    }
  fclose(fd);
#endif
#endif

  return 0;
}


int compare_search(const void *A, const void *B)
{
  if ( *(double*)A >= *(double*)B && *(double*)A < *((double*)B+1))
    return 0;
  if ( *(double*)A < *(double*)B )
    return -1;
  return 1;
}



double interpolate_collapse_time(int ismooth, double l1, double l2, double l3)
{
  int ix, iy;
  double ampl, d, x, y;

  ampl = sqrt(Smoothing.Variance[ismooth]);
  d = (l1+l2+l3)/ampl;
  x = (l1-l2)/ampl;
  y = (l2-l3)/ampl;
  ix = (int)(x/bin_x);
  iy = (int)(y/bin_x);

  /* in the unlikely event the point is beyond the limits */
  if (ix>=CT_NBINS_XY-1)
    ix=CT_NBINS_XY-2;
  if (ix<0)
    ix=0;
  if (iy>=CT_NBINS_XY-1)
    iy=CT_NBINS_XY-2;
  if (iy<0)
    iy=0;

#ifdef ALL_SPLINE

  int ixx, iyy, ixstart, iystart, N=4;
  double xls[4],yls[4],zls[16];
  const gsl_interp2d_type *T = gsl_interp2d_bicubic;
  gsl_interp_accel *xacc = gsl_interp_accel_alloc();
  gsl_interp_accel *yacc = gsl_interp_accel_alloc();

  if (ix==0)
    ixstart=0;
  else if (ix>=CT_NBINS_XY-2)
    ixstart=CT_NBINS_XY-4;
  else
    ixstart=ix-1;

  if (iy==0)
    iystart=0;
  else if (iy>=CT_NBINS_XY-2)
    iystart=CT_NBINS_XY-4;
  else
    iystart=iy-1;

  /* little 4x4 array for interpolation */
  for (ixx=0; ixx<4; ixx++)
    {
      xls[ixx]=(ixx+ixstart)*bin_x;
      yls[ixx]=(ixx+iystart)*bin_x;
    }
  for (ixx=0; ixx<4; ixx++)
    for (iyy=0; iyy<4; iyy++)
      zls[ixx+iyy*4]=my_spline_eval(CT_Spline[ixx+ixstart][iyy+iystart], d, accel);

  /* initialize interpolation */
  gsl_spline2d *LittleSpline = gsl_spline2d_alloc(T, N, N);
  gsl_spline2d_init(LittleSpline, xls, yls, zls, N, N);

  double a=gsl_spline2d_eval(LittleSpline, x, y, xacc, yacc);
  gsl_spline2d_free(LittleSpline);

  return a;

#endif

#ifdef TRILINEAR

  int id;
  if (d<=delta_vector[0])
    id=0;
  else if (d>=delta_vector[CT_NBINS_D-1])
    id = CT_NBINS_D-2;
  else
    id = (double*)bsearch((void*)&d, (void*)delta_vector, (size_t)(CT_NBINS_D-1), sizeof(double), compare_search)-delta_vector;

#ifdef HISTO
  CT_histo[id + ix*CT_NBINS_D + iy*CT_NBINS_D*CT_NBINS_XY]++;
#endif

  double dd = (d-delta_vector[id])/(delta_vector[id+1]-delta_vector[id]);
  double dx = x/bin_x-ix;
  double dy = y/bin_x-iy;

  return ((1.-dd)*(1.-dx)*(1.-dy)*CT_table[id   + (ix  )*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(1.-dx)*(1.-dy)*CT_table[id+1 + (ix  )*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(   dx)*(1.-dy)*CT_table[id   + (ix+1)*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(   dx)*(1.-dy)*CT_table[id+1 + (ix+1)*CT_NBINS_D + (iy  )*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(1.-dx)*(   dy)*CT_table[id   + (ix  )*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(1.-dx)*(   dy)*CT_table[id+1 + (ix  )*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (1.-dd)*(   dx)*(   dy)*CT_table[id   + (ix+1)*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]+
	  (   dd)*(   dx)*(   dy)*CT_table[id+1 + (ix+1)*CT_NBINS_D + (iy+1)*CT_NBINS_D*CT_NBINS_XY]);

#endif

#ifdef BILINEAR_SPLINE

  double dx = x/bin_x-ix;
  double dy = y/bin_x-iy;


  /* if (ismooth==10) */
  /*   printf(" %f %f %f    %f %f %f   %d %d   %f\n",l1,l2,l3,d,x,y,ix,iy,my_spline_eval(CT_Spline[ix  ][iy  ], d, accel)); */

  return ((1.-dx)*(1.-dy)*my_spline_eval(CT_Spline[ix  ][iy  ], d, accel)+
	  (   dx)*(1.-dy)*my_spline_eval(CT_Spline[ix+1][iy  ], d, accel)+
	  (1.-dx)*(   dy)*my_spline_eval(CT_Spline[ix  ][iy+1], d, accel)+
	  (   dx)*(   dy)*my_spline_eval(CT_Spline[ix+1][iy+1], d, accel));

#endif

}


#endif 


double ell(int ismooth, double l1, double l2, double l3)
{

#ifdef ELL_CLASSIC
  double bc=ell_classic(ismooth, l1, l2, l3);
  if (bc>0.0)
    return 1.+InverseGrowingMode(bc,ismooth);
  else
    return 0.0;
#endif

#ifdef ELL_SNG
  double bc=ell_sng(ismooth, l1, l2, l3);
  if (bc>0.0)
    return 1./bc;
  else
    return 0.0;
#endif

}


/* implementation of ellipsoidal collapse following Nadkarni-Ghosh & Singhal (2016) */

#ifdef MOD_GRAV_FR

double ForceModification(double size, double a, double delta)
{

  double ff = 4.*params.OmegaLambda/params.Omega0;
  double thickness = FR0 / params.Omega0 / pow(H_over_c * size,2.0)
  * pow(a,7.)  * pow((1.+delta),-1./3.)
  * ( pow((1.0 + ff) / (1.0 + ff*pow(a,3.)) ,2.0) - pow((1.0 + ff) / (1.0 + delta + ff*pow(a,3.)) ,2.0));
  double F3=(thickness * (3. + thickness * (-3. + thickness)));
  if (F3<0.)
    F3=0.;

  return (F3<1. ? F3/3. : 1./3);
}

#endif

/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */

int sng_system(double t, const double y[], double f[], void *sng_par)
{

  int i,j;
  double sum;
  double omegam=OmegaMatter(1./t-1.);
  double omegal=OmegaLambda(1./t-1.);
  double delta=y[6]+y[7]+y[8];

  for (i=0; i<3; i++) 
    {
      sum=0.;
      for (j=0; j<3; j++) 
	{
	  if (i==j || y[i] == y[j])
	    continue;
	  else
	    sum += (y[j+6] - y[i+6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i+3]) 
					- (1. - y[j]) * (1. - y[j]) * (1. + y[j+3])) 
	      / ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
	}
      f[i]   = (y[i+3]*(y[i]-1.0)) / t;
      f[i+3] = (0.5*(y[i+3]*(omegam - 2.0*omegal - 2.0) 
#ifdef MOD_GRAV_FR
		     - 3.0*omegam*y[i+6] * (1. + ForceModification(*(double*)sng_par, t, delta))
#else
		     - 3.0*omegam*y[i+6] 
#endif
		     - 2.0*y[i+3]*y[i+3])) / t;
      f[i+6] = ((5./6. + y[i+6]) * 
		((3. + y[3] + y[4] + y[5]) 
		 - (1. + delta) / (2.5 + delta) * (y[3] + y[4] + y[5]))
		- (2.5 + delta) * (1. + y[i+3]) + sum) / t;


      //printf("%f %f %f %f\n",t,*(double*)params,delta,ForceModification(*(double*)params, t, delta));
    }
  return GSL_SUCCESS;
}


double ell_sng(int ismooth, double l1, double l2, double l3)
{
  /* call of ODE solution for SNG ellipsoidal collapse */
  double ode_param, hh=1.e-6;
  double amin=1.e-5, amax=5.0;
  double mya = amin;

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9);
  gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);
  gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);
  gsl_odeiv2_system   ode_sys   = {sng_system, jac, 9, (void*)&ode_param};

  /* ICs */
#ifdef SCALE_DEPENDENT
  double D_in = GrowingMode(1./amin-1.,params.k_for_GM/Smoothing.Radius[ismooth]*params.InterPartDist);
#else
  double D_in = GrowingMode(1./amin-1.,0.);
#endif
  double y[9] = {l1*D_in, l2*D_in, l3*D_in,
		 l1*D_in/(l1*D_in - 1.), l2*D_in/(l2*D_in - 1.), l3*D_in/(l3*D_in - 1.),
		 l1*D_in, l2*D_in, l3*D_in};

#ifdef MOD_GRAV_FR
  if (ismooth<Smoothing.Nsmooth-1)
    ode_param=Smoothing.Radius[ismooth];
  else
    ode_param=Smoothing.Radius[ismooth-1];
#endif

  /* Loop to call ODE integrator */
  double olda=mya;
  double oldlam=y[0];

  while (mya<amax)
    {
      int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &mya, amax, &hh, y);
      if (status!=GSL_SUCCESS)
	{
	  printf("ERROR on task %d: integration of cosmological quantities failed\n",ThisTask);
	  fflush(stdout);
	  return -1;
	}
      if (y[0]>=0.99999)
	return olda + (1.-oldlam) * (mya - olda) / (y[0] - oldlam);
    }

  /* in this case the ellipsoid does not collapse */
  return 0;
}


double ell_classic(int ismooth, double l1, double l2, double l3)
{
  /*
    This rouine computes the smallest non-negative solution of the 3rd
    order equation for the ellipsoid, and corrects it to reproduce the
    spherical collapse correctly.
  */

  double a1,a2,a3,den,dis,r,q,sq,t,s1,s2,s3,ell,del,det;

  del=l1+l2+l3;
  det=l1*l2*l3;
  den=det/126.+5.*l1*del*(del-l1)/84.;
  if (fabs(l1)<SMALL)
    ell=-0.1;

  else if (fabs(den)<SMALL)
    {
     if (fabs(del-l1)<SMALL)
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



/* questa funziona, ma e` piu` lenta della spline delle gsl */

/* double cubic_int(double x, double *xv, double *yv, int N, int es) */
/* { */

/*   /\* returns a cubic interpolation of a function yv sampled on values xv */
/*      x: x-value at which interpolation is required */
/*      xv: x array, in increasing order */
/*      yv: y array */
/*      N: length of xv and yv */
/*      es: 1 if x-values in xv are evenly spaced *\/ */

/*   if (x<=xv[0]) */
/*     return yv[0]; */
/*   else if (x>=xv[N-1]) */
/*     return yv[N-1]; */

/*   int ix,ix1; */

/*   if (es) */
/*     ix = (int)( (x-xv[0])/(xv[1]-xv[0]) ); */
/*   else */
/*     ix = (double*)bsearch((void*)&x, (void*)xv, (size_t)N, sizeof(double), compare_search)-xv; */

/*   if (ix==0) */
/*     ix1=0; */
/*   else if (ix==N-2) */
/*     ix1=N-4; */
/*   else */
/*     ix1=ix-1; */

/*   int ix2=ix1+1, ix3=ix1+2, ix4=ix1+3; */

/*   return (  (x      -xv[ix2])*(x      -xv[ix3])*(x      -xv[ix4]) * yv[ix1] */
/* 	  / (xv[ix1]-xv[ix2])/(xv[ix1]-xv[ix3])/(xv[ix1]-xv[ix4]) + */
/* 	    (x      -xv[ix1])*(x      -xv[ix3])*(x      -xv[ix4]) * yv[ix2]  */
/* 	  / (xv[ix2]-xv[ix1])/(xv[ix2]-xv[ix3])/(xv[ix2]-xv[ix4]) + */
/* 	    (x      -xv[ix1])*(x      -xv[ix2])*(x      -xv[ix4]) * yv[ix3] */
/* 	  / (xv[ix3]-xv[ix1])/(xv[ix3]-xv[ix2])/(xv[ix3]-xv[ix4]) + */
/* 	    (x      -xv[ix1])*(x      -xv[ix2])*(x      -xv[ix3]) * yv[ix4] */
/* 	  / (xv[ix4]-xv[ix1])/(xv[ix4]-xv[ix2])/(xv[ix4]-xv[ix3]) ); */

/* } */

