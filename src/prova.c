
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define PI ((double)3.141592)
#define NYQUIST 1

int main(int argc, char **argv)
{
  /* if (argc<3) */
  /*   return 0; */
  /* int N=atoi(argv[1]); */
  /* int i0 = atoi(argv[2]); */
  /* int i; */

  /* for (i=i0; i<N;  */
  /*      i==i0? i=(i0==0) :  */
  /* 	 (i==i0-1? i=i0+1 : i++)) */
  /*   printf("%d\n",i); */


  double BoxSize=100.,ThisRadius=10.;
  int GSglobal_x=128,GSglobal_y=128,GSglobal_z=128;
  int i,j,k;
  double fundamental = 2. * PI / BoxSize;
  double nyquist = NYQUIST * GSglobal_x/2 * fundamental;
  double d3k=pow(fundamental,3.);
  double w, D, kmag;
  double Variance=0.0;

  for (i = GSglobal_x / 2; i <= GSglobal_x / 2; i++)
    {
      for (j = -GSglobal_y / 2; j <= GSglobal_y / 2; j++)
	{
	  for (k = -GSglobal_z / 2; k <= GSglobal_z / 2; k++)
	    {

	      kmag=sqrt((double)(i*i+j*j+k*k))*fundamental;
	      if(kmag > nyquist)
		continue;

	      w = exp(-k*k*ThisRadius*ThisRadius/2.);
	      //D = GrowingMode(1./Time-1.,kmag);
	      D=1.;
	      Variance += 1.e7*pow(kmag,-1.) * D*D * w*w * d3k;

	    }
	}
    }

  return 0;
}
