#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_spline.h>

#define CT_NBINS_XY (40)
#define CT_INTERVAL (0.05)
#define CT_MIN_INTERVAL (0.05)
//#define CT_EXPO (2.0)
#define CT_EXPO (1.0)
#define CT_RANGE (5.0)
#define CT_DELTA0 (-1.0)


double *delta_vector;
  int CT_NBINS_D,ic;


int compare_search(const void *A, const void *B) {
  int ans;

  if ((double*)B == delta_vector+CT_NBINS_D-1)
    ans = 1;
  if (*(double*)A >= *(double*)B && *(double*)A < *((double*)B+1))
    ans = 0;
  if (*(double*)A < *(double*)B )
    ans = -1;
  if (*(double*)A >= *((double*)B+1) )
    ans = 1;

  printf("A=%f, B=%f, answer: %d\n", *(double*)A, *(double*)B, ans);

  return ans;
}


int main(int argc, char **argv) {
  double d,del,interval,bin_delta;

  /* counts the bins in delta, that are not evenly spaced */
  if (CT_EXPO==1) {
	/* in this case the spacing is even */
	CT_NBINS_D=CT_NBINS_XY;
	delta_vector=(double*)malloc(CT_NBINS_D * sizeof(double));
	bin_delta=2.*CT_RANGE/(double)CT_NBINS_D;
	for (ic=0; ic<CT_NBINS_D; ic++) {
	  delta_vector[ic]=ic*bin_delta -CT_RANGE;
	  printf("%d %f\n",ic,delta_vector[ic]);
	}
  }
  else {
	/* sampling is finer around delta=CT_DELTA0 */
	CT_NBINS_D=0;
	bin_delta=0;
	del=-CT_RANGE;
	do  {
	  interval = CT_EXPO * CT_INTERVAL * pow(fabs(del-CT_DELTA0),CT_EXPO-1.0);
	  interval = (interval<CT_MIN_INTERVAL? CT_MIN_INTERVAL : interval);
	  del+=interval;
	  CT_NBINS_D++;
	}  
	while (del<CT_RANGE);
	
	
	/* this contains the delta values */
	delta_vector=(double*)malloc(CT_NBINS_D * sizeof(double));
	ic=0;
	del=-CT_RANGE;
	do {
	  delta_vector[ic]=del;
	  printf("%d %f\n",ic,del);
	  interval = CT_EXPO * CT_INTERVAL * pow(fabs(del-CT_DELTA0),CT_EXPO-1.0);
	  interval = (interval<CT_MIN_INTERVAL? CT_MIN_INTERVAL : interval);
	  del+=interval;
	  ic++;
	}  
	while (del<CT_RANGE);
  }
  
  
  d=atof(argv[1]);
  if (d<delta_vector[0] || d>=delta_vector[CT_NBINS_D-1]) {
	printf("out of range\n");
	return 1;
  }
  
  void *pp;
  pp = bsearch((void*)&d, (void*)delta_vector, (size_t)CT_NBINS_D, sizeof(double), compare_search);
  int id=(double*)pp-delta_vector;
  
  printf("%d [%f,%f], %f\n",id,delta_vector[id],delta_vector[id+1],d);
  
  return 0;
}
