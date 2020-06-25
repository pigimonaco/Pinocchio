 /* TO DO: add interpolation to find time at which the first ld = 1 (collapse)
	-> this can be done by inverting l_d(a) -> a(l_d) and interpolate to find
	a(1), add read cosmo params from file, l1,l2,l3 should somehow be drawn
	from a Doroshkevich distribution */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_ellint.h>

#define littleh 0.67
#define OM0 0.3
#define OL0 0.7
#define OR0 4.2e-5/0.67/0.67
#define OK0 1.0 - OM0 - OL0 - OR0
#define H_0 67.

/* Cosmological functions: H(a), \Omega_m(a), \Omega_L(a) */
double E(double a) {
  return sqrt(OL0 + OK0/a/a + OM0/a/a/a + OR0/a/a/a/a);
}

double omegam(double a) {
  return OM0 / (E(a)*E(a)*a*a*a);
}

double omegal(double a) {
  return OL0 / (E(a)*E(a));
}


double delta_lin(double x, double y, double z) {
  return x+y+z;
}

int jac( double t, const double y[], double *dfdy, double dfdt[], void *params) {
  printf("This integration method should not call this function.\n");
  return GSL_FAILURE;
}

double sum(const double y[], int i) {
  int j;
  double tot=0.;
  for (j=0; j<3; j++) {
	if (j==i)
	  continue;
	if (y[i] == y[j])
	  tot += 0.;
	else {
	  tot += (y[j+6] - y[i+6]) * ((1. - y[i]) * (1. - y[i]) * (1. + y[i+3]) 
								  - (1. - y[j]) * (1. - y[j]) * (1. + y[j+3])) 
		/ ((1. - y[i]) * (1. - y[i]) - (1. - y[j]) * (1. - y[j]));
	}
  }
  return tot;
}

/* System of ODEs: specify f_i(t) = r.h.s. of differential equations */
/* f[i] = d(l_a)/da; f[i+3] = d(l_v)/da; d(l_d)/da                   */
int func(double t, const double y[], double f[], void *params) {
  (void)(t);
  int i;
  for (i=0; i<3; i++) {
	f[i]   = (y[i+3]*(y[i]-1.0)) / t;        
	f[i+3] = (0.5*(y[i+3]*(omegam(t) - 2.0*omegal(t) - 2.0) 
				   - 3.0*omegam(t)*y[i+6] - 2.0*y[i+3]*y[i+3])) / t;
	f[i+6] = ((5./6. + y[i+6]) * ((3. + y[3] + y[4] + y[5]) 
	  - (1. + delta_lin(y[6],y[7],y[8])) / (2.5 + delta_lin(y[6],y[7],y[8]))
								 * (y[3] + y[4] + y[5]))
	  - (2.5 + delta_lin(y[6],y[7],y[8])) * (1. + y[i+3]) + sum(y, i)) / t;
  }
  return GSL_SUCCESS;
}

double ellcoll(double l1, double l2, double l3) {
  int i;
  double ode_param, hh=1.e-6;
  double amin=1.e-3, amax=5., da=(amax - amin) / 250.;
  double scalef = amin;

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step    *ode_s     = gsl_odeiv2_step_alloc(T,9);
  gsl_odeiv2_control *ode_c     = gsl_odeiv2_control_standard_new(1.0e-6, 1.0e-6, 1.0, 1.0);
  gsl_odeiv2_evolve  *ode_e     = gsl_odeiv2_evolve_alloc(9);
  gsl_odeiv2_system   ode_sys   = {func, jac, 9, (void*)&ode_param};

  /* ICs */
  double a_in = amin;
  double y[9] = {l1*a_in, l2*a_in, l3*a_in,
				 l1*a_in/(l1*a_in - 1.), l2*a_in/(l2*a_in - 1.), l3*a_in/(l3*a_in - 1.),
				 l1*a_in, l2*a_in, l3*a_in};

  /* Loop to call ODE integrator */
  printf("Starting integration\n");
  for (i=0; i<=250; i++) {
	double ai = amin + i * da;

	while (scalef<ai) {	
	  int status = gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &scalef, ai, &hh, y);
      if (status != GSL_SUCCESS)
		break;
	}
  }

  /* Add here interpolation to find a(l1=1) */
  return 0;
}

/* Function to be used by qsort, to sort l1>l2>l3 */
int compare(const void * a, const void * b) { 
  double fa = *(const double*) a;
  double fb = *(const double*) b;
  double cc = fb - fa;
  if (cc<0)
	return -1;
  if (cc>0)
	return 1;
  return 0;
} 
  
int main(void) {
  int i, j;
  double lam1, lam2, lam3;
  double delta=3., dx=0.1, dy=0.1; /* d = delta = l1+l2+l3 */
  for(i=0; i<=9; i++) {
	double ix = i*dx;
	for (j=0; j<=9; j++) {
	  double iy = j*dy;
	  lam1 = (delta + 2. * ix + 1. * iy) / 3.;
	  lam2 = (delta - 1. * ix + 1. * iy) / 3.;
	  lam3 = (delta - 1. * ix - 2. * iy) / 3.;

	  double values[] = {lam1, lam2, lam3};
	  qsort(values, 3, sizeof(double), compare);
	  printf("%.5e %.5e %.5e %.5e\n", lam1, lam2, lam3, ellcoll(lam1,lam2,lam3));
	}
  }
  return 0;
}	




/* NB: must find a(lam1=1), a=ai, lam1=y[0] */
/* #define NBIN 200 */
/* static gsl_spline *splineColl = 0x0; */
/* splineColl = gsl_spline_alloc(gsl_interp_cspline, NBIN); */
/* gsl_spline_init(splineColl, X, Y, NBIN); */
/* gsl_spline_eval (splineColl, X, ?); */
/* gsl_spline_free (splineColl); */



/* double ell_int =  gsl_sf_ellint_RD(y1, y2, y3, gsl_mode_t  GSL_PREC_DOUBLE); */
/* printf("%.5e\n", ell_int); */
