#ifndef MY_CUBIC_SPLINE_H
#define MY_CUBIC_SPLINE_H
#include <gsl/gsl_spline.h>

typedef struct {
    int size;
    double *x;
    double *y;
    double *d2y_data;
    double *coeff_a; 
    double *coeff_b;  
    double *coeff_c; 
    double *coeff_d;  
} CubicSpline;


CubicSpline **my_spline;
gsl_spline **SPLINE;
gsl_interp_accel **ACCEL;

void calculate_second_derivatives(CubicSpline* spline);
void cubic_spline_coefficients(CubicSpline* spline);

// Custom names for cubic spline functions

CubicSpline* custom_cubic_spline_alloc(int size);

void custom_cubic_spline_init(CubicSpline* spline, double* x_data, double* y_data, int size);

double custom_cubic_spline_eval(CubicSpline* spline, double x);
#pragma omp declare target(custom_cubic_spline_eval)

double use_custom_spline(double, int);
#pragma omp declare target(my_spline)
#pragma omp declare target(use_custom_spline)

void custom_cubic_spline_free(CubicSpline* spline);


#ifndef SCALE_DEPENDENT
#define NkBINS 1
#else
#define NkBINS 10
#define LOGKMIN ((double)-3.0)
#define DELTALOGK ((double)0.5)
#endif

#define SP_TIME     0
#define SP_INVTIME  1
#define SP_COMVDIST 2
#define SP_DIAMDIST 3
#define SP_INVGROW  4
#define SP_MASSVAR  5
#define SP_DISPVAR  6
#define SP_RADIUS   7
#define SP_DVARDR   8
#define NSPLINES (9+8*NkBINS+3)
#define NBINS 1000
#define NBB 10

#endif