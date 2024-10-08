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
// gsl_spline **SPLINE;
// gsl_interp_accel **ACCEL;

void calculate_second_derivatives(CubicSpline* spline);
void cubic_spline_coefficients(CubicSpline* spline);

// Custom names for cubic spline functions

CubicSpline* custom_cubic_spline_alloc(int size);

void custom_cubic_spline_init(CubicSpline* spline, double* x_data, double* y_data, int size);

double custom_cubic_spline_eval(CubicSpline* spline, double x);

double use_custom_spline(double, int);

void custom_cubic_spline_free(CubicSpline* spline);

#endif
