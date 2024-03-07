#ifndef MY_CUBIC_SPLINE_H
#define MY_CUBIC_SPLINE_H


#include <omp.h> 

typedef struct {
    int size;
    double *x_data;
    double *y_data;
    double *d2y_data;
    double *coeff_a; 
    double *coeff_b;  
    double *coeff_c; 
    double *coeff_d;  
} CubicSpline;


CubicSpline* cubic_spline_alloc(int size);

// Function to initialize cubic spline on the GPU
void custom_cubic_spline_gpu(double* x_data, double* y_data, int size, CubicSpline *spline);

// Declaration of the function to evaluate cubic spline on the GPU
// #pragma omp declare target
double custom_cubic_spline_eval(CubicSpline *spline, double x);
#pragma omp declare target(custom_cubic_spline_eval)

#endif