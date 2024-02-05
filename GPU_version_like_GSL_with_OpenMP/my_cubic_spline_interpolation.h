#ifndef MY_CUBIC_SPLINE_H
#define MY_CUBIC_SPLINE_H

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

void calculate_second_derivatives(CubicSpline spline);
void cubic_spline_coefficients(CubicSpline spline);

// Custom names for cubic spline functions

CubicSpline* custom_cubic_spline_alloc(int size);

void custom_cubic_spline_init(CubicSpline* spline, double* x_data, double* y_data, int size);

#pragma omp declare target
double custom_cubic_spline_eval(CubicSpline spline, double x);
#pragma omp end declare target

void custom_cubic_spline_free(CubicSpline spline);

#endif