#ifndef MY_CUBIC_SPLINE_H
#define MY_CUBIC_SPLINE_H

// Data structure for holding spline coefficients
typedef struct {
    double a;
    double b;
    double c;
    double d;
} SplineCoefficients;

typedef struct {
    SplineCoefficients* coefficients;
    int size;
} CubicSpline;

// #pragma omp declare target
// void calculate_second_derivatives(double* x_data, double* y_data, double* d2y_data, int size);
// #pragma omp end declare target

// #pragma omp declare target
// void cubic_spline_coefficients(double* x_data, double* y_data, double* d2y_data, SplineCoefficients* coefficients, int size);
// #pragma omp end declare target

// #pragma omp declare target
// double cubic_spline_interpolate(double x, double* x_data, SplineCoefficients* coefficients, int size);
// #pragma omp end declare target

// #endif 

#pragma acc routine
void calculate_second_derivatives(double* x_data, double* y_data, double* d2y_data, int size);

#pragma acc routine
void cubic_spline_coefficients(double* x_data, double* y_data, double* d2y_data, SplineCoefficients* coefficients, int size);

#pragma acc routine
double cubic_spline_interpolate(double x, double* x_data, SplineCoefficients* coefficients, int size);

#endif