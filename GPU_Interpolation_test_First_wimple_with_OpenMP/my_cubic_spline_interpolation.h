#ifndef MY_CUBIC_SPLINE_H
#define MY_CUBIC_SPLINE_H

typedef struct {
    double a;
    double b;
    double c;
    double d;
    int size;          // New field to store the size of the array
    double* d2y_data;  // Change this to a pointer
} SplineCoefficients;


void calculate_second_derivatives(double* x_data, double* y_data, double* d2y_data, int size);
void cubic_spline_coefficients(double* x_data, double* y_data, double* d2y_data, SplineCoefficients* coefficients, int size);


// Custom names for cubic spline functions

SplineCoefficients* custom_cubic_spline_alloc(int size);

void custom_cubic_spline_init(SplineCoefficients* coefficients, double* x_data, double* y_data, int size);

#pragma omp declare target
double custom_cubic_spline_eval(double x, SplineCoefficients* coefficients, double* x_data, int size);
#pragma omp end declare target

void custom_cubic_spline_free(SplineCoefficients* coefficients);

#endif