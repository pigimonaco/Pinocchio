#include "interpolation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void calculate_second_derivatives(CubicSpline* spline) {
    // Compute second derivatives at data points using finite differences
    //  #pragma omp target teams distribute parallel for map(to: spline->x[:spline->size], spline->y[:spline->size]) map(from: spline->d2y_data[1:spline->size-1])
    // #pragma omp target teams distribute parallel for
    for (int i = 1; i < spline->size - 1; ++i) {
        double h1 = spline->x[i] - spline->x[i - 1];
        double h2 = spline->x[i + 1] - spline->x[i];
        double dy1 = spline->y[i] - spline->y[i - 1];
        double dy2 = spline->y[i + 1] - spline->y[i];

        spline->d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
    }
}

void cubic_spline_coefficients(CubicSpline* spline) {
    // Calculate the coefficients for the cubic spline
    // #pragma omp target teams distribute parallel for map(to: spline->x[:spline->size], spline->y[:spline->size], spline->d2y_data[1:spline->size-1]) map(from: spline->coeff_a[:spline->size-1], spline->coeff_b[:spline->size-1], spline->coeff_c[:spline->size-1], spline->coeff_d[:spline->size-1])
    // #pragma omp target teams distribute parallel for
    for (int i = 0; i < spline->size - 1; ++i) {
        double h = spline->x[i + 1] - spline->x[i];
        spline->coeff_a[i] = spline->y[i];
        spline->coeff_b[i] =
            (spline->y[i + 1] - spline->y[i]) / h - h / 6.0 * (spline->d2y_data[i + 1] + 2.0 * spline->d2y_data[i]);
        spline->coeff_c[i] = spline->d2y_data[i] / 2.0;
        spline->coeff_d[i] = (spline->d2y_data[i + 1] - spline->d2y_data[i]) / (6.0 * h);
    }
}

void custom_cubic_spline_init(CubicSpline* spline, double* x_data, double* y_data, int size) {


    spline->size = size;
    spline->x = (double*)malloc(size * sizeof(double));
    memcpy(spline->x, x_data, size * sizeof(double));

    spline->y = (double*)malloc(size * sizeof(double));
    memcpy(spline->y, y_data, size * sizeof(double));

    spline->d2y_data = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_a  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_b  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_c  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_d  = (double*)malloc((size - 1) * sizeof(double));

    
    // Calculate second derivatives and cubic spline coefficients directly within the struct
    calculate_second_derivatives(spline);
    cubic_spline_coefficients(spline);
}

CubicSpline* custom_cubic_spline_alloc(int size) {
    CubicSpline* spline = (CubicSpline*)malloc(sizeof(CubicSpline));

    spline->size = size;
    spline->x        = (double*)malloc(size * sizeof(double));
    spline->y        = (double*)malloc(size * sizeof(double));
    spline->d2y_data = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_a  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_b  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_c  = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_d  = (double*)malloc((size - 1) * sizeof(double));

    return spline;
}


double custom_cubic_spline_eval(CubicSpline *spline, double x) {
    if (x < spline->x[0]) {
        // Linear extrapolation at the lower bound
        double dx = x - spline->x[0];
        double a = spline->coeff_a[0];
        double b = spline->coeff_b[0];
        return a + b * dx; // Linear extrapolation using the slope at the first interval
    } else if (x > spline->x[spline->size - 1]) {
        // Linear extrapolation at the upper bound
        int k = spline->size - 2;
        double dx = x - spline->x[spline->size - 1];
        double a = spline->coeff_a[k];
        double b = spline->coeff_b[k];
        return a + b * dx; // Linear extrapolation using the slope at the last interval
    } else {
        // Find the interval [x_k, x_{k+1}] containing x
        int k = 0;
        while (k < spline->size - 1 && x > spline->x[k + 1]) {
            k++;
        }
        // Evaluate the cubic spline at the point x using local coefficients
        double dx = x - spline->x[k];
        double a = spline->coeff_a[k];
        double b = spline->coeff_b[k];
        double c = spline->coeff_c[k];
        double d = spline->coeff_d[k];

        return a + b * dx + c * pow(dx, 2) + d * pow(dx, 3);
    }
}

void custom_cubic_spline_free(CubicSpline* spline) {
    free(spline->x);
    free(spline->y);
    free(spline->d2y_data);
    free(spline->coeff_a);
    free(spline->coeff_b);
    free(spline->coeff_c);
    free(spline->coeff_d);

}