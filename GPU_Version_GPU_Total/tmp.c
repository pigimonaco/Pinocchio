#include "my_cubic_spline.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void calculate_second_derivatives(CubicSpline *spline) {
    #pragma omp target teams distribute parallel for map(to: spline->x_data[:spline->size], spline->y_data[:spline->size]) map(from: spline->d2y_data[1:spline->size-1])
    for (int i = 1; i < spline->size - 1; ++i) {
        double h1 = spline->x_data[i] - spline->x_data[i - 1];
        double h2 = spline->x_data[i + 1] - spline->x_data[i];
        double dy1 = spline->y_data[i] - spline->y_data[i - 1];
        double dy2 = spline->y_data[i + 1] - spline->y_data[i];

        spline->d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
    }
}

void cubic_spline_coefficients(CubicSpline *spline) {
    #pragma omp target teams distribute parallel for map(to: spline->x_data[:spline->size], spline->y_data[:spline->size], spline->d2y_data[1:spline->size-1]) map(from: spline->coeff_a[:spline->size-1], spline->coeff_b[:spline->size-1], spline->coeff_c[:spline->size-1], spline->coeff_d[:spline->size-1])
    for (int i = 0; i < spline->size - 1; ++i) {
        double h = spline->x_data[i + 1] - spline->x_data[i];
        spline->coeff_a[i] = spline->y_data[i];
        spline->coeff_b[i] =
            (spline->y_data[i + 1] - spline->y_data[i]) / h - h / 6.0 * (spline->d2y_data[i + 1] + 2.0 * spline->d2y_data[i]);
        spline->coeff_c[i] = spline->d2y_data[i] / 2.0;
        spline->coeff_d[i] = (spline->d2y_data[i + 1] - spline->d2y_data[i]) / (6.0 * h);
    }
}

void custom_cubic_spline_init_gpu(double* x_data, double* y_data, int size, CubicSpline *spline) {
    // Allocate memory on GPU
    #pragma omp target enter data map(alloc: spline->x_data[0:size], spline->y_data[0:size], spline->d2y_data[0:(size-1)], \
                                      spline->coeff_a[0:(size-1)], spline->coeff_b[0:(size-1)], \
                                      spline->coeff_c[0:(size-1)], spline->coeff_d[0:(size-1)]) \
                                      map(to: x_data[0:size], y_data[0:size])
    {   
        #pragma omp target teams distribute parallel for
        for (int i = 1; i < spline->size - 1; ++i) {
            double h1 = spline->x_data[i] - spline->x_data[i - 1];
            double h2 = spline->x_data[i + 1] - spline->x_data[i];
            double dy1 = spline->y_data[i] - spline->y_data[i - 1];
            double dy2 = spline->y_data[i + 1] - spline->y_data[i];
            spline->d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
        }

        #pragma omp target distribute parallel for
        for (int i = 0; i < spline->size - 1; ++i) {
            double h = spline->x_data[i + 1] - spline->x_data[i];
            spline->coeff_a[i] = spline->y_data[i];
            spline->coeff_b[i] = (spline->y_data[i + 1] - spline->y_data[i]) / h - h / 6.0 * (spline->d2y_data[i + 1] + 2.0 * spline->d2y_data[i]);
            spline->coeff_c[i] = spline->d2y_data[i] / 2.0;
            spline->coeff_d[i] = (spline->d2y_data[i + 1] - spline->d2y_data[i]) / (6.0 * h);
        }   
    }

#pragma omp declare target
double custom_cubic_spline_eval(CubicSpline *spline, double x) {
    double result = 0.0;

    #pragma omp target teams distribute parallel for map(to: spline->x_data[:spline->size], spline->coeff_a[:spline->size-1], spline->coeff_b[:spline->size-1], spline->coeff_c[:spline->size-1], spline->coeff_d[:spline->size-1]) reduction(+:result)
    for (int k = 0; k < spline->size - 1; ++k) {
        if (x <= spline->x_data[k + 1]) {
            // Evaluate the cubic spline at the point x using local coefficients
            double dx = x - spline->x_data[k];
            double a = spline->coeff_a[k];
            double b = spline->coeff_b[k];
            double c = spline->coeff_c[k];
            double d = spline->coeff_d[k];

            result += a + b * dx + c * dx * dx + d * dx * dx * dx;
            break;  // Exit loop once the interval is found
        }
    }

    return result;
}
#pragma omp end declare target

void custom_cubic_spline_free(CubicSpline *spline) {
    // Free GPU memory
    #pragma omp target exit data map(delete: spline->x_data[0:spline->size], spline->y_data[0:spline->size], spline->d2y_data[0:(spline->size-1)], \
                                      spline->coeff_a[0:(spline->size-1)], spline->coeff_b[0:(spline->size-1)], \
                                      spline->coeff_c[0:(spline->size-1)], spline->coeff_d[0:(spline->size-1)])
}