#include "my_cubic_spline_interpolation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void calculate_second_derivatives(CubicSpline spline) {
    // Compute second derivatives at data points using finite differences
    for (int i = 1; i < spline.size - 1; ++i) {
        double h1 = spline.x_data[i] - spline.x_data[i - 1];
        double h2 = spline.x_data[i + 1] - spline.x_data[i];
        double dy1 = spline.y_data[i] - spline.y_data[i - 1];
        double dy2 = spline.y_data[i + 1] - spline.y_data[i];

        spline.d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
    }
}

void cubic_spline_coefficients(CubicSpline spline) {
    // Calculate the coefficients for the cubic spline
    for (int i = 0; i < spline.size - 1; ++i) {
        double h = spline.x_data[i + 1] - spline.x_data[i];
        spline.coeff_a[i] = spline.y_data[i];
        spline.coeff_b[i] =
            (spline.y_data[i + 1] - spline.y_data[i]) / h - h / 6.0 * (spline.d2y_data[i + 1] + 2.0 * spline.d2y_data[i]);
        spline.coeff_c[i] = spline.d2y_data[i] / 2.0;
        spline.coeff_d[i] = (spline.d2y_data[i + 1] - spline.d2y_data[i]) / (6.0 * h);
    }
}


void custom_cubic_spline_init(CubicSpline* spline, double* x_data, double* y_data, int size) {
    spline->size = size;
    spline->x_data = (double*)malloc(size * sizeof(double));
    memcpy(spline->x_data, x_data, size * sizeof(double));

    spline->y_data = (double*)malloc(size * sizeof(double));
    memcpy(spline->y_data, y_data, size * sizeof(double));

    spline->d2y_data = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_a = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_b = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_c = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_d = (double*)malloc((size - 1) * sizeof(double));

    // Calculate second derivatives and cubic spline coefficients directly within the struct
    calculate_second_derivatives(*spline);
    cubic_spline_coefficients(*spline);
}

CubicSpline* custom_cubic_spline_alloc(int size) {
    CubicSpline* spline = (CubicSpline*)malloc(sizeof(CubicSpline));
    if (spline == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for cubic spline\n");
        exit(EXIT_FAILURE);
    }

    spline->size = size;
    spline->x_data = (double*)malloc(size * sizeof(double));
    spline->y_data = (double*)malloc(size * sizeof(double));
    spline->d2y_data = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_a = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_b = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_c = (double*)malloc((size - 1) * sizeof(double));
    spline->coeff_d = (double*)malloc((size - 1) * sizeof(double));

    if (spline->x_data == NULL || spline->y_data == NULL || spline->d2y_data == NULL ||
        spline->coeff_a == NULL || spline->coeff_b == NULL || spline->coeff_c == NULL || spline->coeff_d == NULL) {
        fprintf(stderr, "Error: Unable to allocate memory for cubic spline data\n");
        free(spline);
        exit(EXIT_FAILURE);
    }

    return spline;
}

double custom_cubic_spline_eval(CubicSpline spline, double x) {
    // Find the interval [x_k, x_{k+1}] containing x
    int k = 0;
    while (k < spline.size - 1 && x > spline.x_data[k + 1]) {
        k++;
    }

    // Evaluate the cubic spline at the point x using local coefficients
    double dx = x - spline.x_data[k];
    double a = spline.coeff_a[k];
    double b = spline.coeff_b[k];
    double c = spline.coeff_c[k];
    double d = spline.coeff_d[k];

    return a + b * dx + c * pow(dx, 2) + d * pow(dx, 3);
}


void custom_cubic_spline_free(CubicSpline spline) {
    free(spline.x_data);
    free(spline.y_data);
    free(spline.d2y_data);
    free(spline.coeff_a);
    free(spline.coeff_b);
    free(spline.coeff_c);
    free(spline.coeff_d);
}