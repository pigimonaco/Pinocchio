#include "my_cubic_spline_interpolation.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


void calculate_second_derivatives(double* x_data, double* y_data, double* d2y_data, int size) {
    // Compute second derivatives at data points using finite differences
    for (int i = 1; i < size - 1; ++i) {
        double h1 = x_data[i] - x_data[i - 1];
        double h2 = x_data[i + 1] - x_data[i];
        double dy1 = y_data[i] - y_data[i - 1];
        double dy2 = y_data[i + 1] - y_data[i];

        d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
    }
}

void cubic_spline_coefficients(double* x_data, double* y_data, double* d2y_data, SplineCoefficients* coefficients, int size) {
    // Calculate the coefficients for the cubic spline
    for (int i = 0; i < size - 1; ++i) {
        double h = x_data[i + 1] - x_data[i];
        coefficients[i].a = y_data[i];
        coefficients[i].b = (y_data[i + 1] - y_data[i]) / h - h / 6.0 * (d2y_data[i + 1] + 2.0 * d2y_data[i]);
        coefficients[i].c = d2y_data[i] / 2.0;
        coefficients[i].d = (d2y_data[i + 1] - d2y_data[i]) / (6.0 * h);
    }
}


SplineCoefficients* custom_cubic_spline_alloc(int size) {
    SplineCoefficients* coefficients = (SplineCoefficients*)malloc((size - 1) * sizeof(SplineCoefficients));
    if (coefficients == NULL) {
        // Handle allocation failure
        fprintf(stderr,"Error: Unable to allocate memory for coefficients\n");
        exit(EXIT_FAILURE);
    }

    coefficients->size = size - 1;
    coefficients->d2y_data = (double*)malloc((size - 1) * sizeof(double));
    if (coefficients->d2y_data == NULL) {
        // Handle allocation failure
        fprintf(stderr, "Error: Unable to allocate memory for d2y_data\n");
        free(coefficients);
        exit(EXIT_FAILURE);
    }
    return coefficients;
}

void custom_cubic_spline_init(SplineCoefficients* coefficients, double* x_data, double* y_data, int size) {
    coefficients->size = size - 1; 
    coefficients->d2y_data = (double*)malloc((size - 1) * sizeof(double));

    calculate_second_derivatives(x_data, y_data, coefficients->d2y_data, size);
    cubic_spline_coefficients(x_data, y_data, coefficients->d2y_data, coefficients, size);
}


double custom_cubic_spline_eval(double x, SplineCoefficients* coefficients, double* x_data, int size) {
    // Local variables to store coefficients for this thread
    SplineCoefficients local_coefficients;

    // Find the interval [x_k, x_{k+1}] containing x
    int k = 0;
    while (k < size - 1 && x > x_data[k + 1]) {
        k++;
    }

    // Copy coefficients to local variables
    local_coefficients = coefficients[k];

    // Evaluate the cubic spline at the point x using local coefficients
    double dx = x - x_data[k];
    return local_coefficients.a + local_coefficients.b * dx + local_coefficients.c * pow(dx, 2) + local_coefficients.d * pow(dx, 3);
}

void custom_cubic_spline_free(SplineCoefficients* coefficients) {
    free(coefficients);
    free(coefficients -> d2y_data);
}