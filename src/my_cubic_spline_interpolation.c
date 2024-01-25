#include "my_cubic_spline_interpolation.h"
#include <math.h>


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


double cubic_spline_interpolate(double x, double* x_data, SplineCoefficients* coefficients, int size) {
      // Find the interval [x_k, x_{k+1}] containing x
    int k = 0;
    while (k < size - 1 && x > x_data[k + 1]) {
        k++;
    }

    // Evaluate the cubic spline at the point x
    double dx = x - x_data[k];
    return coefficients[k].a + coefficients[k].b * dx + coefficients[k].c * pow(dx, 2) + coefficients[k].d * pow(dx, 3);
}