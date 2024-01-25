#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef NO_GPU
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#endif
#include <time.h> 

#define N 120000 // Number of data points

// Data structure for holding spline coefficients
typedef struct {
    double a;
    double b;
    double c;
    double d;
} SplineCoefficients;

double dummy_function(double x) {
    return sin(x);  // You can replace this with any function you want to interpolate
}


void calculate_second_derivatives(double* x_data, double* y_data, double* d2y_data, int size) {
    // Compute second derivatives at data points using finite differences
    #pragma acc parallel loop
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
    #pragma acc parallel loop
    for (int i = 0; i < size - 1; ++i) {
        double h = x_data[i + 1] - x_data[i];
        coefficients[i].a = y_data[i];
        coefficients[i].b = (y_data[i + 1] - y_data[i]) / h - h / 6.0 * (d2y_data[i + 1] + 2.0 * d2y_data[i]);
        coefficients[i].c = d2y_data[i] / 2.0;
        coefficients[i].d = (d2y_data[i + 1] - d2y_data[i]) / (6.0 * h);
    }
}

#pragma acc routine
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

// Comparison function for qsort
int compare_function(const void* a, const void* b) {
    return (*(double*)a > *(double*)b) - (*(double*)a < *(double*)b);
}

int main() {
    // Sample data points
    double x_data[N];
    double y_data[N];
    double d2y_data[N - 1]; // Second derivatives at data points
    SplineCoefficients coefficients[N - 1];

    // Generate non-equally spaced x values
    for (int i = 0; i < N; ++i) {
        x_data[i] = ((double)rand() / RAND_MAX) * 30.0;  // Example: Random values between 0 and 10
    }

#ifdef NO_GPU
    // Sort x_data
    gsl_sort(x_data, 1, N);
#else    
    qsort(x_data, N, sizeof(double), compare_function);
#endif

    // Populate y_data with dummy function values
    for (int i = 0; i < N; ++i) {
        y_data[i] = dummy_function(x_data[i]);
    }

    // Calculate second derivatives
    calculate_second_derivatives(x_data, y_data, d2y_data, N);

    // Calculate cubic spline coefficients
    cubic_spline_coefficients(x_data, y_data, d2y_data, coefficients, N);


    // Timing variables
    clock_t start_time, end_time;
    double cpu_time_used;

    // Interpolate at all x values using cubic spline interpolation
    double interpolated_cubic_spline[N];

    start_time = clock();
    #pragma acc parallel loop
    for (int i = 0; i < N; ++i) {
        interpolated_cubic_spline[i] = cubic_spline_interpolate(x_data[i], x_data, coefficients, N - 1);
    }
    
    end_time = clock();
    cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

#ifdef NO_GPU
    printf("CPU_version Time for custom interpolation: %f seconds\n", cpu_time_used);
#else
    printf("GPU_version Time for custom interpolation: %f seconds\n", cpu_time_used);
#endif

#ifdef NO_GPU
    // Interpolate at all x values using GSL spline interpolation
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline, x_data, y_data, N);

    double interpolated_gsl_spline[N];
    for (int i = 0; i < N; ++i) {
        interpolated_gsl_spline[i] = gsl_spline_eval(spline, x_data[i], acc);
    }

    // Save data to a file
    FILE* file = fopen("interpolation_data_NO_GPU.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    fprintf(file, "x_data\ty_data\tcubic_spline\tgsl_spline\n");
    for (int i = 0; i < N; ++i) {
        fprintf(file, "%.6f\t%.6f\t%.6f\t%.6f\n", x_data[i], y_data[i], interpolated_cubic_spline[i], interpolated_gsl_spline[i]);
    }
#else   
    // Save data to a file
    FILE* file = fopen("interpolation_data_GPU.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    fprintf(file, "x_data\ty_data\tcubic_spline\n");
    for (int i = 0; i < N; ++i) {
        fprintf(file, "%.6f\t%.6f\t%.6f\n", x_data[i], y_data[i], interpolated_cubic_spline[i]);
    }
#endif

    fclose(file);

#ifdef NO_GPU
    // Free GSL resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
#endif

    return 0;
}