#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_cubic_spline_interpolation.h"
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <omp.h>
#include <sys/time.h>


#define N 200000

double dummy_function(double x) {
    // Replace this with your actual function
    return sin(x);
}

int compare_function(const void* a, const void* b) {
    return (*(double*)a - *(double*)b);
}

int main() {

     // Timing variables
    clock_t start_time, end_time;
    double cpu_time_used;
    // Measure start time
    start_time = omp_get_wtime();

    // Sample data points
    double x_data[N];
    double y_data[N];

    // Generate non-equally spaced x values
    for (int i = 0; i < N; ++i) {
        x_data[i] = ((double)rand() / RAND_MAX) * 30.0;  // Example: Random values between 0 and 10
    }

// #ifdef NO_GPU
    // Sort x_data
    gsl_sort(x_data, 1, N);
// #else    
    qsort(x_data, N, sizeof(double), compare_function);
// #endif


    // Populate y_data with dummy function values
    for (int i = 0; i < N; ++i) {
        y_data[i] = dummy_function(x_data[i]);
    }

    // Calculate second derivatives using your custom function
    SplineCoefficients* coefficients_custom = custom_cubic_spline_alloc(N);

    custom_cubic_spline_init(coefficients_custom, x_data, y_data, N);

    // Interpolate at all x values using your custom cubic spline interpolation
    double interpolated_cubic_spline[N];
   

    #pragma omp target map(to: x_data[0:N], coefficients_custom[0:N-1]) map(from: interpolated_cubic_spline[0:N])
    #pragma omp teams distribute parallel for
    for (int i = 0; i < N; ++i) {
        interpolated_cubic_spline[i] = custom_cubic_spline_eval(x_data[i], coefficients_custom, x_data, N - 1);
    }
    
    end_time = omp_get_wtime();
    cpu_time_used = end_time - start_time;


    printf("Time for custom interpolation: %e seconds\n", cpu_time_used);

#ifdef SI_GSL
    // Interpolate at all x values using GSL spline interpolation
    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_cspline, N);
    gsl_spline_init(spline, x_data, y_data, N);

    double interpolated_gsl_spline[N];
    for (int i = 0; i < N; ++i) {
        interpolated_gsl_spline[i] = gsl_spline_eval(spline, x_data[i], acc);
    }
#endif 
    // Save data to a file
    FILE* file = fopen("interpolation_data_si_GPU.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }
#ifdef SI_GSL
    fprintf(file, "x_data\ty_data\tcubic_spline\tgsl_spline\n");
    for (int i = 0; i < N; ++i) {
        fprintf(file, "%.6f\t%.6f\t%.6f\t%.6f\n", x_data[i], y_data[i], interpolated_cubic_spline[i], interpolated_gsl_spline[i]);
    }
    fclose(file);
#else
     fprintf(file, "x_data\ty_data\tcubic_spline\tgsl_spline\n");
    for (int i = 0; i < N; ++i) {
        fprintf(file, "%.6f\t%.6f\t%.6f\n", x_data[i], y_data[i], interpolated_cubic_spline[i]);
    }
    fclose(file);
#endif

#ifdef SI_GSL
    // Free GSL resources
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
#endif

    return 0;
}