#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "interpolation.h"
#include <omp.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include <sys/time.h>

#define N 15

double dummy_function(double x) {
    // Replace this with your actual function
    return pow(x, 2);
    // return sin(x);
}

CubicSpline **my_spline;
gsl_spline **SPLINE;
gsl_interp_accel **ACCEL;


double init_spline() {

#ifdef SI_GSL
    int i;

    SPLINE = (gsl_spline **)calloc(NSPLINES, sizeof(gsl_spline *));
    for (i = 0; i < NSPLINES - 3; i++)
        if (i != SP_COMVDIST && i != SP_DIAMDIST)
            SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS);
        else
            SPLINE[i] = gsl_spline_alloc(gsl_interp_cspline, NBINS - NBB);

    ACCEL = (gsl_interp_accel **)calloc(NSPLINES, sizeof(gsl_interp_accel *));
    for (i = 0; i < NSPLINES; i++)
        ACCEL[i] = gsl_interp_accel_alloc();
#endif


    my_spline = (CubicSpline**)calloc(NSPLINES, sizeof(CubicSpline*));

    for (int i = 0; i < NSPLINES - 3; i++) {
        if (i != SP_COMVDIST && i != SP_DIAMDIST) {
            my_spline[i] = custom_cubic_spline_alloc(NBINS);
        } else {
            my_spline[i] = custom_cubic_spline_alloc(NBINS - NBB);
        }
    }

    double* a, * b;

    a = (double*)malloc(NBINS * sizeof(double));
    b = (double*)malloc(NBINS * sizeof(double));

    for (int i = 0; i < NBINS; i++) {
        a[i] = ((double)rand() / RAND_MAX) * 30.0;
        b[i] = dummy_function(a[i]);
    }

   FILE* spline_data_file = fopen("spline_data.txt", "w"); // Open a file for spline data
    if (spline_data_file == NULL) {
        perror("Error opening spline data file");
        return 1;
    }

    gsl_sort(a, 1, NBINS);  // Sort the 'a' array before initializing the spline
    gsl_sort(b, 1, NBINS); 
    
    for (int i = 0; i < NBINS; i++) {
        // Write spline data to the file
        fprintf(spline_data_file, "%.6f\t%.6f\n", a[i], b[i]);
    }

    fclose(spline_data_file); // Close the file

    custom_cubic_spline_init(my_spline[SP_INVGROW], a, b, NBINS);

    // Allocate memory on the GPU for my_spline and its components
    // #pragma omp target enter data map(alloc: my_spline[0:NSPLINES])
    // for (int i = 0; i < NSPLINES - 3; i++) {
    //     #pragma omp target enter data map(alloc: my_spline[i]->x[0:my_spline[i]->size], \
    //                                        my_spline[i]->y[0:my_spline[i]->size], \
    //                                        my_spline[i]->d2y_data[0:(my_spline[i]->size - 1)], \
    //                                        my_spline[i]->coeff_a[0:(my_spline[i]->size - 1)], \
    //                                        my_spline[i]->coeff_b[0:(my_spline[i]->size - 1)], \
    //                                        my_spline[i]->coeff_c[0:(my_spline[i]->size - 1)], \
    //                                        my_spline[i]->coeff_d[0:(my_spline[i]->size - 1)])
    // }

#ifdef SI_GSL
    gsl_spline_init(SPLINE[SP_INVGROW], a, b, NBINS);
#endif  

    return 0;
}

double custom_spline_eval(CubicSpline* my_spline, double x) {
    // Perform linear extrapolation beyond the x-range limits
    if (x < my_spline->x[0])
        return my_spline->y[0] + (x - my_spline->x[0]) * (my_spline->y[1] - my_spline->y[0]) / (my_spline->x[1] - my_spline->x[0]);
    else if (x > my_spline->x[my_spline->size - 1])
        return my_spline->y[my_spline->size - 1] + (x - my_spline->x[my_spline->size - 1]) *
               (my_spline->y[my_spline->size - 1] - my_spline->y[my_spline->size - 2]) /
               (my_spline->x[my_spline->size - 1] - my_spline->x[my_spline->size - 2]);
    else {

        return custom_cubic_spline_eval(my_spline, x);

    }
}

#ifdef SI_GSL

double GSL_spline_eval(gsl_spline* spline, double x, gsl_interp_accel* accel) {
    /* this function performs linear extrapolation beyond the x-range limits,
    and calls the spline evaluation in between */
    if (x < spline->x[0])
        return spline->y[0] + (x - spline->x[0]) * (spline->y[1] - spline->y[0]) / (spline->x[1] - spline->x[0]);
    else if (x > spline->x[spline->size - 1])
        return spline->y[spline->size - 1] + (x - spline->x[spline->size - 1]) *
               (spline->y[spline->size - 1] - spline->y[spline->size - 2]) / (spline->x[spline->size - 1] - spline->x[spline->size - 2]);
    else {
        return gsl_spline_eval(spline, x, accel);
    }
}
#endif

double use_custom_spline(double a, int ismooth) {
    return custom_cubic_spline_eval(my_spline[SP_INVGROW], a);

}

double use_GSL_spline(double a, int ismooth){
    return GSL_spline_eval(SPLINE[SP_INVGROW], a , ACCEL[SP_INVGROW]);
}

int compare_function(const void* a, const void* b) {
    return (*(double*)a - *(double*)b);
}

int main() {

    // Timing variables
    double start_time, end_time, cpu_time_used;

    // Measure start time
    // start_time = omp_get_wtime();

    // Evaluation values
    double x_data[N];

    // Generate non-equally spaced x values
    for (int i = 0; i < N; ++i) {
        x_data[i] = ((double)rand() / RAND_MAX) * 50.0;  // Example: Random values between 0 and 10
    }
    
    // Calculate actual function values at random x values
    double actual_function_values[N];
    for (int i = 0; i < N; ++i) {
        actual_function_values[i] = dummy_function(x_data[i]);
    }

    // gsl_sort(x_data, 1, N);
    init_spline();

    // Interpolate at all x values using your custom cubic spline interpolation
    double interpolated_cubic_spline[N];
    
    // #pragma omp target map(to: x_data[0:N], my_spline[0:NSPLINES]) map(from: interpolated_cubic_spline[0:N])
    // #pragma omp teams distribute parallel for
    for (int i = 0; i < N; ++i) {
        interpolated_cubic_spline[i] = use_custom_spline(x_data[i], 1);
    }

#ifdef SI_GSL
    // Interpolate at all x values using GSL spline interpolation
    double interpolated_gsl_spline[N];
    for (int i = 0; i < N; ++i) {
        interpolated_gsl_spline[i] = use_GSL_spline(x_data[i], 1);
    }
#endif


    // Measure end time
    // end_time = omp_get_wtime();
    // cpu_time_used = end_time - start_time;

    // gsl_sort(x_data, 1, N);
    // gsl_sort(interpolated_cubic_spline, 1, N);
    // gsl_sort(interpolated_gsl_spline, 1, N);
    // printf("Time for custom interpolation: %e seconds\n", cpu_time_used);

    // Save data to a file
    FILE* file = fopen("interpolation_data.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return 1;
    }

    // fprintf(file, "x_data\tcubic_spline\tgsl_spline\n");
    // for (int i = 0; i < N; ++i) {
    //     fprintf(file, "%.6f\t%.6f\t", x_data[i], interpolated_cubic_spline[i]);
    
    fprintf(file, "x_data\tactual_function\tcubic_spline\n");
    for (int i = 0; i < N; ++i) {
        fprintf(file, "%.6f\t%.6f\t%.6f\t", x_data[i], actual_function_values[i], interpolated_cubic_spline[i]);
    
#ifdef SI_GSL
        fprintf(file, "%.6f", interpolated_gsl_spline[i]);
#endif

        fprintf(file, "\n");
    }

    fclose(file);

#ifdef SI_GSL
    // Free GSL resources
    for (int i = 0; i < NSPLINES - 3; i++) {
        gsl_spline_free(SPLINE[i]);
        gsl_interp_accel_free(ACCEL[i]);
    }
#endif

    return 0;
}