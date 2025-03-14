#ifdef GPU_OMP_FULL

#include "cubic_spline_interpolation.h"
#include "pinocchio.h"
#include <omp.h>
#include <assert.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>

CubicSpline *custom_cubic_spline_alloc(const int size)
{
  CubicSpline *spline = (CubicSpline *)malloc(sizeof(CubicSpline));
  spline->size     = size;
  spline->x        = NULL;
  spline->y        = NULL;
  spline->d2y_data = malloc(size * sizeof(double));
  spline->coeff_a  = malloc((size - 1) * sizeof(double));
  spline->coeff_b  = malloc((size - 1) * sizeof(double));
  spline->coeff_c  = malloc((size - 1) * sizeof(double));
  spline->coeff_d  = malloc((size - 1) * sizeof(double));

  return spline;
}

void custom_cubic_spline_init(CubicSpline  *const restrict spline, 
                              const double *const restrict x_data, 
                              const double *const restrict y_data, 
                              const int size)
{    
    // Copy data
    #pragma omp target device(devID)
    {
      spline->x                  = x_data;
      spline->y                  = y_data;
      spline->size               = size;
      spline->d2y_data[0]        = 0.0;  // Boundary condition
      spline->d2y_data[size - 1] = 0.0;  // Boundary condition
    }
    
    // Compute second derivatives
    #pragma omp target teams distribute parallel for device(devID)
    for (int i = 1; i < size - 1; i++)
      {
        const double h1  = (spline->x[i + 0] - spline->x[i - 1]);
        const double h2  = (spline->x[i + 1] - spline->x[i - 0]);
        const double dy1 = (spline->y[i + 0] - spline->y[i - 1]);
        const double dy2 = (spline->y[i + 1] - spline->y[i - 0]);

        spline->d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
      }

    // Compute coefficients
    #pragma omp target teams distribute parallel for device(devID)
    for (int i = 0; i < size - 1; i++)
      {
        double h           = spline->x[i + 1] - spline->x[i];
        spline->coeff_a[i] = spline->y[i];
        spline->coeff_b[i] = (spline->y[i + 1] - spline->y[i]) / h - h / 6.0 * (spline->d2y_data[i + 1] + 2.0 * spline->d2y_data[i]);
        spline->coeff_c[i] = spline->d2y_data[i] * 0.5;
        spline->coeff_d[i] = (spline->d2y_data[i + 1] - spline->d2y_data[i]) / (6.0 * h);
      }

    return;
}

double custom_cubic_spline_eval(CubicSpline *const spline, const double x, int size) {
   
    // Handle extrapolation
    if (x < spline->x[0]) {
        return spline->coeff_a[0] + spline->coeff_b[0] * (x - spline->x[0]);
    } else if (x > spline->x[size - 1]) {
        return spline->coeff_a[size - 2] + spline->coeff_b[size - 2] * (x - spline->x[size - 1]);
    }
    
    // Find the interval
    int k = 0;
    while (k < size - 1 && x > spline->x[k + 1]) {
        k++;
    }

    //  Step 2: Binary search for the interval index 
    // int k = 0;
    // int step = size / 2;
    // while (step > 0) {
    //     int mid = k + step;
    //     if (mid < size - 1 && x > spline->x[mid]) {
    //         k = mid;
    //     }
    //     step /= 2;
    // }

    // Evaluate spline
    double dx = x - spline->x[k];
    return spline->coeff_a[k] + spline->coeff_b[k] * dx + spline->coeff_c[k] * dx * dx + spline->coeff_d[k] * dx * dx * dx;
}

void custom_cubic_spline_free(CubicSpline *spline)
{
  if (spline->x != NULL)
    free(spline->x);

  if (spline->y != NULL)
    free(spline->y);

  if (spline->d2y_data != NULL)
    free(spline->d2y_data);

  if (spline->coeff_a != NULL)
    free(spline->coeff_a);

  if (spline->coeff_b != NULL)
    free(spline->coeff_b);

  if (spline->coeff_c != NULL)
    free(spline->coeff_c);

  if (spline->coeff_d != NULL)
    free(spline->coeff_d);

  if (spline != NULL)
    free(spline);
}

#endif
