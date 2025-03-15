#if defined(GPU_OMP) || defined(CUSTOM_INTERPOLATION)

#include "cubic_spline_interpolation.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

void calculate_second_derivatives(CubicSpline *const spline)
{
  // Compute second derivatives at data points using finite differences
  assert(spline != NULL);

  for (int i=1 ; i<(spline->size - 1) ; i++)
    {
      const double h1  = spline->x[i + 0] - spline->x[i - 1];
      const double h2  = spline->x[i + 1] - spline->x[i - 0];
      const double dy1 = spline->y[i + 0] - spline->y[i - 1];
      const double dy2 = spline->y[i + 1] - spline->y[i - 0];

      spline->d2y_data[i] = 6.0 / (h1 + h2) * ((dy2 / h2) - (dy1 / h1));
    }
    spline->d2y_data[0]             = 0.0;  // Boundary condition
    spline->d2y_data[spline->size]  = 0.0;  // Boundary condition
}

void cubic_spline_coefficients(CubicSpline *const spline)
{
  // Calculate the coefficients for the cubic spline
  assert(spline != NULL);

  for (int i=0 ; i<(spline->size - 1) ; ++i)
    {
      const double h     = (spline->x[i + 1] - spline->x[i]);
      
      spline->coeff_a[i] = spline->y[i];
      
      spline->coeff_b[i] = (spline->y[i + 1] - spline->y[i]) / h - h / 6.0 * (spline->d2y_data[i + 1] + 2.0 * spline->d2y_data[i]);

      spline->coeff_c[i] = (spline->d2y_data[i] * 0.5);

      spline->coeff_d[i] = (spline->d2y_data[i + 1] - spline->d2y_data[i]) / (6.0 * h);
  }
}

void custom_cubic_spline_init(      CubicSpline *const restrict spline,
			      const double      *const restrict x_data,
			      const double      *const restrict y_data,
			      const int                         size)
{
  assert((spline != NULL) && (spline->size == size));
  
  memcpy(spline->x, x_data, size * sizeof(double));
  memcpy(spline->y, y_data, size * sizeof(double));

  // Calculate second derivatives and cubic spline coefficients directly within the struct
  calculate_second_derivatives(spline);
  cubic_spline_coefficients(spline);
}

CubicSpline* custom_cubic_spline_alloc(const int size)
{
  CubicSpline *spline = (CubicSpline *)malloc(sizeof(CubicSpline));
  assert(spline != NULL);

  spline->size     = size;
  spline->x        = (double*)malloc(size * sizeof(double));
  spline->y        = (double*)malloc(size * sizeof(double));
  spline->d2y_data = (double*)malloc(size * sizeof(double));
  spline->coeff_a  = (double*)malloc((size - 1) * sizeof(double));
  spline->coeff_b  = (double*)malloc((size - 1) * sizeof(double));
  spline->coeff_c  = (double*)malloc((size - 1) * sizeof(double));
  spline->coeff_d  = (double*)malloc((size - 1) * sizeof(double));

  assert((spline->x        != NULL) &&
	 (spline->y        != NULL) &&
	 (spline->d2y_data != NULL) &&
	 (spline->coeff_a  != NULL) &&
	 (spline->coeff_b  != NULL) &&
	 (spline->coeff_c  != NULL) &&
	 (spline->coeff_d  != NULL));
  
  return spline;
}

double custom_cubic_spline_eval(CubicSpline *const spline,
				const double       x)
{
  if (x < spline->x[0])
    {
      // Linear extrapolation at the lower bound
      const double dx = x - spline->x[0];
      const double a  = spline->coeff_a[0];
      const double b  = spline->coeff_b[0];

      const double ret = a + (b * dx); // Linear extrapolation using the slope at the first interval

      return ret;
    }
  else if (x > spline->x[spline->size - 1])
    {
      // Linear extrapolation at the upper bound
      const double dx = x - spline->x[spline->size - 1];
      const double a  = spline->coeff_a[spline->size - 2];
      const double b  = spline->coeff_b[spline->size - 2];

      const double ret = a + (b * dx); // Linear extrapolation using the slope at the last interval
      
      return ret;
    }
  else
    {
      // Find the interval [x_k, x_{k+1}] containing x
      int k = 0;
      while ((k < (spline->size - 1)) && (x > spline->x[k + 1]))
	{
	  k++;
	}
      // Evaluate the cubic spline at the point x using local coefficients
      const double dx = x - spline->x[k];
      const double a = spline->coeff_a[k];
      const double b = spline->coeff_b[k];
      const double c = spline->coeff_c[k];
      const double d = spline->coeff_d[k];

      const double ret = (a + b * dx + c * (dx * dx) + d * (dx * dx * dx));
      
      return ret;
    }
}

void custom_cubic_spline_free(CubicSpline *spline)
{
  if (spline->x)
    free(spline->x);

  if (spline->y)
    free(spline->y);

  if (spline->d2y_data)
    free(spline->d2y_data);

  if (spline->coeff_a)
    free(spline->coeff_a);

  if (spline->coeff_b)
    free(spline->coeff_b);

  if (spline->coeff_c)
    free(spline->coeff_c);

  if (spline->coeff_d)
    free(spline->coeff_d);

  return;
}

#endif // GPU_OMP