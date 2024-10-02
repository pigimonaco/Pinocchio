#ifndef CUBIC_SPLINE
#define CUBIC_SPLINE

typedef struct GPU_CubicSpline
{
  int size;
  double *x;
  double *y;
  double *d2y_data;
  double *coeff_a; 
  double *coeff_b;  
  double *coeff_c; 
  double *coeff_d;  
} CubicSpline;

extern CubicSpline *host_spline;

double custom_cubic_spline_eval(CubicSpline *const spline, const double x);

#if defined(GPU_OMP)
extern CubicSpline gpu_spline;
#pragma omp declare target(gpu_spline, custom_cubic_spline_eval)
#endif // GPU_OMP

void custom_cubic_spline_init(      CubicSpline *const restrict spline,
			      const double      *const restrict x_data,
			      const double      *const restrict y_data,
			      const int                         size);

CubicSpline* custom_cubic_spline_alloc(const int size);
void calculate_second_derivatives(CubicSpline *const spline);
void cubic_spline_coefficients(CubicSpline *const spline);
void custom_cubic_spline_free(CubicSpline *spline);

#endif // CUBIC_SPLINE
