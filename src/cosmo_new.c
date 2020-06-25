/*****************************************************************
 *                        PINOCCHIO  V4.1                        *
 *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
 *****************************************************************
 
 This code was written by
 Pierluigi Monaco
 Copyright (C) 2016
 
 web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
 
 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "pinocchio.h"
#include "def_splines.h"


#define NVAR (1+NkBINS*8)
#define NBINS 210
#define NBB 10
#ifdef NORADIATION
#define OMEGARAD_H2 ((double)0.0)
#else
#define OMEGARAD_H2 ((double)4.2e-5)
#endif
#define TOLERANCE ((double)1.e-4)
#define UnitLength_in_cm ((double)3.085678e24)
#define HUBBLETIME_GYR ((double)3.085678e24/(double)1.e7/(double)3.1558150e16)
#define DELTA_C ((double)1.686)
//#define FOMEGA_GAMMA 0.554

static int Today;
static int WhichSpectrum, NPowerTable=0, NtabEoS=0;
static double PkNorm, MatterDensity, OmegaK, OmegaRad;

/* declaration of gsl quantities */
gsl_function Function;

#ifdef SCALE_DEPENDENT
double kmin,kmax;
#endif

int system_of_ODEs(double, const double [], double*, void* );
int read_TabulatedEoS(void);
int initialize_PowerSpectrum(void);
int normalize_PowerSpectrum(void);
int read_Pk_from_file(void);
double IntegrandForEoS(double, void*);
double DE_EquationOfState(double);
double IntegrandComovingDistance(double, void*);
double ComputeMassVariance(double);
double ComputeDisplVariance(double);


/**************************/
/* INITIALIZATION SECTION */
/**************************/

int initialize_cosmology()
{

  /*
    Computes the following functions:
    Scale factor, growth first, second and third-order LPT, cosmic time, 
    comoving and diameter distance on a grid of values to be interpolated
  */

  double ode_param;
  double y[NVAR], x1, x2, hh, norm, result, error, SqrtOmegaK, R0, k;
  int status=GSL_SUCCESS, i, j;
#ifdef SCALE_DEPENDENT
  int ik;
#endif
  char filename[BLENGTH];
  FILE *fd;
  double log_amin=-4.,dloga=-log_amin/(double)(NBINS-NBB);

  double *scalef, *cosmtime, *grow1, *grow2, *IntEoS, *comvdist, *diamdist,
    *fomega1, *fomega2, *grow31, *grow32, *fomega31, *fomega32;

  const gsl_odeiv2_step_type *T = gsl_odeiv2_step_rkf45;
  gsl_odeiv2_step *ode_s = gsl_odeiv2_step_alloc(T,NVAR);
  gsl_odeiv2_control *ode_c = gsl_odeiv2_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
  gsl_odeiv2_evolve *ode_e = gsl_odeiv2_evolve_alloc(NVAR);
  gsl_odeiv2_system ode_sys = {system_of_ODEs, jac, NVAR, (void*)&ode_param};

#ifdef MOD_GRAV_FR
  H_over_c = 100. / SPEEDOFLIGHT;
#endif
  OmegaRad = OMEGARAD_H2 / params.Hubble100 / params.Hubble100;
  OmegaK = 1.0-params.Omega0 - params.OmegaLambda - OmegaRad;
  SqrtOmegaK = sqrt(fabs(OmegaK));
  MatterDensity = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0;
  if (params.DEw0==-1 && params.DEwa==0 && !strcmp(params.TabulatedEoSfile,"no"))
    params.simpleLambda=1;
  else
    params.simpleLambda=0;

  /* allocation of (most) splines */
  SPLINE = (gsl_spline **)calloc(NSPLINES, sizeof(gsl_spline *));
  for (i=0; i<NSPLINES-3; i++)
    if (i!=SP_COMVDIST && i!=SP_DIAMDIST)
      SPLINE[i] = gsl_spline_alloc (gsl_interp_cspline, NBINS);
    else
      SPLINE[i] = gsl_spline_alloc (gsl_interp_cspline, NBINS-NBB);
  ACCEL = (gsl_interp_accel **)calloc(NSPLINES, sizeof(gsl_interp_accel *));
  for (i=0; i<NSPLINES; i++)
    ACCEL[i] = gsl_interp_accel_alloc();

  /* if needed, read the tabulated power spectrum
     and initialize its spline */
  if (strcmp(params.FileWithInputSpectrum,"no"))
    if (read_Pk_from_file())
      return 1;

  /* if needed, read the tabulated Equation of State of the dark energy 
     and initialize its spline, then compute the integrand for the DE EoS */
  if (!params.simpleLambda)
    {
      if (strcmp(params.TabulatedEoSfile,"no"))
	{
	  if (read_TabulatedEoS())  /* this allocates and initializes the spline */
	    return 1;
	}
      else
	NtabEoS=0;

      SPLINE[SP_INTEOS] = gsl_spline_alloc(gsl_interp_cspline, NBINS);

      scalef  = (double*)malloc(NBINS * sizeof(double));
      IntEoS  = (double*)malloc(NBINS * sizeof(double));
      Function.function = &IntegrandForEoS;
      for (i=0; i<NBINS; i++)
	{
	  x2=pow(10., log_amin + (i+1)*dloga);
	  gsl_integration_qags(&Function, x2, 1.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  scalef[i] = log_amin + (i+1)*dloga;
	  IntEoS[i] = result;
	}
      
      gsl_spline_init(SPLINE[SP_INTEOS], scalef, IntEoS, NBINS);

      free(IntEoS);
      free(scalef);
    }


#ifdef SCALE_DEPENDENT
  kmin=pow(10.,LOGKMIN);
  kmax=pow(10.,LOGKMIN+(NkBINS-1)*DELTALOGK);
#endif


  /* allocation of vectors for interpolation */
  scalef    = (double*)malloc(NBINS * sizeof(double));
  cosmtime  = (double*)malloc(NBINS * sizeof(double));
  comvdist  = (double*)malloc((NBINS-NBB) * sizeof(double));
  diamdist  = (double*)malloc((NBINS-NBB)* sizeof(double));
  grow1     = (double*)malloc(NBINS * NkBINS * sizeof(double));
  grow2     = (double*)malloc(NBINS * NkBINS * sizeof(double));
  grow31    = (double*)malloc(NBINS * NkBINS * sizeof(double));
  grow32    = (double*)malloc(NBINS * NkBINS * sizeof(double));
  fomega1   = (double*)malloc(NBINS * NkBINS * sizeof(double));
  fomega2   = (double*)malloc(NBINS * NkBINS * sizeof(double));
  fomega31  = (double*)malloc(NBINS * NkBINS * sizeof(double));
  fomega32  = (double*)malloc(NBINS * NkBINS * sizeof(double));


  /* ICs for the runge-kutta integration */
  x1   = pow(10., log_amin -2.);   /*  initial value of scale factor a */
  y[0] = 2. / 3. * pow(x1, 1.5);   /*  initial value of t(a)*Hubble0   */

  /* this is valid for scale-independent and scale-dependent functions */
  for (j=0; j<NkBINS; j++) {
    y[1+j*8] = 1.0;                  /*  initial value of dD1/da   */
    y[2+j*8] = x1;                   /*  initial value of D1(a)    */
    y[3+j*8] = -6./7. * x1;          /*  initial value of dD2/da   */
    y[4+j*8] = -3./7. * x1*x1;       /*  initial value of D2(a)    */
    y[5+j*8] = -x1*x1;               /*  initial value of dD^3a/da */
    y[6+j*8] = -1./3. * x1*x1*x1;    /*  initial value of D^3a     */
    y[7+j*8] = 10./7. * x1*x1;       /*  initial value of dD^3b/da */
    y[8+j*8] = 10./21. * x1*x1*x1;   /*  initial value of D^3b     */
  }

  hh=x1/10.;                       /*  initial guess of time-step */

  /* this function will be integrated within the loop */
  Function.function = &IntegrandComovingDistance;

  /***********************************************/
  /* ODE integration of time-dependent functions */
  /***********************************************/
  for (i=0, Today=0; i<NBINS; i++)
    {
      x2=pow(10., log_amin + i*dloga);
      if (fabs(log_amin + i*dloga)<dloga/10.)
	x2=1.0;

      /* integration of ODE system */
      while (x1<x2 && status==GSL_SUCCESS)
	{
	  status=gsl_odeiv2_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
	  if (status!=GSL_SUCCESS)
	    {
	      printf("ERROR on task %d: integration of cosmological quantities failed\n",ThisTask);
	      fflush(stdout);
	      return 1;
	    }
	}

      scalef[i]   = x2;
      cosmtime[i] = log10(y[0]*HUBBLETIME_GYR/params.Hubble100);
      if (!Today && x2>=1.)
	Today=i;

      for (j=0; j<NkBINS; j++) /* First-order growth rate */
	grow1[i+j*NBINS]  = y[2+j*8];
      for (j=0; j<NkBINS; j++) /* Second-order growth rate */
	grow2[i+j*NBINS]  = -y[4+j*8];
      for (j=0; j<NkBINS; j++) /* Third-order first growth rate */
	grow31[i+j*NBINS] = -y[6+j*8] / 3.;
      for (j=0; j<NkBINS; j++) /* Third-order second growth rate */
	grow32[i+j*NBINS] = y[8+j*8] / 4.;

      for (j=0; j<NkBINS; j++) /* First-order f(Omega) */
	fomega1[i+j*NBINS]  = x2 * y[1+j*8] / y[2+j*8];
      for (j=0; j<NkBINS; j++) /* Second-order f(Omega) */
	fomega2[i+j*NBINS]  = x2 * y[3+j*8] / y[4+j*8];
      for (j=0; j<NkBINS; j++) /* Third-order first f(Omega) */
	fomega31[i+j*NBINS] = x2 * y[5+j*8] / y[6+j*8];
      for (j=0; j<NkBINS; j++) /* Third-order second f(Omega) */
	fomega32[i+j*NBINS] = x2 * y[7+j*8] / y[8+j*8];

      /* calculation of comoving distance (Mpc) for a generic cosmology (i.e. flat, open or closed) 
	 and generic equation of state for the DE component  */
      if (i<NBINS-NBB)
	{
	  gsl_integration_qags(&Function, 0.0, 1./x2-1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  comvdist[i] = SPEEDOFLIGHT*result;
	  if (fabs(OmegaK) < 1.e-4)
	    diamdist[i] = x2 * comvdist[i];
	  else if (OmegaK <0)
	    {
	      R0 = SPEEDOFLIGHT/params.Hubble100/100. / SqrtOmegaK;
	      diamdist[i] = x2 * R0 * sin(comvdist[i]/R0);
	    }
	  else
	    {
	      R0 = SPEEDOFLIGHT/params.Hubble100/100. / SqrtOmegaK;
	      diamdist[i] = x2 * R0 * sinh(comvdist[i]/R0);
	    }
	}

      if (!params.simpleLambda)
	gsl_integration_qags(&Function, x2, 1.0, 0.0, TOLERANCE, NWINT, workspace, &IntEoS[i], &error);

      /* closing the loop on integrations */
      x1=x2;
    }

  /* normalization of the first- and second- order growth rate */
  /* this is valid for LambdaCDM; for scale-dependent growth due to modified gravity, 
     the power spectrum is given as the LambdaCDM P(k) extrapolated at z=0, but this
     is valid only at high redshift; this normalization is still correct 
     when the k=0 growth rate (identical to LambdaCDM) is used at all scales */
  norm =  grow1[Today];
  for (i=0; i<NBINS*NkBINS; i++)
    {
      grow1[i] /= norm;
      grow2[i] /= norm * norm;
      grow31[i] /= norm * norm * norm;
      grow32[i] /= norm * norm * norm;
    }

  /* these quantities will be interpolated logarithmically */
  for (i=0; i<NBINS; i++)
    scalef[i]   = log10(scalef[i]);

  for (j=0; j<NkBINS; j++) /* First-order growth rate */
    for (i=0; i<NBINS; i++)
      {
	grow1[i+j*NBINS]  = log10(grow1[i+j*NBINS]);
	grow2[i+j*NBINS]  = log10(grow2[i+j*NBINS]);
	grow31[i+j*NBINS] = log10(grow31[i+j*NBINS]);
	grow32[i+j*NBINS] = log10(grow32[i+j*NBINS]);
      }

  /* initialization of spline interpolations of time-dependent quantities */
  gsl_spline_init(SPLINE[SP_TIME], scalef, cosmtime, NBINS);
  gsl_spline_init(SPLINE[SP_INVTIME], cosmtime, scalef, NBINS);
  gsl_spline_init(SPLINE[SP_COMVDIST], scalef, comvdist, NBINS-NBB);
  gsl_spline_init(SPLINE[SP_DIAMDIST], scalef, diamdist, NBINS-NBB);
#ifdef ELL_CLASSIC
  gsl_spline_init(SPLINE[SP_INVGROW], grow1, scalef, NBINS);  // ????
#endif

  for (j=0; j<NkBINS; j++)
    {
      gsl_spline_init(SPLINE[SP_GROW1+j],  scalef, grow1+j*NBINS,  NBINS);
      gsl_spline_init(SPLINE[SP_GROW2+j],  scalef, grow2+j*NBINS,  NBINS);
      gsl_spline_init(SPLINE[SP_GROW31+j], scalef, grow31+j*NBINS, NBINS);
      gsl_spline_init(SPLINE[SP_GROW32+j], scalef, grow32+j*NBINS, NBINS);

      gsl_spline_init(SPLINE[SP_FOMEGA1+j],  scalef, fomega1+j*NBINS,  NBINS);
      gsl_spline_init(SPLINE[SP_FOMEGA2+j],  scalef, fomega2+j*NBINS,  NBINS);
      gsl_spline_init(SPLINE[SP_FOMEGA31+j], scalef, fomega31+j*NBINS, NBINS);
      gsl_spline_init(SPLINE[SP_FOMEGA32+j], scalef, fomega32+j*NBINS, NBINS);
    }

  /* deallocation of vectors for interpolation */
  free(fomega32);
  free(fomega31);
  free(fomega2);
  free(fomega1);
  free(grow32);
  free(grow31);
  free(grow2);
  free(grow1);
  if (!params.simpleLambda)
    free(IntEoS);
  free(diamdist);
  free(comvdist);
  free(cosmtime);
  free(scalef);

  /* initialization of power spectrum */
  if (initialize_PowerSpectrum())
    return 1;

  /* initialization of mass variance with Gaussian filter */
  WindowFunctionType=0;
  if (initialize_MassVariance())
    return 1;

  /* write out cosmological quantities */
  if (!ThisTask)
    {

      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".cosmology.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Cosmological quantities used in PINOCCHIO (h=%f)\n",params.Hubble100);
      fprintf(fd,"# TIME-DEPENDENT QUANTITIES\n");
      fprintf(fd,"# 1: scale factor\n");
      fprintf(fd,"# 2: cosmic time (Gyr)\n");
      fprintf(fd,"# 3: comoving distance (Mpc)\n");
      fprintf(fd,"# 4: diameter distance (Mpc)\n");
      fprintf(fd,"# 5: Omega matter\n");
      fprintf(fd,"# 6: dark energy EOS\n");
      fprintf(fd,"# 7: linear growth rate\n");
      fprintf(fd,"# 8: 2nd-order growth rate\n");
      fprintf(fd,"# 9: first 3rd-order growth rate\n");
      fprintf(fd,"#10: second 3rd-order growth rate\n");
      fprintf(fd,"#11: linear d ln D/d ln a\n");
      fprintf(fd,"#12: 2nd-order d ln D/d ln a\n");
      fprintf(fd,"#13: first 3rd-order d ln D/d ln a\n");
      fprintf(fd,"#14: second 3rd-order d ln D/d ln a\n");
      fprintf(fd,"# SCALE-DEPENDENT QUANTITIES\n");
      fprintf(fd,"#15: smoothing scale (Mpc)\n");
      fprintf(fd,"#16: mass variance\n");
      fprintf(fd,"#17: variance of displacements\n");
      fprintf(fd,"#18: d Log sigma^2 / d Log R\n");
      fprintf(fd,"# POWER SPECTRUM\n");
      fprintf(fd,"#19: k (true Mpc^-1)\n");
      fprintf(fd,"#20: P(k)\n");
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  k=pow(10.,-4.0+(double)i/(double)NBINS*6.0);
	  fprintf(fd," %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg %12lg   %12lg %12lg %12lg %12lg   %12lg %12lg\n",
		  pow(10.,SPLINE[SP_TIME]->x[i]),
		  pow(10.,SPLINE[SP_TIME]->y[i]),
		  (i<NBINS-NBB ? SPLINE[SP_COMVDIST]->y[i] : 0.0),
		  (i<NBINS-NBB ? SPLINE[SP_DIAMDIST]->y[i] : 0.0),
		  OmegaMatter(1./pow(10.,SPLINE[SP_TIME]->x[i])-1.0),
		  (!params.simpleLambda ? -1 : (NtabEoS? SPLINE[SP_EOS]->y[i] : DE_EquationOfState(pow(10.,SPLINE[SP_TIME]->x[i])))),
		  pow(10.,SPLINE[SP_GROW1]->y[i]),
		  pow(10.,SPLINE[SP_GROW2]->y[i]),
		  pow(10.,SPLINE[SP_GROW31]->y[i]),
		  pow(10.,SPLINE[SP_GROW32]->y[i]),
		  SPLINE[SP_FOMEGA1]->y[i],
		  SPLINE[SP_FOMEGA2]->y[i],
		  SPLINE[SP_FOMEGA31]->y[i],
		  SPLINE[SP_FOMEGA32]->y[i],
		  pow(10.,SPLINE[SP_MASSVAR]->x[i]),
		  pow(10.,SPLINE[SP_MASSVAR]->y[i]),
		  pow(10.,SPLINE[SP_DISPVAR]->y[i]),
		  SPLINE[SP_DVARDR]->y[i],
		  k, PowerSpectrum(k));
	  }
      fclose(fd);
#ifdef SCALE_DEPENDENT
      strcpy(filename,"pinocchio.");
      strcat(filename,params.RunFlag);
      strcat(filename,".scaledep.out");

      fd=fopen(filename,"w");

      fprintf(fd,"# Scale-dependent growth rates\n");
      fprintf(fd,"# Scales considered: ");
      for (ik=0; ik<NkBINS; ik++)
	{
	  if (!ik)
	    k=0.0;
	  else
	    k = pow(10., LOGKMIN + ik*DELTALOGK);
	  if (ik==NkBINS-1)
	    fprintf(fd,"%d) k=%8.5f\n",ik+1,k);
	  else
	    fprintf(fd,"%d) k=%8.5f, ",ik+1,k);
	}

      fprintf(fd,"# 1: scale factor\n");
      fprintf(fd,"# %d-%d: linear growth rate\n",                     2,  NkBINS+1);
      fprintf(fd,"# %d-%d: 2nd-order growth rate\n",           NkBINS+2,2*NkBINS+1);
      fprintf(fd,"# %d-%d: first 3rd-order growth rate\n",   2*NkBINS+2,3*NkBINS+1);
      fprintf(fd,"# %d-%d: second 3rd-order growth rate\n",  3*NkBINS+2,4*NkBINS+1);
      fprintf(fd,"# %d-%d: linear d ln D/d ln a\n",          4*NkBINS+2,5*NkBINS+1);
      fprintf(fd,"# %d-%d: 2nd-order d ln D/d ln a\n",       5*NkBINS+2,6*NkBINS+1);
      fprintf(fd,"# %d-%d: first 3rd-order d ln D/d ln a\n", 6*NkBINS+2,7*NkBINS+1);
      fprintf(fd,"# %d-%d: second 3rd-order d ln D/d ln a\n",7*NkBINS+2,8*NkBINS+1);
      fprintf(fd,"#\n");

      for (i=0; i<NBINS; i++)
	{
	  fprintf(fd," %12lg",pow(10.,SPLINE[SP_TIME]->x[i]));
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW1+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW2+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW31+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_GROW32+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_FOMEGA1+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_FOMEGA2+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_FOMEGA31+ik]->y[i]));
	  fprintf(fd,"   ");
	  for (ik=0; ik<NkBINS; ik++)
	    fprintf(fd," %12lg",pow(10.,SPLINE[SP_FOMEGA32+ik]->y[i]));
	  fprintf(fd,"\n");
	}
      fclose(fd);
#endif
    }  

  return 0;
}

#ifdef MOD_GRAV_FR
/* scale-dependent functions for 2LPT term */

double mu(double a, double k) {
  double B1, B2, emme;
  B1   = params.Omega0 / pow(a, 3.) + 4. * params.OmegaLambda;
  B2   = params.Omega0 + 4. * params.OmegaLambda;
  emme = 0.5 * H_over_c * H_over_c * pow(B1, 3.) / (B2 * B2 * FR0); 
  return 1. + k * k / 3. / (k * k + a * a * emme);
}

#endif


int system_of_ODEs(double x, const double y[], double *dydx, void *param)
{
  
  double a1, b1, Esq, de_eos, de_Esq, de_a1;

  /*
    0 : cosmic time t
    then loop over k bins:
    1 : dD(1) / da
    2 : D(1)
    3 : dD(2) / da
    4 : D(2)
    5 : dD(3a) / da
    6 : D(3a)
    7 : dD(3b) / da
    8 : D(3b)
  */

  if (params.simpleLambda) 
    {
      de_Esq = 1.;
      de_a1  = 0.;
    }
  else
    {
      de_eos = my_spline_eval(SPLINE[SP_INTEOS], log10(x), ACCEL[SP_INTEOS]);
      de_Esq = exp(3.*de_eos)/pow(x,3.0);
      de_a1 = - 3*((1+DE_EquationOfState(x))/pow(x,4)) * params.OmegaLambda*exp(3*de_eos);
    }

  Esq = params.Omega0/pow(x,3.0)    /* adimensional Hubble factor squared */
    + OmegaK/(x*x)
    + OmegaRad/pow(x,4.0)
    + params.OmegaLambda * de_Esq;

  /* coefficients in the growth rate equations */
  a1   = -(3./x + 0.5 *((- 3.*params.Omega0/pow(x,4.0)
			 - 2*OmegaK/pow(x,3.0)
			 - 4*OmegaRad/pow(x,5.0)
			 + de_a1)/Esq));
  b1   = 1.5 *params.Omega0/(Esq*pow(x,5.0));


  /* cosmic time */
  dydx[0] = 1.0/x/sqrt(Esq);

#ifdef FOMEGA_GAMMA
  /* this is a toy model: D(t) is the one obtained by forcing 
     the gamma RSD parameter to a certain value. 
     Higher orders are obtained with Bouchet's fits */

  dydx[1] = 0.;
  dydx[2] =  (pow(params.Omega0/pow(x,3.0)/Esq,FOMEGA_GAMMA))/x * y[1];

  for (i=3; i<9; i++)
    dydx[i]=0.;

  return GSL_SUCCESS;
#endif

#if !defined(SCALE_DEPENDENT) && !defined(MOD_GRAV_FR)

  dydx[1] = a1*y[1] + b1*y[2];                                           /* d^2 D1/ dt^2 */
  dydx[2] = y[1];                                                        /* d D1/ dt */
  dydx[3] = a1*y[3] + b1*y[4] - b1*y[2]*y[2];                            /* d^2 D2/ dt^2 */
  dydx[4] = y[3];                                                        /* d D2/ dt */
  dydx[5] = a1*y[5] + b1*y[6] - 2.*b1*y[2]*y[2]*y[2];                    /* d^2 D31/ dt^2 */
  dydx[6] = y[5];                                                        /* d D31/ dt */
  dydx[7] = a1*y[7] + b1*y[8] - 2.*b1*y[2]*y[4] + 2.*b1*y[2]*y[2]*y[2];  /* d^2 D32/ dt^2 */
  dydx[8] = y[7];                                                        /* d D32/ dt */

  return GSL_SUCCESS;

#endif

#ifdef MOD_GRAV_FR

  /* equations for D2(k,a) as in Moretti et al. (2019) */

  double B1,B2,kkk,PI1,PI2,M2;
  B1  = params.Omega0 + 4. * params.OmegaLambda;
  B2  = params.Omega0 / pow(x, 3.) + 4. * params.OmegaLambda;

  for (int ik=0; ik<NkBINS; ik++)
    {
      if (!ik)
	kkk=0.0;
      else
	kkk = pow(10., LOGKMIN + ik*DELTALOGK);
      PI1 = kkk * kkk / x / x      + 0.5 * H_over_c * H_over_c * pow(B2, 3.) / (B1 * B1 * FR0);
      PI2 = kkk * kkk / x / x / 2. + 0.5 * H_over_c * H_over_c + pow(B2, 3.) / (B1 * B1 * FR0);
      M2  = params.Omega0 * H_over_c * H_over_c * kkk * kkk 
	* (1.5 * H_over_c / FR0) * (1.5 * H_over_c / FR0)
	* pow(B2, 5.) / pow(B1, 4.) / (9. * pow(x, 5.));
      
      dydx[1 + 8*ik] = a1 * y[1 + 8*ik] + mu(x, kkk) * b1 * y[2 + 8*ik];
      dydx[2 + 8*ik] = y[1 + 8*ik];
      dydx[3 + 8*ik] = a1 * y[3 + 8*ik] + mu(x, kkk) * b1 * y[4 + 8*ik]
	- (mu(x, kkk) - M2 / PI1 / PI2 / PI2) * b1 * y[2 + 8*ik] * y[2 + 8*ik] ; 
      dydx[4 + 8*ik] = y[3 + 8*ik];
      
      /* third-order growth is not used, it is set to the LCDM one */
      dydx[5 + 8*ik] = a1*y[5 + 8*ik] + b1*y[6 + 8*ik]
	- 2.*b1*y[2 + 8*ik]*y[2 + 8*ik]*y[2 + 8*ik];              /* d^2 D31/ dt^2 */
      dydx[6 + 8*ik] = y[5 + 8*ik];                               /* d D31/ dt */
      dydx[7 + 8*ik] = a1*y[7 + 8*ik] + b1*y[8 + 8*ik] 
	- 2.*b1*y[2 + 8*ik]*y[4 + 8*ik] 
	+ 2.*b1*y[2 + 8*ik]*y[2 + 8*ik]*y[2 + 8*ik];              /* d^2 D32/ dt^2 */
      dydx[8 + 8*ik] = y[7 + 8*ik];                               /* d D32/ dt */

    }

  return GSL_SUCCESS;

#endif

  return GSL_FAILURE;
}


int jac( double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  printf("This integration method should not call this function.\n");
  return GSL_FAILURE;
}

/* integrand for the computation of comoving distance */
double IntegrandComovingDistance(double z, void *param)
{
  return 1./Hubble(z);
}

/***********************************/
/* EQUATION OF STATE OF DE SECTION */
/***********************************/

int read_TabulatedEoS(void)
{
  /* Reads the tabulated Equation of State of dark energy from a file */

  int i;
  FILE *fd;
  double a, w;
  double *scalef,*EoS;

  if (!ThisTask)
    {
      if(!(fd = fopen(params.TabulatedEoSfile, "r")))
	{
	  printf("ERROR on task 0: can't open tabulated EoS in file '%s'\n", params.TabulatedEoSfile);
	  fflush(stdout);
	  return 1;
	}

      NtabEoS=0;
      while (1)
	{
	  if(fscanf(fd, " %lg %lg ", &a, &w) == 2)
	    NtabEoS++;
	  else
	    break;
	}

      fclose(fd);

      if (!NtabEoS)
	{
	  printf("ERROR on task 0: can't read tabulated EoS in file '%s'\n", params.TabulatedEoSfile);
	  fflush(stdout);
	  return 1;
	}

      printf("Found %d pairs of values in tabulated EoS file\n", NtabEoS);
      fflush(stdout);
    }

  /* Task 0 communicates the number of lines to all tasks */
  MPI_Bcast(&NtabEoS, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  SPLINE[SP_EOS] = gsl_spline_alloc(gsl_interp_cspline, NtabEoS);

  scalef = (double*)malloc(NtabEoS * sizeof(double));
  EoS    = (double*)malloc(NtabEoS * sizeof(double));

  if (!ThisTask)
    {
      fd = fopen(params.TabulatedEoSfile, "r");

      for (i=0; i<NtabEoS; i++)
	{
	  if(fscanf(fd, " %lg %lg ", &a, &w) == 2)
	    {
	      scalef[i]=log10(a);
	      EoS[i]=w;
	    }
	}

      fclose(fd);
    }

  /* Task 0 broadcasts the Pk and logk vectors to all tasks */
  MPI_Bcast(scalef, NtabEoS*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(EoS,    NtabEoS*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  gsl_spline_init(SPLINE[SP_EOS], scalef, EoS, NtabEoS);

  free(EoS);
  free(scalef);

  return 0;
}

/* parametric equation of state for Dark Energy */
double DE_EquationOfState(double a)
{
  if (!NtabEoS)
    return params.DEw0+(1-a)*params.DEwa;
  else
    return my_spline_eval(SPLINE[SP_EOS], log10(a), ACCEL[SP_EOS]);
}

double IntegrandForEoS(double a, void *param)
{
  return DE_EquationOfState(a)/a;
}



/**************************/
/* POWER SPECTRUM SECTION */
/**************************/

/* Part of this code is taken from N-GenIC (see GenIC.c) */
double PowerSpec_Tabulated(double);
/* double PowerSpec_Efstathiou(double); */
double PowerSpec_EH(double);
double transf_EH(double);
double T0(double,double,double);

//#define BODE_ET_AL_A9

double PowerSpectrum(double k)
{
  /* input: k (1/Mpc) */
  /* output: P(k) (Mpc^3) */

  double power, alpha, Tf;

  switch (WhichSpectrum)
    {
    case 1:
      power = PowerSpec_EH(k);
      break;

    case 2:
      power = PowerSpec_Tabulated(k);
      break;

    /* case 3: */
    /*   power = PowerSpec_Efstathiou(k); */
    /*   break; */

/* #ifdef SCALE_DEPENDENT */
/*     case 4: */
/*       power = PowerSpec_Tabulated(k); */
/*       break; */

/*     case 5: */
/*       power = PowerSpec_CAMBtable(k); */
/*       break; */
/* #endif */
      
    default:
      power = 0.0;
      break;

    }

  if(params.WDM_PartMass_in_kev > 0.)
    {
#ifdef BODE_ET_AL_A9
      /* Eqn. (A9) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
      /* NB: this is a length in Mpc/h*/
      alpha =
	0.048 * pow((params.Omega0 - params.OmegaBaryon) / 0.4, 0.15)
	* pow(params.Hubble100 / 0.65, 1.3) * pow(1.0 / params.WDM_PartMass_in_kev, 1.15);
      Tf = pow(1 + pow(alpha * k / params.Hubble100 * (3.085678e24 / UnitLength_in_cm), 2 * 1.2), -5.0 / 1.2);
#else
      /* Eqn. just after (A7) in Bode, Ostriker & Turok (2001), assuming gX=1.5  */
      /* NB: this is a length in Mpc/h*/
      alpha =
	0.05 * pow((params.Omega0 - params.OmegaBaryon) / 0.4, 0.15)
	* pow(params.Hubble100 / 0.65, 1.3) * pow(1.0 / params.WDM_PartMass_in_kev, 1.15);

      Tf = pow(1 + pow(alpha * k / params.Hubble100 * (3.085678e24 / UnitLength_in_cm), 2), -5.0);
#endif
      power *= Tf * Tf;
    }

  return PkNorm * power;
}

int initialize_PowerSpectrum(void)
{

//#ifndef SCALE_DEPENDENT
  if (!strcmp(params.FileWithInputSpectrum,"no"))
    WhichSpectrum=1;
  else
    WhichSpectrum=2;

//#else

  /* if (!strcmp(params.camb.MatterFile,"no")) */
  /*   WhichSpectrum=4; */
  /* else */
  /*   WhichSpectrum=5; */
//#endif

/*   double *logk,*Pk; */
/*   int i; */

/*   WhichSpectrum=2; */
/*   params.Sigma8=0.0; */

/*   /\* the interpolation is already loaded, the spline is initialized here  *\/ */
/*   NPowerTable=params.camb.Nkbins; */

/*   accelPk = gsl_interp_accel_alloc(); */
/*   splinePk = gsl_spline_alloc(gsl_interp_cspline, NPowerTable); */
/*   logk = (double*)malloc(NPowerTable * sizeof(double)); */
/*   Pk   = (double*)malloc(NPowerTable * sizeof(double)); */

/*   for (i=0; i<NPowerTable; i++) */
/*     { */
/*       logk[i] = log10(exp(params.camb.Logk[i])); */
/*       Pk  [i] = log10(exp(params.camb.LogPkref[i])) + 3.*logk[i]; */
/*     } */

/*   gsl_spline_init(splinePk, logk, Pk, NPowerTable); */

/*   free(Pk); */
/*   free(logk); */

/*   /\* it re-initializes the linear growing mode *\/ */
/*   accelCAMB = gsl_interp_accel_alloc(); */
/*   splineGrowCAMB = gsl_spline_alloc (gsl_interp_cspline, params.camb.NCAMB); */
/*   gsl_spline_init(splineGrowCAMB, params.camb.Scalef, params.camb.RefGM, params.camb.NCAMB); */
/*   splineInvGrowCAMB = gsl_spline_alloc (gsl_interp_cspline, params.camb.NCAMB); */
/*   gsl_spline_init(splineInvGrowCAMB, params.camb.RefGM, params.camb.Scalef, params.camb.NCAMB); */

/* #ifdef OUTPUT_GM */
/*   if (!ThisTask) */
/*     { */
/*       FILE *fd; */
/*       fd=fopen("TwoGrowingModes","w"); */
/*       fprintf(fd,"# 1: scale factor\n# 2: D as computed from equations\n# 3: from CAMB\n"); */
/*       double a=params.camb.Scalef[0]/1.05; */
/*       do */
/* 	{ */
/* 	  a*=1.05; */
/* 	  if (a>1) */
/* 	    a=1.0; */
/* 	  fprintf(fd," %g %g %g\n", */
/* 		  a, */
/* 		  pow(10.,my_spline_eval(splineGrow, log10(a), accel)), */
/* 		  exp(my_spline_eval(splineGrowCAMB, a, accelCAMB)) */
/* 		  ); */
/* 	} */
/*       while (a<1.0); */
/*       fclose(fd); */
/*     } */
/* #endif */

//#endif

  if (normalize_PowerSpectrum())
    return 1;

  return 0;
}


int normalize_PowerSpectrum(void)
{
  double tmp;

  WindowFunctionType=2;
  PkNorm=1.0;
  if (params.Sigma8!=0.0)
    {
      tmp = params.Sigma8 * params.Sigma8 / ComputeMassVariance(8.0/params.Hubble100);
      PkNorm = tmp;
      if (!ThisTask)
	printf("Normalization constant for the power spectrum: %g\n",PkNorm);
    }
  else
    {
      params.Sigma8 = sqrt(ComputeMassVariance(8.0/params.Hubble100));
#ifdef SCALE_DEPENDENT_GROWTH
      params.Sigma8 *= sqrt(params.camb.D2ref);
#endif
      if (!ThisTask)
	printf("Normalization of the provided P(k): Sigma8=%f\n",params.Sigma8);
    }

  return 0;
}


int read_Pk_from_file(void)
{
  /* This is adapted from N-GenIC */

  int i, oldread=0;
  FILE *fd;
  double k, p;
  double *logk,*Pk;

  if (!ThisTask)
    {
      if(!(fd = fopen(params.FileWithInputSpectrum, "r")))
	{
	  printf("ERROR on task 0: can't open input spectrum in file '%s' on task %d\n", 
		 params.FileWithInputSpectrum, ThisTask);
	  fflush(stdout);
	  return 1;
	}

      NPowerTable = 0;
      do
	{
	  if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
	    NPowerTable++;
	  else
	    break;
	}
      while(1);

      fclose(fd);

      if (!NPowerTable)
	{
	  printf("ERROR on task 0: can't read data from input spectrum in file '%s'\n", 
		 params.FileWithInputSpectrum);
	  fflush(stdout);
	  return 1;
	}

      printf("Found %d pairs of values in input spectrum table\n", NPowerTable);
      fflush(stdout);
    }

  /* Task 0 communicates the number of lines to all tasks */
  MPI_Bcast(&NPowerTable, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);

  SPLINE[SP_PK] = gsl_spline_alloc(gsl_interp_cspline, NPowerTable);

  logk = (double*)malloc(NPowerTable * sizeof(double));
  Pk   = (double*)malloc(NPowerTable * sizeof(double));

  if (!ThisTask)
    {
      fd = fopen(params.FileWithInputSpectrum, "r");

      for (i=0; i<NPowerTable; i++)
	{
	  if(fscanf(fd, " %lg %lg ", &k, &p) == 2)
	    {

	      if (!i)
		{
		  if (k<0.0)
		    {
		      oldread=1;
		      printf("I am assuming that file %s contains Log k - Log k^3 P(k)\n",params.FileWithInputSpectrum);
		    }
		  else
		    {
		      oldread=0;
		      printf("I am assuming that file %s contains k - P(k)\n",params.FileWithInputSpectrum);
		    }
		}

	      if (oldread)
		{
		  logk[i] = k;
		  Pk[i]   = p;
		}
	      else
		{
		  logk[i] = log10(k);
		  Pk[i]   = log10(p*k*k*k);
		}

	      /* translates into Htrue */
	      logk[i] += log10(params.Hubble100);
	      /* fixes the units if necessary */
	      if (params.InputSpectrum_UnitLength_in_cm!=0.0)
		logk[i] += log10(params.InputSpectrum_UnitLength_in_cm/UnitLength_in_cm);
	    }
	}

      fclose(fd);
    }

  /* Task 0 broadcasts the Pk and logk vectors to all tasks */
  MPI_Bcast(Pk,   NPowerTable*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(logk, NPowerTable*sizeof(double), MPI_BYTE, 0, MPI_COMM_WORLD);

  gsl_spline_init(SPLINE[SP_PK], logk, Pk, NPowerTable);

  free(Pk);
  free(logk);

  return 0;
}

double PowerSpec_Tabulated(double k)
{
  return pow(10.,my_spline_eval(SPLINE[SP_PK], log10(k), ACCEL[SP_PK]))/k/k/k;
}

/* #define ShapeGamma ((double)0.21) */
/* double PowerSpec_Efstathiou(double k) */
/* { */
/*   return pow(k,params.PrimordialIndex) / pow(1 + pow(6.4 / ShapeGamma * k + pow(3.0 / ShapeGamma * k, 1.5) + pow(1.7 / ShapeGamma,2.0) * k * k, 1.13), 2 / 1.13); */
/* } */

double PowerSpec_EH(double k)
{
  return pow(k,params.PrimordialIndex)*pow(transf_EH(k),2.);
}

double transf_EH(double fk)
{
  static double Teta_27=1.0104;
  double q,Omegac,Oh2,b1,b2,zd,Rd,zeq,Req,keq,s,ks,alc,bec,f,Tc,beb,bno,kst,ksi,Tb,Tr,y,alb,Ob2,OB;

  /* Eisenstein & Hu fit of the transfer function */

  OB   = (params.OmegaBaryon > 1.e-6 ? params.OmegaBaryon : 1.e-6);
  Omegac = params.Omega0-OB;
  Oh2 = params.Omega0*params.Hubble100*params.Hubble100;
  Ob2 = OB*params.Hubble100*params.Hubble100;
  b1  = 0.313*pow(Oh2,-0.419)*(1+0.607*pow(Oh2,0.674));
  b2  = 0.238*pow(Oh2,0.223);
  zd  = 1291.*pow(Oh2,0.251)*(1.+b1*pow(Ob2,b2))/(1.+0.659*pow(Oh2,0.828));
  Rd  = 31.5*Ob2/(pow(Teta_27,4.0)*0.001*zd);
  zeq = 2.5e4*Oh2/pow(Teta_27,4.0);
  Req = 31.5*Ob2/(pow(Teta_27,4.0)*0.001*zeq);
  keq = 7.46e-2*Oh2/Teta_27/Teta_27;
  s   = 1.633*log((sqrt(1.+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))/(keq*sqrt(Req)); /* 2sqrt(6)/3=1.633 */
  ks  = fk*s;
  q   = fk*Teta_27*Teta_27/Oh2;
  alc = pow(pow(46.9*Oh2,0.670)*(1.+pow(32.1*Oh2,-0.532)),-OB/params.Omega0) *
    pow(pow(12.0*Oh2,0.424)*(1.+pow(45.0*Oh2,-0.582)),-pow(OB/params.Omega0,3.0));
  bec = 1./(1.+(0.944/(1.+pow(458.*Oh2,-0.708)))*(pow(Omegac/params.Omega0,pow(0.395*Oh2,-0.0266))-1.));
  f   = 1./(1+pow(ks/5.4,4.0));
  Tc  = f*T0(q,1.,bec)+(1.-f)*T0(q,alc,bec);
  beb = 0.5+OB/params.Omega0+(3.-2.*OB/params.Omega0)*sqrt(pow(17.2*Oh2,2.0)+1.);
  bno = 8.41*pow(Oh2,0.435);
  kst = ks/pow(1.+pow(bno/ks,3.0),0.3333);
  ksi = 1.6*pow(Ob2,0.52)*pow(Oh2,0.73)*(1.+pow(10.4*Oh2,-0.95));
  y   = (1.+zeq)/(1+zd);
  alb = 2.07*keq*s*pow(1.0+Rd,-0.75)*(y*(-6.*sqrt(1.+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.))));
  Tb  = (T0(q,1.,1.)/(1.+pow(ks/5.2,2.0))+alb/(1.+pow(beb/ks,3.0))*exp(-pow(fk/ksi,1.4)))*sin(kst)/kst;
  Tr  = (OB*Tb+Omegac*Tc)/params.Omega0;

  return Tr;
}


double T0(double q,double a,double b)
{
  double ll,C;

  ll = log(exp(1.)+1.8*b*q);
  C  = 14.2/a+386./(1.+69.9*pow(q,1.08));

  return ll/(ll+C*q*q);
}

/*************************/
/* MASS VARIANCE SECTION */
/*************************/

double ThisRadius;
double IntegrandForMassVariance(double, void *);
double IntegrandForDisplVariance(double, void *);

int initialize_MassVariance(void)
{
  int i;
  double r;
  double rmin=-6.0, dr=0.04;   /* NB: rmax=2.0 */
  double *rv, *massvar, *dmvdr, *massvarneg, *displv;

  rv         = (double*)malloc(NBINS * sizeof(double));
  massvar    = (double*)malloc(NBINS * sizeof(double));
  dmvdr      = (double*)malloc(NBINS * sizeof(double));
  massvarneg = (double*)malloc(NBINS * sizeof(double));
  displv     = (double*)malloc(NBINS * sizeof(double));

  for (i=NBINS-1; i>=0; i--)
    {
      rv[i]=(rmin+i*dr);
      r=pow(10.,rv[i]);
      massvar[i]=log10(ComputeMassVariance(r));
      /* To allow interpolation of reverse function, 
	 the variance must increase with decreasing radius */
      if (i<NBINS-1 && massvar[i]-massvar[i+1]<1.e-6)
	massvar[i]=massvar[i+1]+1.e-6;

      massvarneg[i]=-massvar[i];
      displv[i]=log10(ComputeDisplVariance(r));
    }

  for (i=0; i<NBINS; i++)
    {
      if (i==0)
        dmvdr[i]=(massvar[i+1]-massvar[i])/(rv[i+1]-rv[i]);
      else if (i==NBINS-1)
        dmvdr[i]=(massvar[i]-massvar[i-1])/(rv[i]-rv[i-1]);
      else
        dmvdr[i]=(massvar[i+1]-massvar[i-1])/(rv[i+1]-rv[i-1]);
    }

  /* initialization of splines for interpolation */
  gsl_spline_init(SPLINE[SP_MASSVAR], rv, massvar, NBINS);
  gsl_spline_init(SPLINE[SP_RADIUS],  massvarneg, rv, NBINS);
  gsl_spline_init(SPLINE[SP_DVARDR],  rv, dmvdr, NBINS);
  gsl_spline_init(SPLINE[SP_DISPVAR], rv, displv, NBINS);

  free(displv    );
  free(massvarneg);
  free(dmvdr     );
  free(massvar   );
  free(rv        );

  return 0;
}

double ComputeMassVariance(double R)
{
  double result, error;

  ThisRadius=R;
  Function.function = &IntegrandForMassVariance;
  gsl_integration_qags (&Function, -10., log(500.0/R) , 0.0, TOLERANCE, NWINT, workspace, &result, &error);

  return result;
}

double IntegrandForMassVariance(double logk, void *param)
{
  double w, k;
  k = exp(logk);
  w = WindowFunction(k * ThisRadius);
  return PowerSpectrum(k) * w*w * k*k*k / (2.*PI*PI);
}

double ComputeDisplVariance(double R)
{
  double result, error;

  ThisRadius=R;
  Function.function = &IntegrandForDisplVariance;
  gsl_integration_qags (&Function, -10., log(500.0/R) , 0.0, TOLERANCE, NWINT, workspace, &result, &error);

  return result;
}

double IntegrandForDisplVariance(double logk, void *param) {
  double w, k;
  k = exp(logk);
  w = WindowFunction(k * ThisRadius);
  return PowerSpectrum(k) * w*w * k / (2.*PI*PI);
}

double WindowFunction(double kr) {
  /* Window function:

     WindowFunctionType = 0: Gaussian smoothing
     WindowFunctionType = 1: SKS smoothing
     WindowFunctionType = 2: top-hat smoothing

     DIMENSIONLESS
  */
  double window, kr2;

  switch (WindowFunctionType)
    {
    case 0:                       /* Gaussian */
      window=exp(-kr*kr/2.);
      break;

    case 1:                       /* sharp k-space */
      window=(kr<1? 1.0 : 0.0);
      break;

    case 2:                       /* top-hat */
      if (kr < 1.e-5)
	window=1.0;
      else
	{
	  kr2=kr*kr;
	  window=3.*(sin(kr)/kr2/kr-cos(kr)/kr2);
	}
      break;

    default:
      window=1;
      break;
    }

  return window;
}

double MassVariance(double R)
{
  return pow(10.,my_spline_eval(SPLINE[SP_MASSVAR], log10(R), ACCEL[SP_MASSVAR]));
}


double dMassVariance_dr(double R)
{
  return my_spline_eval(SPLINE[SP_DVARDR], log10(R), ACCEL[SP_DVARDR]);
}

double DisplVariance(double R)
{
  return pow(10.,my_spline_eval(SPLINE[SP_DISPVAR], log10(R), ACCEL[SP_DISPVAR]));
}

double Radius(double Var)
{
  return pow(10.,my_spline_eval(SPLINE[SP_RADIUS], -log10(Var), ACCEL[SP_RADIUS]));
}


/**********************************/
/* COSMOLOGICAL FUNCTIONS SECTION */
/**********************************/

double OmegaMatter(double z)
{
  /* Cosmological mass density parameter as a function of redshift
     DIMENSIONLESS */

  return params.Omega0 * pow(1.+z,3.) * pow(Hubble(z)/100./params.Hubble100, -2);
}

double OmegaLambda(double z)
{
  /* Cosmological mass density parameter as a function of redshift
     DIMENSIONLESS */

  return params.OmegaLambda * pow(Hubble(z)/100./params.Hubble100, -2);
}


double Hubble(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: km/s/Mpc  */
  double Esq,de_eos;

  if (params.simpleLambda)
    Esq = params.Omega0*pow(1.+z,3.)+OmegaK*pow(1.+z,2.)+params.OmegaLambda;
  else
    {
      de_eos=my_spline_eval(SPLINE[SP_INTEOS], -log10(1.+z), ACCEL[SP_INTEOS]);
      Esq = params.Omega0*pow(1.+z,3.)+OmegaK*pow(1.+z,2.)+params.OmegaLambda*pow(1.+z,3.) * exp(3.*de_eos);
    }

  return 100.*params.Hubble100*sqrt(Esq);
}


double Hubble_Gyr(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: Gyr^-1 */

  return Hubble(z)/HUBBLETIME_GYR/100.;
}


double InterpolateGrowth(double z, double k, int pointer)
{
  /* This function interpolates the table for all scale-dependent growth functions */

#ifdef SCALE_DEPENDENT
  int kk;
  double dk;
  if (k<kmin)
    return my_spline_eval(SPLINE[pointer], -log10(1.+z), ACCEL[pointer]);
  else if (k>kmax)
    return my_spline_eval(SPLINE[pointer+NkBINS-1], -log10(1.+z), ACCEL[pointer+NkBINS-1]);
  else
    {
      dk = (log10(k)-LOGKMIN)/DELTALOGK;
      kk=(int)dk;
      dk-=kk;

      return dk * my_spline_eval(SPLINE[pointer+kk+1], -log10(1.+z), ACCEL[pointer+kk+1]) +
	 (1-dk) * my_spline_eval(SPLINE[pointer+kk  ], -log10(1.+z), ACCEL[pointer+kk  ]);
    }
#else
  /* scale-independent case, just return the interpolation */

  return my_spline_eval(SPLINE[pointer], -log10(1.+z), ACCEL[pointer]);
#endif

}


double fomega(double z, double k)
{
  /* Peebles' f(Omega) function, dlogD/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z,k,SP_FOMEGA1);
}


double fomega_2LPT(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z,k,SP_FOMEGA2);
}


double fomega_3LPT_1(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z,k,SP_FOMEGA31);
}

double fomega_3LPT_2(double z, double k)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return InterpolateGrowth(z,k,SP_FOMEGA32);
}

double GrowingMode(double z, double k)
{
  /* linear growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10.,InterpolateGrowth(z,k,SP_GROW1));
}


double GrowingMode_2LPT(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10.,InterpolateGrowth(z,k,SP_GROW2));
}

double GrowingMode_3LPT_1(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return -pow(10.,InterpolateGrowth(z,k,SP_GROW31));
}

double GrowingMode_3LPT_2(double z, double k)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10.,InterpolateGrowth(z,k,SP_GROW32));
}


#ifdef ELL_CLASSIC
double InverseGrowingMode(double D)
{
  /* redshift corresponding to a linear growing mode, interpolation on the grid
     DIMENSIONLESS */

  return 1./pow(10.,my_spline_eval(SPLINE[SP_INVGROW], log10(D), ACCEL[SP_INVGROW])) -1.;
}
#endif

double CosmicTime(double z)
{
  /* cosmic time, interpolation on the grid
     Gyr */

  return pow(10.,my_spline_eval(SPLINE[SP_TIME], -log10(1.+z), ACCEL[SP_TIME]));
}

double InverseCosmicTime(double t)
{
  /* scale factor corresponding to a cosmic time, interpolation on the grid
     Gyr */

  return pow(10.,my_spline_eval(SPLINE[SP_INVTIME], log10(t), ACCEL[SP_INVTIME]));
}

double ComovingDistance(double z)
{
  /* comoving distance, interpolation on the grid
     Mpc */

  return my_spline_eval(SPLINE[SP_COMVDIST], -log10(1.+z), ACCEL[SP_COMVDIST]);
}

double DiameterDistance(double z)
{
  /* diameter distance, interpolation on the grid
     Mpc */

  return my_spline_eval(SPLINE[SP_DIAMDIST], -log10(1.+z), ACCEL[SP_DIAMDIST]);
}

double SizeForMass(double m)
{
  /* Radius corresponding to mass m
     m: M_sun 

     DIMENSION: Mpc
  */

  switch(WindowFunctionType)
    {
    case 0:
      return pow(m/pow(2.0*PI,1.5)/MatterDensity,0.3333333333333333);
      break;
    case 1:
      return pow(m/(6.0*PI*PI*MatterDensity),0.3333333333333333);
      break;
    case 2:
      return pow(m/(4.0*PI*MatterDensity/3.0),0.3333333333333333);
      break;
    default:
      return 0.;
      break;
    }
}

double MassForSize(double size)
{

  switch(WindowFunctionType)
    {
    case 0:
      return MatterDensity * pow(2.0*PI,1.5)*pow(size,3.0);
      break;
    case 1:
      return MatterDensity * 6.*PI*PI*pow(size,3.0);
      break;
    case 2:
      return MatterDensity * 4.*PI/3.*pow(size,3.0);
      break;
    default:
      return 0;
      break;
    }
}

/* double TypicalCollapsingMass(double z) */
/* { */
/*   /\* NB: this is never used...  */
/*      for SCALE_INDEPENDENT_GROWTH one should specify the scale *\/ */
/*   return MassForSize(Radius(pow(DELTA_C/GrowingMode(z),2.0))); */
/* } */

/**********************************/
/* ANALYTIC MASS FUNCTION SECTION */
/**********************************/

#define SQRT2PI ((double)0.39894228)
#define ALPHA ((double)0.569558118758974)

double dOmega_dVariance(double v,double z)
{

  /*
    dOmega/dLambda

    DIMENSIONLESS

    params.AnalyticMassFunction = 0:  Press  & Schechter (1974)
    params.AnalyticMassFunction = 1:  Sheth & Tormen (2001)
    params.AnalyticMassFunction = 2:  Jenkins et al. (2001)
    params.AnalyticMassFunction = 3:  Warren et al. (2006)
    params.AnalyticMassFunction = 4:  Reed et al. (2007)
    params.AnalyticMassFunction = 5:  Crocce et al. (2010)
    params.AnalyticMassFunction = 6:  Tinker et al. (2008)
    params.AnalyticMassFunction = 7:  Courtin et al. (2010)
    params.AnalyticMassFunction = 8:  Angulo et al. (2012)
    params.AnalyticMassFunction = 9:  Watson et al. (2013)
    params.AnalyticMassFunction =10:  Crocce et al. (2010), with forced universality

  */

  double sv, ni, ni2, onepz;

  sv=sqrt(v);
  ni=DELTA_C/sv;

  switch (params.AnalyticMassFunction)
    {
    case 0:       // Press & Schechter
      return 2.* exp(-0.5*ni*ni) * ni * SQRT2PI;
      break;

    case 1:       // Sheth & Tormen
      ni2=sqrt(0.707)*ni;
      return 2.* 0.3222 * SQRT2PI * ni2 * exp(-0.5*ni2*ni2) * (1.0+ 1.0/pow(ni2,0.6));
      break;

    case 2:       // Jenkins et al.
      return 0.315 * exp(-pow(fabs(-log(sv)+0.61),3.8));
      break;

    case 3:       // Warren et al. (2006)
      return 0.7234 * (pow(sv,-1.625)+0.2538)*exp(-1.1982/v);
      break;

    case 4:       // Reed et al. (2007)
      ni2=sqrt(0.707)*ni;
      return 2.*0.3222 * SQRT2PI * ni2 * exp(-0.54*ni2*ni2)*
	(1.0+ 1.0/pow(ni2,0.6) + 0.2 * exp(-(pow(-log(sv)-0.4,2.0)/0.72)));
      break;

    case 5:       // Crocce et al. (2010)
      onepz = (z<1 ? 1.0+z : 2.0);
      return 0.58* pow(onepz,-0.13)
	* (pow( sv,-1.37 * pow(onepz,-0.15) ) + 0.3 * pow(onepz,-0.084))
	* exp(-1.036 * pow(onepz,-0.024)/v);
      break;

    case 6:       // Tinker et al. (2010)
      onepz = (z<2.5 ? 1.0+z : 3.5);
      return 0.186*pow(onepz,-0.14)*(pow(2.57*pow(onepz,-0.569558118758974)/sv,
					 1.47*pow(onepz,-0.06))+1.)*exp(-1.19/v);
      break;

    case 7:       // Courtin et al. (2010)
      ni2=sqrt(0.695)*1.673/sv;
      return 0.348 * 2.*SQRT2PI * ni2 * (1.+pow(1./ni2/ni2,0.1))*exp(-ni2*ni2/2.);
      break;

    case 8:       // Angulo et al. (2012)
      return 0.201*(pow(ni*2.08/DELTA_C,1.7) + 1.0) *exp(-1.172 *ni*ni/DELTA_C/DELTA_C);
      break;

    case 9:       // Watson et al. (2013)
      return 0.282*(pow(ni*1.406/DELTA_C,2.163) + 1.0) *exp(-1.210 *ni*ni/DELTA_C/DELTA_C);
      break;

    case 10:       // Crocce et al. (2010) universal
      onepz = 1.0;
      return 0.58* pow(onepz,-0.13)
	* (pow( sv,-1.37 * pow(onepz,-0.15) ) + 0.3 * pow(onepz,-0.084))
	* exp(-1.036 * pow(onepz,-0.024)/v);
      break;

    default:
      return 0.0;
      break;

    }

}

double AnalyticMassFunction(double mass, double z)
{
  double r,D;

  r=SizeForMass(mass);
#ifdef SCALE_DEPENDENT_GROWTH
  /* This function still must be adapted to scale-dependent growing mode */
  SDGM.flag=0;       /* mass growing mode */
  SDGM.radius=r;
#endif
  D=GrowingMode(z,0.);

  return MatterDensity * dOmega_dVariance(MassVariance(r)*D*D,z)
    * fabs(dMassVariance_dr(r) /6.0) /mass /mass;

}

double my_spline_eval(gsl_spline *spline, double x, gsl_interp_accel *accel)
{
  /* this function performs linear extrapolation beyond the x-range limits,
     and calls the spline evaluation in between */
  if (x<spline->x[0])
    return spline->y[0]+(x-spline->x[0])*(spline->y[1]-spline->y[0])/(spline->x[1]-spline->x[0]);
  else if (x>spline->x[spline->size-1])
    return spline->y[spline->size-1]+(x-spline->x[spline->size-1])*
      (spline->y[spline->size-1]-spline->y[spline->size-2])/(spline->x[spline->size-1]-spline->x[spline->size-2]);
  else
    return gsl_spline_eval(spline,x,accel);

}

