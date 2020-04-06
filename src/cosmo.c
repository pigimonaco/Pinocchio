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

#define NVAR 5
#define NBIN 210
#ifdef NORADIATION
#define OMEGARAD_H2 ((double)0.0)
#else
#define OMEGARAD_H2 ((double)4.2e-5)
#endif
#define TOLERANCE ((double)1.e-4)
#define UnitLength_in_cm ((double)3.085678e24)
#define ShapeGamma ((double)0.21)
#define HUBBLETIME_GYR ((double)3.085678e24/(double)1.e7/(double)3.1558150e16)
#define DELTA_C ((double)1.686)
//#define FOMEGA_GAMMA 0.554
//#define OUTPUT_GM

static int WhichSpectrum, NPowerTable=0, NtabEoS=0;;
static double PkNorm, MatterDensity, OmegaK, OmegaRad;

/* declaration of gsl quantities */
gsl_function Function;

static gsl_spline *splineGrow=0x0, *splineGrow2=0x0,*splineGrow31=0x0,*splineGrow32=0x0, *splineTime=0x0, 
  *splineEoS=0x0, *splinePk=0x0, *splineVar=0x0, *splineInvGrow=0x0, 
  *splineInvTime=0x0, *splineDist=0x0, *splineInvDist=0x0, *splinedrdz=0x0, 
  *splineRadius=0x0, *splinefO=0x0, *splinefO2=0x0,*splinefO31=0x0,*splinefO32=0x0, *splinedldr=0x0, 
  *splineIntEoS=0x0, *splineDVar=0x0;
static gsl_interp_accel *accel=0x0, *accelPk=0x0, *accelEoS=0x0;

#ifdef SCALE_DEPENDENT_GROWTH
static gsl_spline *splineGrowCAMB=0x0, *splineInvGrowCAMB=0x0;
static gsl_interp_accel *accelCAMB=0x0;
#endif

int system_of_ODEs(double, const double [], double*, void* );
int jac(double, const double [], double *, double [], void *);
int read_TabulatedEoS(void);
int initialize_PowerSpectrum(void);
int normalize_PowerSpectrum(void);
double IntegrandForEoS(double, void*);
double DE_EquationOfState(double);
double IntegrandProperDistance(double, void*);
double ComputeMassVariance(double);
double ComputeDisplVariance(double);


/**************************/
/* INITIALIZATION SECTION */
/**************************/

int initialize_cosmology()
{

  /*
    Computes the following functions:
    Scale factor, growth first and second order LPT, cosmic time, proper distance
    on a grid of values to be interpolated
  */

  double ode_param;
  double y[NVAR], x1, x2, hh, h1, norm, k, result, error, SqrtOmegaK;
  int status=GSL_SUCCESS, i, i1, i2;
  char filename[BLENGTH];
  FILE *fd;
  double log_amin=-4.,log_amax=0.2, dloga;

  double *scalef, *cosmtime, *grow1, *grow2, *IntEoS, *prdist, *prdistneg, *dprdist,
    *fomegav, *fomega2v, *grow31, *grow32, *fomega31v, *fomega32v;

  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step *ode_s = gsl_odeiv_step_alloc(T,NVAR);
  gsl_odeiv_control *ode_c = gsl_odeiv_control_standard_new(1.0e-8, 1.0e-8, 1.0, 1.0);
  gsl_odeiv_evolve *ode_e = gsl_odeiv_evolve_alloc(NVAR);
  gsl_odeiv_system ode_sys = {system_of_ODEs, jac, NVAR, (void*)&ode_param};

  gsl_spline *spline=0x0;

  dloga=(log_amax-log_amin)/(double)NBIN;

  OmegaRad = OMEGARAD_H2 / params.Hubble100 / params.Hubble100;
  OmegaK = 1.0-params.Omega0 - params.OmegaLambda - OmegaRad;
  SqrtOmegaK = sqrt(fabs(OmegaK));
  MatterDensity = 2.775499745e11 * params.Hubble100 * params.Hubble100 * params.Omega0;
  if (params.DEw0==-1 && params.DEwa==0 && strcmp(params.TabulatedEoSfile,"no"))
    params.simpleLambda=1;
  else
    params.simpleLambda=0;

  /* allocation of splines */
  splineGrow    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineInvGrow = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineGrow2   = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineGrow31  = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineGrow32  = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinefO      = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinefO2     = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinefO31    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinefO32    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  if (!params.simpleLambda) 
    splineIntEoS  = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineTime    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineInvTime = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineDist    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineInvDist = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinedrdz    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineVar     = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineDVar    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splineRadius  = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  splinedldr    = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  accel         = gsl_interp_accel_alloc ();

  /* read the tabulated Equation of State of the dark energy */
  if (strcmp(params.TabulatedEoSfile,"no"))
    if (read_TabulatedEoS())
      return 1;

  /* allocation of vectors for interpolation */
  scalef     = (double*)malloc(NBIN * sizeof(double));
  cosmtime   = (double*)malloc(NBIN * sizeof(double));
  grow1      = (double*)malloc(NBIN * sizeof(double));
  grow2      = (double*)malloc(NBIN * sizeof(double));
  grow31     = (double*)malloc(NBIN * sizeof(double));
  grow32     = (double*)malloc(NBIN * sizeof(double));
  prdist     = (double*)malloc(NBIN * sizeof(double));
  prdistneg  = (double*)malloc(NBIN * sizeof(double));
  dprdist    = (double*)malloc(NBIN * sizeof(double));
  fomegav    = (double*)malloc(NBIN * sizeof(double));
  fomega2v   = (double*)malloc(NBIN * sizeof(double));
  fomega31v  = (double*)malloc(NBIN * sizeof(double));
  fomega32v  = (double*)malloc(NBIN * sizeof(double));
  if (!params.simpleLambda)
    IntEoS  = (double*)malloc(NBIN * sizeof(double));


  /* compute the integrand for the DE EoS if necessary */
  if (!params.simpleLambda)
    {
      Function.function = &IntegrandForEoS;
      for (i=0; i<NBIN; i++)
	{
	  x2=pow(10., log_amin + (i+1)*dloga);
	  gsl_integration_qags(&Function, x2, 1.0, 0.0, TOLERANCE, NWINT, workspace, &result, &error);
	  scalef[i] = log_amin + (i+1)*dloga;
	  IntEoS[i] = result;
	}
      gsl_spline_init(splineIntEoS, scalef, IntEoS, NBIN);
    }

  /* ICs for the runge-kutta integration */
  x1=1.e-6;                   /*  initial value of scale factor a */
  y[0]=1.0;                   /*  initial value of dD/da */
  y[1]=x1;                    /*  initial value of D(a) */
  y[2]=2.0/3.0*pow(x1,1.5);   /*  initial value of t(a)*Hubble0 */
  y[3]=-6./7.*x1;             /*  initial value of the dD2/da */
  y[4]=-3./7.*x1*x1;          /*  intial value of D2(a) */
  h1=1.e-6;                   /*  initial guess of time-step */
  hh=h1;

  /* this function will be integrated within the loop */
  Function.function = &IntegrandProperDistance;

  /* loop on bins for calling the ODE integrator */
  for (i=0; i<NBIN; i++)
    {
      x2=pow(10., log_amin + i*dloga);

      /* integration of ODE system */
      while (x1<x2 && status==GSL_SUCCESS)
	{
	  status=gsl_odeiv_evolve_apply(ode_e, ode_c, ode_s, &ode_sys, &x1, x2, &hh, y);
	  if (status!=GSL_SUCCESS)
	    {
	      printf("ERROR on task %d: integration of cosmological quantities failed\n",ThisTask);
	      fflush(stdout);
	      return 1;
	    }
	}

      scalef[i]   = x2;
      cosmtime[i] = y[2]*HUBBLETIME_GYR/params.Hubble100;
      grow1[i]    = y[1];
      grow2[i]    = -y[4];

      /* calculation of proper distance (Mpc) for a generic cosmology (i.e. flat, open or closed) 
	 and generic equation of state for the DE component  */

      gsl_integration_qags(&Function, 0.0, 1./x2-1., 0.0, TOLERANCE, NWINT, workspace, &result, &error);

      if (fabs(OmegaK) < 1.e-4)
	prdist[i] = SPEEDOFLIGHT*result/(params.Hubble100*100.);

      else if(OmegaK <0)
	prdist[i] = SPEEDOFLIGHT/(params.Hubble100*100.)/SqrtOmegaK*sin(SqrtOmegaK*result);

      else
        prdist[i]=SPEEDOFLIGHT/(params.Hubble100*100.)/SqrtOmegaK*sinh(SqrtOmegaK*result);

      x1=x2;

    }

  /* normalization of the first- and second- order growth rate */
  spline = gsl_spline_alloc (gsl_interp_cspline, NBIN);
  gsl_spline_init(spline, scalef, grow1, NBIN);
  norm=my_spline_eval(spline, 1.00, accel);
  for (i=0; i<NBIN; i++)
    {
      grow1[i] = grow1[i] / norm;
      grow2[i] = grow2[i] / norm/norm;
      grow31[i] = grow1[i]*grow1[i]*grow1[i]*pow(params.Omega0,-4./275.)/9.; //the "-" will be set later
      grow32[i] = grow1[i]*grow1[i]*grow1[i]*pow(params.Omega0,-269./17875.)*5./42.;
    }
  gsl_spline_free(spline);

  /* f(Omega) at first and second order */
  for (i=0; i<NBIN; i++)
    {
      i1=(i>0 ? i-1 : 0);
      i2=(i<NBIN-1 ? i+1 : NBIN-1);
      fomegav[i] = (grow1[i2]-grow1[i1]) / (scalef[i2]-scalef[i1]) * scalef[i] / grow1[i];
      fomega2v[i]= (grow2[i2]-grow2[i1]) / (scalef[i2]-scalef[i1]) * scalef[i] / grow2[i];
      fomega31v[i] = (grow31[i2]-grow31[i1]) / (scalef[i2]-scalef[i1]) * scalef[i] / grow31[i];
      fomega32v[i]= (grow32[i2]-grow32[i1]) / (scalef[i2]-scalef[i1]) * scalef[i] / grow32[i];
    }

  /* derivative of proper distance */
  for (i=0; i<NBIN; i++)
    {
      prdistneg[i]=-prdist[i];
      i1=(i>0 ? i-1 : 0);
      i2=(i<NBIN-1 ? i+1 : NBIN-1);
      dprdist[i]=(prdist[i2]-prdist[i1])/(1./scalef[i2]-1./scalef[i1]);
    }

  /* these quantities will be interpolated logarithmically */
  for (i=0; i<NBIN; i++)
    {
      scalef[i]   = log10(scalef[i]  );
      cosmtime[i] = log10(cosmtime[i]);
      grow1[i]    = log10(grow1[i]   );
      grow2[i]    = log10(grow2[i]   );
      grow31[i]    = log10(grow31[i]   );
      grow32[i]    = log10(grow32[i]   );
    }

  /* initialization of spline interpolations of time-dependent quantities */
  gsl_spline_init(splineGrow, scalef, grow1, NBIN);
  gsl_spline_init(splineInvGrow, grow1, scalef, NBIN);
  gsl_spline_init(splineGrow2, scalef, grow2, NBIN);
  gsl_spline_init(splineGrow31, scalef, grow31, NBIN);
  gsl_spline_init(splineGrow32, scalef, grow32, NBIN);
  gsl_spline_init(splinefO, scalef, fomegav, NBIN);
  gsl_spline_init(splinefO2, scalef, fomega2v, NBIN);
  gsl_spline_init(splinefO31, scalef, fomega31v, NBIN);
  gsl_spline_init(splinefO32, scalef, fomega32v, NBIN);
  gsl_spline_init(splineTime, scalef, cosmtime, NBIN);
  gsl_spline_init(splineInvTime, cosmtime, scalef, NBIN);
  gsl_spline_init(splineDist, scalef, prdist, NBIN);
  gsl_spline_init(splineInvDist, prdistneg, scalef, NBIN);
  gsl_spline_init(splinedrdz, scalef, dprdist, NBIN);


  /* deallocation of vectors for interpolation */
  if (!params.simpleLambda)
    free(IntEoS  );
  free(fomega2v);
  free(fomegav );
  free(dprdist );
  free(prdistneg);
  free(prdist  );
  free(grow2   );
  free(grow1   );
  free(grow32  );
  free(grow31  );
  free(cosmtime);
  free(scalef  );

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
      fprintf(fd,"# 3: growth factor\n");
      fprintf(fd,"# 4: 2nd-order growth factor\n");
      fprintf(fd,"# 5: dark energy EOS\n");
      fprintf(fd,"# 6: proper distance (true Mpc)\n");
      fprintf(fd,"# 7: d/dz of proper distance (true Mpc)\n");
      fprintf(fd,"# SCALE-DEPENDENT QUANTITIES\n");
      fprintf(fd,"# 8: scale (true Mpc)\n");
      fprintf(fd,"# 9: mass variance\n");
      fprintf(fd,"#10: d Log sigma^2 / d Log R\n");
      fprintf(fd,"#11: k (true Mpc^-1)\n");
      fprintf(fd,"#12: P(k)\n");
      fprintf(fd,"#\n");

      for (i=0; i<NBIN; i++)
	{
	  k=pow(10.,-4.0+(double)i/(double)NBIN*6.0);
	  fprintf(fd," %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",
		  pow(10.,splineGrow->x[i]),
		  pow(10.,splineTime->y[i]),
		  pow(10.,splineGrow->y[i]),
		  -pow(10.,splineGrow2->y[i]),
		  (!params.simpleLambda ? -1 : splineEoS->y[i]),
		  splineDist->y[i],
		  splinedrdz->y[i],
		  pow(10.,splineVar->x[i]),
		  pow(10.,splineVar->y[i]),
		  pow(10.,splinedldr->y[i]),
		  k,PowerSpectrum(k));
	  }
      fclose(fd);
    }  

  return 0;
}


int system_of_ODEs(double x, const double y[], double *dydx, void *param)
{
  
  double a1, b1, Esq, eqds;

  /* NB: qui si potrebbe compattare tutto in un'unico caso */


  if (params.simpleLambda) 
    {
      Esq = params.Omega0/pow(x,3.0)
        + OmegaK/(x*x)
        + params.OmegaLambda
	+ OmegaRad/pow(x,4.0);
      a1   = -(3./x + 1./2.*((- 3.*params.Omega0/pow(x,4.0)
			      - 2*OmegaK/pow(x,3.0)
			      - 4*OmegaRad/pow(x,5.0))/Esq));
      b1   = 3./2.*params.Omega0/(Esq*pow(x,5.0));
    }
  else
    /* in this case Dark Energy is not a simple cosmological constant */
    {
      eqds=my_spline_eval(splineIntEoS, log10(x), accel);

      Esq = params.Omega0/pow(x,3.0)
	+ OmegaK/(x*x)
	+ params.OmegaLambda/pow(x,3.0) * exp(3.*eqds)
	+ OmegaRad/pow(x,4.0);
      a1   = -(3./x+1./2.*((- 3.*params.Omega0/pow(x,4.0)
			    - 2*OmegaK/pow(x,3.0)
			    - 4*OmegaRad/pow(x,5.0)
			    - 3*((1+DE_EquationOfState(x))/pow(x,4))
			    * params.OmegaLambda*exp(3*eqds))/Esq));
      b1   = 3./2.*params.Omega0/(Esq*pow(x,5.0));
    }

#ifndef FOMEGA_GAMMA
  dydx[0] = a1*y[0]+b1*y[1];
  dydx[1] = y[0];
#else
  dydx[0] = 0.;
  dydx[1] =  (pow(params.Omega0/pow(x,3.0)/Esq,FOMEGA_GAMMA))/x * y[1];
#endif
  dydx[2] = 1.0/x/sqrt(Esq);
  dydx[3] = a1*y[3]+b1*y[4]-b1*y[1]*y[1];
  dydx[4] = y[3];

  return GSL_SUCCESS;
}


int jac( double t, const double y[], double *dfdy, double dfdt[], void *params)
{
  printf("This integration method should not call this function.\n");
  return GSL_FAILURE;
}

/* integrand for the computation of proper distance */
double IntegrandProperDistance(double z, void *param)
{
  return 1./Hubble(z)*params.Hubble100*100.;
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

  accelEoS = gsl_interp_accel_alloc();
  splineEoS = gsl_spline_alloc(gsl_interp_cspline, NtabEoS);

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

  gsl_spline_init(splineEoS, scalef, EoS, NtabEoS);

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
    return my_spline_eval(splineEoS, log10(a), accelEoS);
}

double IntegrandForEoS(double a, void *param)
{
  return DE_EquationOfState(a)/a;
}



/**************************/
/* POWER SPECTRUM SECTION */
/**************************/

/* Part of this code is taken from N-GenIC (see GenIC.c) */
int read_Pk_table(void);
double PowerSpec_Tabulated(double);
double PowerSpec_Efstathiou(double);
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

    default:
      power = PowerSpec_Efstathiou(k);
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

#ifndef SCALE_DEPENDENT_GROWTH
  if (!strcmp(params.FileWithInputSpectrum,"no"))
    WhichSpectrum=1;
  else
    WhichSpectrum=2;

  if (strcmp(params.FileWithInputSpectrum,"no"))
    if (read_Pk_table())
      return 1;
#else
  double *logk,*Pk;
  int i;

  WhichSpectrum=2;
  params.Sigma8=0.0;

  /* the interpolation is already loaded, the spline is initialized here  */
  NPowerTable=params.camb.Nkbins;

  accelPk = gsl_interp_accel_alloc();
  splinePk = gsl_spline_alloc(gsl_interp_cspline, NPowerTable);
  logk = (double*)malloc(NPowerTable * sizeof(double));
  Pk   = (double*)malloc(NPowerTable * sizeof(double));

  for (i=0; i<NPowerTable; i++)
    {
      logk[i] = log10(exp(params.camb.Logk[i]));
      Pk  [i] = log10(exp(params.camb.LogPkref[i])) + 3.*logk[i];
    }

  gsl_spline_init(splinePk, logk, Pk, NPowerTable);

  free(Pk);
  free(logk);

  /* it re-initializes the linear growing mode */
  accelCAMB = gsl_interp_accel_alloc();
  splineGrowCAMB = gsl_spline_alloc (gsl_interp_cspline, params.camb.NCAMB);
  gsl_spline_init(splineGrowCAMB, params.camb.Scalef, params.camb.RefGM, params.camb.NCAMB);
  splineInvGrowCAMB = gsl_spline_alloc (gsl_interp_cspline, params.camb.NCAMB);
  gsl_spline_init(splineInvGrowCAMB, params.camb.RefGM, params.camb.Scalef, params.camb.NCAMB);

#ifdef OUTPUT_GM
  if (!ThisTask)
    {
      FILE *fd;
      fd=fopen("TwoGrowingModes","w");
      fprintf(fd,"# 1: scale factor\n# 2: D as computed from equations\n# 3: from CAMB\n");
      double a=params.camb.Scalef[0]/1.05;
      do
	{
	  a*=1.05;
	  if (a>1)
	    a=1.0;
	  fprintf(fd," %g %g %g\n",
		  a,
		  pow(10.,my_spline_eval(splineGrow, log10(a), accel)),
		  exp(my_spline_eval(splineGrowCAMB, a, accelCAMB))
		  );
	}
      while (a<1.0);
      fclose(fd);
    }
#endif

#endif

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


int read_Pk_table(void)
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

  accelPk = gsl_interp_accel_alloc();
  splinePk = gsl_spline_alloc(gsl_interp_cspline, NPowerTable);

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

  gsl_spline_init(splinePk, logk, Pk, NPowerTable);

  free(Pk);
  free(logk);

  return 0;
}

double PowerSpec_Tabulated(double k)
{
  return pow(10.,my_spline_eval(splinePk, log10(k), accelPk))/k/k/k;
}

double PowerSpec_Efstathiou(double k)
{
  return pow(k,params.PrimordialIndex) / pow(1 + pow(6.4 / ShapeGamma * k + pow(3.0 / ShapeGamma * k, 1.5) + pow(1.7 / ShapeGamma,2.0) * k * k, 1.13), 2 / 1.13);
}

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
  s   = 1.633*log((sqrt(1.+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)))/(keq*sqrt(Req)); // 2sqrt(6)/3=1.633
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
  double *rv, *lamv, *dldr, *lamvneg, *displv;

  rv        = (double*)malloc(NBIN * sizeof(double));
  lamv      = (double*)malloc(NBIN * sizeof(double));
  dldr      = (double*)malloc(NBIN * sizeof(double));
  lamvneg   = (double*)malloc(NBIN * sizeof(double));
  displv    = (double*)malloc(NBIN * sizeof(double));

  for (i=NBIN-1; i>=0; i--)
    {
      rv[i]=(rmin+i*dr);
      r=pow(10.,rv[i]);
      lamv[i]=log10(ComputeMassVariance(r));
      /* To allow interpolation of reverse function, 
	 the variance must increase with decreasing radius */
      if (i<NBIN-1 && lamv[i]-lamv[i+1]<1.e-6)
	lamv[i]=lamv[i+1]+1.e-6;   

      lamvneg[i]=-lamv[i];
      displv[i]=log10(ComputeDisplVariance(r));
    }

  for (i=0; i<NBIN; i++)
    {
      if (i==0)
        dldr[i]=(lamv[i+1]-lamv[i])/(rv[i+1]-rv[i]);
      else if (i==NBIN-1)
        dldr[i]=(lamv[i]-lamv[i-1])/(rv[i]-rv[i-1]);
      else
        dldr[i]=(lamv[i+1]-lamv[i-1])/(rv[i+1]-rv[i-1]);
    }

  /* initialization of splines for interpolation */
  gsl_spline_init(splineVar, rv, lamv, NBIN);
  gsl_spline_init(splineRadius, lamvneg, rv, NBIN);
  gsl_spline_init(splinedldr, rv, dldr, NBIN);
  gsl_spline_init(splineDVar, rv, displv, NBIN);

  free(displv  );
  free(lamvneg );
  free(dldr    );
  free(lamv    );
  free(rv      );

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
  k=exp(logk);
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

double IntegrandForDisplVariance(double logk, void *param)
{
  double w, k;
  k=exp(logk);
  w = WindowFunction(k * ThisRadius);
  return PowerSpectrum(k) * w*w * k / (2.*PI*PI);
}

double WindowFunction(double kr)
{
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
  return pow(10.,my_spline_eval(splineVar, log10(R), accel));
}


double dMassVariance_dr(double R)
{
  return my_spline_eval(splinedldr, log10(R), accel);
}

double DisplVariance(double R)
{
  return pow(10.,my_spline_eval(splineDVar, log10(R), accel));
}

double Radius(double Var)
{
  return pow(10.,my_spline_eval(splineRadius, -log10(Var), accel));
}


/**********************************/
/* COSMOLOGICAL FUNCTIONS SECTION */
/**********************************/

double Omega(double z)
{
  /* Cosmological mass density parameter as a function of redshift
     DIMENSIONLESS */

  return params.Omega0 * pow(1.+z,3.) * pow(Hubble(z)/100./params.Hubble100, -2);
}


double Hubble(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: km/s/Mpc  */
  double Esq,eqds;

  if (params.simpleLambda)
    Esq = params.Omega0*pow(1.+z,3.)+OmegaK*pow(1.+z,2.)+params.OmegaLambda;
  else
    {
      eqds=my_spline_eval(splineIntEoS, -log10(1.+z), accel);
      Esq = params.Omega0*pow(1.+z,3.)+OmegaK*pow(1.+z,2.)+params.OmegaLambda*pow(1.+z,3.) * exp(3.*eqds);
    }

  return 100.*params.Hubble100*sqrt(Esq);
}


double Hubble_Gyr(double z)
{
  /* Hubble parameter as a function of redshift
     DIMENSION: Gyr^-1 */

  return Hubble(z)/HUBBLETIME_GYR/100.;
}


double fomega(double z)
{
  /* Peebles' f(Omega) function, dlogD/dloga
     DIMENSIONLESS */

#ifndef SCALE_DEPENDENT_GROWTH
  return my_spline_eval(splinefO, -log10(1.+z), accel);
#else
  /* This is the implementation of scale-dependent growing mode */
  switch (SDGM.flag)
    {
    case -1:
      return my_spline_eval(splinefO, -log10(1.+z), accel);
      break;
    case 0:
    case 1:
      return VelfOmega(z);
      break;
    default:
      return 0;
      break;
    }
#endif
}

double fomega_2LPT(double z)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

#ifndef SCALE_DEPENDENT_GROWTH
  return my_spline_eval(splinefO2, -log10(1.+z), accel);
#else
  /* This is the implementation of scale-dependent growing mode */
  switch (SDGM.flag)
    {
    case -1:
      return my_spline_eval(splinefO2, -log10(1.+z), accel);
      break;
    case 0:
    case 1:
      return VelfOmega2(z);
      break;
    default:
      return 0;
      break;
    }
#endif
}

double fomega_3LPT_1(double z)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return my_spline_eval(splinefO31, -log10(1.+z), accel);
}

double fomega_3LPT_2(double z)
{
  /* second-order f(Omega) function, dlogD2/dloga
     DIMENSIONLESS */

  return my_spline_eval(splinefO32, -log10(1.+z), accel);
}

double GrowingMode(double z)
{
  /* linear growing mode, interpolation on the grid
     DIMENSIONLESS */

#ifndef SCALE_DEPENDENT_GROWTH
  return pow(10.,my_spline_eval(splineGrow, -log10(1.+z), accel));
#else
  /* This is the implementation of scale-dependent growing mode */
  switch (SDGM.flag)
    {
    case -1:
      return exp(my_spline_eval(splineGrowCAMB, 1./(1.+z), accelCAMB));
      break;
    case 0:
      return MatterGrowingMode(z);
      break;
    case 1:
      return VelGrowingMode(z);
      break;
    default:
      return 0;
      break;
    }
#endif
}


double GrowingMode_2LPT(double z)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

#ifndef SCALE_DEPENDENT_GROWTH
  return pow(10.,my_spline_eval(splineGrow2, -log10(1.+z), accel));
#else
  /* This is a rough implementation of scale-dependent 2nd-order growing mode */
  switch (SDGM.flag)
    {
    case -1:
      return pow(10.,my_spline_eval(splineGrow2, -log10(1.+z), accel));
      break;
    case 0:
    case 1:
      return 3./7.*pow(GrowingMode(z),2.0) * pow(Omega(z),-0.007);
      break;
    default:
      return 0;
      break;
    }
#endif
}

double GrowingMode_3LPT_1(double z)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return -pow(10.,my_spline_eval(splineGrow31, -log10(1.+z), accel));
}

double GrowingMode_3LPT_2(double z)
{
  /* second-order growing mode, interpolation on the grid
     DIMENSIONLESS */

  return pow(10.,my_spline_eval(splineGrow32, -log10(1.+z), accel));
}

double InverseGrowingMode(double D)
{
  /* redshift corresponding to a linear growing mode, interpolation on the grid
     DIMENSIONLESS */

#ifndef SCALE_DEPENDENT_GROWTH
  return 1./pow(10.,my_spline_eval(splineInvGrow, log10(D), accel)) -1.;
#else
  /* This is the implementation of scale-dependent growing mode */
  switch (SDGM.flag)
    {
    case -1:
      return 1./my_spline_eval(splineInvGrowCAMB, log(D), accel) -1.;
      break;
    case 0:
      return InverseMatterGrowingMode(D);
      break;
    case 1:
      return 0; //VelocityInverseGrowingMode(z);
      break;
    default:
      return 0;
      break;
    }
#endif
}


double CosmicTime(double z)
{
  /* cosmic time, interpolation on the grid
     Gyr */

  return pow(10.,my_spline_eval(splineTime, -log10(1.+z), accel));
}

double InverseCosmicTime(double t)
{
  /* scale factor corresponding to a cosmic time, interpolation on the grid
     Gyr */

  return pow(10.,my_spline_eval(splineInvTime, log10(t), accel));
}

double ProperDistance(double z)
{
  /* proper distance, interpolation on the grid
     Mpc */

  return my_spline_eval(splineDist, -log10(1.+z), accel);
}

double InverseProperDistance(double d)
{
  /* redshift corresponding to a proper distance, interpolation on the grid
     DIMENSIONLESS */

  return 1./pow(10.,my_spline_eval(splineInvDist, -d, accel))-1;
}

double dProperDistance_dz(double z)
{
  /* cosmic time, interpolation on the grid
     Gyr */

  return my_spline_eval(splinedrdz, -log10(1.+z), accel);
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

double TypicalCollapsingMass(double z)
{
  /* NB: this is never used... 
     for SCALE_INDEPENDENT_GROWTH one should specify the scale */
  return MassForSize(Radius(pow(DELTA_C/GrowingMode(z),2.0)));
}

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
  SDGM.flag=0;       // mass growing mode
  SDGM.radius=r;
#endif
  D=GrowingMode(z);

  return MatterDensity * dOmega_dVariance(MassVariance(r)*D*D,z)
    * fabs(dMassVariance_dr(r) /6.0) /mass /mass;

}

/*
  Table of cosmo.F90 and cosmo.c cosmological functions

  C version                     F90 version

  Omega(z)                      Omega(z)
  Hubble(z)                     Hubble(z)
  Hubble_Gyr(z)                 Hubble_Gyr(z)
  fomega(z)                     fomega(z)
  fomega_2LPT(z)                fomega_2LPT(z)
  CosmicTime(z)                 time(z) or cosmic_time(z)
  InverseCosmicTime(t)          scale(t)
  GrowingMode(z)                grow(z)
  GrowingMode_2LPT(z)           grow_2LPT(z)
  InverseGrowingMode(D)         inverse_grow(D)
  ProperDistance(z)             proper_dist(z)
  InverseProperDistance(r)      inverse_proper_dist(r)
  dProperDistance_dz(z)         dproper_dist_dz(z)
  PowerSpectrum(k)              power(k)
  MassVariance(r)               aLambda(r)
  dMassVariance_dr(r)           dLambda(r)
  DisplVariance(r)              --
  Radius(v)                     Radius(v)
  SizeForMass(m)                Size(m)
  MassForSize(r)                --
  TypicalCollapsingMass(z)      aM_star(z)
  dOmega_dVariance(v,z)         dodl(v,z)
  AnalyticMassFunction(m,z)     an_h(v,z)
  WindowFunction(kr)            window(kr)

  MassVariance(SizeForMass(m))  Variance(m)
  MassForSize(Radius(v))        aMass(v)

 */

double my_spline_eval(gsl_spline *spline, double x, gsl_interp_accel *accel)
{
  if (x<spline->x[0])
    return spline->y[0]+(x-spline->x[0])*(spline->y[1]-spline->y[0])/(spline->x[1]-spline->x[0]);
  else if (x>spline->x[spline->size-1])
    return spline->y[spline->size-1]+(x-spline->x[spline->size-1])*
      (spline->y[spline->size-1]-spline->y[spline->size-2])/(spline->x[spline->size-1]-spline->x[spline->size-2]);
  else
    return gsl_spline_eval(spline,x,accel);
}

