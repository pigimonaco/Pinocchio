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

#ifdef SCALE_DEPENDENT_GROWTH

#include "pinocchio.h"

#define ONLY_CDM_TRANSF
//#define OUTPUT_GM
//#define NOBARYONS

static double *StoredLogK, *CAMBRedshifts, *CAMBScalefac, *RefGrowingMode,
  *StoredLogTotalPowerSpectrum;

static gsl_interp_accel *accelGrow=0x0, *accel=0x0;
static gsl_spline **splineGrowMatter, **splineGrowVel, **splineInvGrowMatter, **splineFomega, **splineFomega2, *spline=0x0;

int ThisOutput, ThisRadius;

double IntegrandForSDMassVariance(double, void *);
double IntegrandForSDVelVariance(double, void *);
double PowerFromCAMB(double);

int read_power_table_from_CAMB()
{
  /* this routine reads P(k) at various redshifts from a series of CAMB outputs 
     and stores them in memory */

  int i,j,dummy,ff;
  double kappa,Pk;
  char filename[BLENGTH],buffer[BLENGTH],*ugo;
  FILE *fd;
#ifdef ONLY_CDM_TRANSF
  gsl_spline *splineTransf=0x0;
  gsl_interp_accel *accel=0x0;
  int Nktransf;
  double *Ktransf,*Ttransf,Tcdm,Ttot,Tbar;
  FILE *fd2;
#endif

  if (!ThisTask)
    {
      /* counts the number of CAMB files */
      params.camb.NCAMB=0;
      sprintf(filename,"%s_%s_%03d.dat",params.camb.RunName,params.camb.MatterFile,params.camb.NCAMB);
      while ((fd=fopen(filename,"r"))!=0x0)
	{
	  if (!params.camb.NCAMB)
	    /* if this is the first opened file, count the number of data lines */
	    {
	      params.camb.Nkbins=0;
	      while(!feof(fd))
		{
		  ugo=fgets(buffer,BLENGTH,fd);
		  i=sscanf(buffer,"%lf",&kappa);
		  if (i && ugo!=0x0)
		    params.camb.Nkbins++;
		}
	    }

#ifdef ONLY_CDM_TRANSF
	  /* open the transfer file to construct the P(k) of matter */
	  sprintf(filename,"%s_%s_000.dat",params.camb.RunName,params.camb.TransferFile);
	  if ((fd2=fopen(filename,"r"))!=0x0)
	    {
	      if (!params.camb.NCAMB)
		/* if this is the first opened file, count again the number of data lines */
		{
		  Nktransf=0;
		  while(!feof(fd2))
		    {
		      ugo=fgets(buffer,BLENGTH,fd2);
		      i=sscanf(buffer,"%lf",&kappa);
		      if (i && ugo!=0x0)
			Nktransf++;
		    }
		}
	      fclose(fd2);
	    }
	  else
	    {
	      printf("Error on Task 0: CAMB transfer file %s not found\n",filename);
	      return 1;
	    }
#endif
	  fclose(fd);
	  params.camb.NCAMB++;
	  sprintf(filename,"%s_%s_%03d.dat",params.camb.RunName,params.camb.MatterFile,params.camb.NCAMB);
	}

      if (!params.camb.NCAMB)
	{
	  printf("Error on Task 0: CAMB file %s not found\n",filename);
	  return 1;
	}
      else if (!params.camb.Nkbins)
	{
	  sprintf(filename,"%s_%s_%03d.dat",params.camb.RunName,params.camb.MatterFile,0);
	  printf("Error on Task 0: problem in reading CAMB file %s\n",filename);
	  return 1;
	}
#ifdef ONLY_CDM_TRANSF
      else if (!Nktransf)
	{
	  printf("Error on Task 0: no lines found in transfer function file\n");
	  return 1;
	}
#endif
      else if (params.camb.ReferenceOutput>params.camb.NCAMB-1)
	{
	  printf("Error on Task 0: ReferenceOutput is larger than NCAMB-1\n");
	  return 1;
	}

      printf("Found %d CAMB matter power files with %d lines each\n",
	     params.camb.NCAMB,params.camb.Nkbins);
#ifdef ONLY_CDM_TRANSF
      printf("Number of data lines in transfer function file: %d\n",Nktransf);
#endif
    }

  /* Task 0 broadcasts NCAMB and Nkbins */
  MPI_Bcast(&params.camb.NCAMB,  sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&params.camb.Nkbins, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
#ifdef ONLY_CDM_TRANSF
  MPI_Bcast(&Nktransf, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* allocates the needed memory for the power spectra */
  StoredLogTotalPowerSpectrum=(double*)malloc(params.camb.NCAMB *params.camb.Nkbins * sizeof(double));
  StoredLogK=(double*)malloc(params.camb.Nkbins * sizeof(double));
  CAMBRedshifts=(double*)malloc(params.camb.NCAMB * sizeof(double));
  CAMBScalefac=(double*)malloc(params.camb.NCAMB * sizeof(double));
  RefGrowingMode=(double*)malloc(params.camb.NCAMB * sizeof(double));

  /* Task 0 reads all the CAMB files and stores the power spectra */
  if (!ThisTask)
    {

#ifdef ONLY_CDM_TRANSF
      /* this is used only by Task 0 */
      accel = gsl_interp_accel_alloc();
      splineTransf = gsl_spline_alloc(gsl_interp_cspline, Nktransf);
      Ktransf=(double*)malloc(Nktransf*sizeof(double));
      Ttransf=(double*)malloc(Nktransf*sizeof(double));
#endif

      for (i=0; i<params.camb.NCAMB; i++)
	{
#ifdef ONLY_CDM_TRANSF
	  sprintf(filename,"%s_%s_%03d.dat",params.camb.RunName,params.camb.TransferFile,i);
	  fd=fopen(filename,"r");
	  for (j=0; j<Nktransf; j++)
	    {
	      ugo=fgets(buffer,BLENGTH,fd);
	      sscanf(buffer,"%lf %lf %lf %*f %*f %*f %lf",&kappa,&Tcdm,&Tbar,&Ttot);
	      Ktransf[j]=log(kappa);
#ifdef NOBARYONS
	      Ttransf[j]=2.*log(Tcdm/Ttot);
#else
	      Ttransf[j]=2.*log(((params.Omega0-params.OmegaBaryon)*Tcdm+params.OmegaBaryon*Tbar)/(params.Omega0)/Ttot);
#endif
	    }
	  gsl_spline_init(splineTransf, Ktransf, Ttransf, Nktransf);
	  fclose(fd);
#endif
	  sprintf(filename,"%s_%s_%03d.dat",params.camb.RunName,params.camb.MatterFile,i);
	  fd=fopen(filename,"r");
	  for (j=0; j<params.camb.Nkbins; j++)
	    {
	      ff=fscanf(fd,"%lf %lf",&kappa,&Pk);
	      StoredLogTotalPowerSpectrum[i*params.camb.Nkbins +j]=
		log(Pk/pow(params.Hubble100,3.))
#ifdef ONLY_CDM_TRANSF
		/* here it corrects the power spectrum by the square of the ratio 
		   of the matter and total transfer functions */
		+my_spline_eval(splineTransf, log(kappa), accel)
#endif
		;
	      if (!i)
		StoredLogK[j]=log(kappa*params.Hubble100);
	    }
	  fclose(fd);
	}
#ifdef ONLY_CDM_TRANSF
      free(Ktransf);
      free(Ttransf);
      gsl_spline_free(splineTransf);
      gsl_interp_accel_free(accel);
#endif

      if ((fd=fopen(params.camb.RedshiftsFile,"r"))==0x0)
	{
	  printf("Error: Redshift file %s not found\n",params.camb.RedshiftsFile);
	  return 1;
	}
      for (i=0; i<params.camb.NCAMB; i++)
	ff=fscanf(fd,"%d %lf",&dummy,CAMBRedshifts+i);
      fclose(fd);

      /* The last redshift MUST be z=0 */
      if (CAMBRedshifts[params.camb.NCAMB-1]!=0.0)
	{
	  printf("ERROR on Task 0: last CAMB redsbift must be 0.0\n");
	  return 1;
	}
    }

  /* broadcast of loaded and computed quantities */
  MPI_Bcast(StoredLogK, params.camb.Nkbins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(StoredLogTotalPowerSpectrum, params.camb.NCAMB*params.camb.Nkbins, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(CAMBRedshifts, params.camb.NCAMB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  for (i=0; i<params.camb.NCAMB; i++)
    {
      CAMBScalefac[i]=1./(CAMBRedshifts[i]+1.);
      RefGrowingMode[i]=0.5*(StoredLogTotalPowerSpectrum[i * params.camb.Nkbins + params.camb.ReferenceScale] - 
			     StoredLogTotalPowerSpectrum[params.camb.ReferenceOutput * params.camb.Nkbins + params.camb.ReferenceScale]);
    }
  params.camb.ReferenceRedshift=CAMBRedshifts[params.camb.ReferenceOutput];

  /* make P(k) at the reference redshift available to cosmo.c */
  params.camb.Logk = StoredLogK;
  params.camb.LogPkref = StoredLogTotalPowerSpectrum +
    params.camb.ReferenceOutput * params.camb.Nkbins;
  params.camb.D2ref = exp(StoredLogTotalPowerSpectrum[(params.camb.NCAMB-1) * params.camb.Nkbins + params.camb.ReferenceScale] - 
			  StoredLogTotalPowerSpectrum[params.camb.ReferenceOutput * params.camb.Nkbins + params.camb.ReferenceScale]);
  params.camb.Scalef = CAMBScalefac;
  params.camb.RefGM = RefGrowingMode;

  return 0;
}


int initialize_ScaleDependentGrowth(void)
{

  /* Here it computes the scale-dependent growing modes at first and second order,
     for density and velocity, and the relative fomegas */

  int i, i1, i2;
  double *VarMatter, *VarVel, *ThisMatter, *ThisVel, *ThisfO, *ThisfO2;
  double result,error;
  gsl_function FuncMatter, FuncVel;
  double tolerance=1.e-4;

  /* functions to integrate */
  FuncMatter.function = &IntegrandForSDMassVariance;
  FuncVel.function = &IntegrandForSDVelVariance;

  /* initialization of splines for growing modes */
  accelGrow           = gsl_interp_accel_alloc();
  splineGrowMatter    = (gsl_spline**)calloc(Smoothing.Nsmooth,sizeof(gsl_spline*));
  splineGrowVel       = (gsl_spline**)calloc(Smoothing.Nsmooth,sizeof(gsl_spline*));
  splineInvGrowMatter = (gsl_spline**)calloc(Smoothing.Nsmooth,sizeof(gsl_spline*));
  splineFomega        = (gsl_spline**)calloc(Smoothing.Nsmooth,sizeof(gsl_spline*));
  splineFomega2       = (gsl_spline**)calloc(Smoothing.Nsmooth,sizeof(gsl_spline*));
  for (i=0; i<Smoothing.Nsmooth; i++)
    {
      splineGrowMatter   [i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
      splineGrowVel      [i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
      splineInvGrowMatter[i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
      splineFomega       [i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
      splineFomega2      [i] = gsl_spline_alloc(gsl_interp_cspline, params.camb.NCAMB);
    }

  /* computes the needed growing modes and stores them */
  VarMatter =(double*)malloc(Smoothing.Nsmooth*sizeof(double));
  VarVel    =(double*)malloc(Smoothing.Nsmooth*sizeof(double));
  ThisMatter=(double*)malloc(params.camb.NCAMB*sizeof(double));
  ThisVel   =(double*)malloc(params.camb.NCAMB*sizeof(double));
  ThisfO    =(double*)malloc(params.camb.NCAMB*sizeof(double));
  ThisfO2   =(double*)malloc(params.camb.NCAMB*sizeof(double));

  /* This is needed by PowerFromCAMB */
  accel = gsl_interp_accel_alloc ();
  spline = gsl_spline_alloc (gsl_interp_cspline, params.camb.Nkbins);

  /* values at the reference redshift */
  ThisOutput=params.camb.ReferenceOutput;
  for (ThisRadius=0; ThisRadius<Smoothing.Nsmooth; ThisRadius++)
    {
      gsl_integration_qags (&FuncMatter, StoredLogK[0], StoredLogK[params.camb.Nkbins-1], 0.0, tolerance, NWINT, workspace, &result, &error);
      VarMatter[ThisRadius]=result;
      gsl_integration_qags (&FuncVel, StoredLogK[0], StoredLogK[params.camb.Nkbins-1], 0.0, tolerance, NWINT, workspace, &result, &error);
      VarVel[ThisRadius]=result;
    }

#ifdef OUTPUT_GM
  SDGM.flag=-1;
  FILE *fd;
  if (!ThisTask) 
    fd=fopen("GrowingModes","w");
#endif

  /* computation of growing modes for matter and velocity as a function of z,
     for all smoothing radii */
  for (ThisRadius=0; ThisRadius<Smoothing.Nsmooth; ThisRadius++)
    {
      for (ThisOutput=0; ThisOutput<params.camb.NCAMB; ThisOutput++)
	{
	  gsl_integration_qags (&FuncMatter, StoredLogK[0], StoredLogK[params.camb.Nkbins-1], 0.0, tolerance, NWINT, workspace, &result, &error);
	  ThisMatter[ThisOutput]=sqrt(result/VarMatter[ThisRadius]);
	  gsl_integration_qags (&FuncVel, StoredLogK[0], StoredLogK[params.camb.Nkbins-1], 0.0, tolerance, NWINT, workspace, &result, &error);
	  ThisVel[ThisOutput]=sqrt(result/VarVel[ThisRadius]);
	}
      /* initialization of splines */
      gsl_spline_init(splineGrowMatter[ThisRadius],    CAMBScalefac, ThisMatter,   params.camb.NCAMB);
      gsl_spline_init(splineGrowVel[ThisRadius],       CAMBScalefac, ThisVel,      params.camb.NCAMB);
      gsl_spline_init(splineInvGrowMatter[ThisRadius], ThisMatter,   CAMBScalefac, params.camb.NCAMB);

      /* fomega */
      for (ThisOutput=0; ThisOutput<params.camb.NCAMB; ThisOutput++)
	{


	  if (ThisOutput<params.camb.NCAMB-1)
	    {
	      i1=(ThisOutput>0 ? ThisOutput-1 : 0);
	      i2=ThisOutput+1;  //(ThisOutput<params.camb.NCAMB-1 ? ThisOutput+1 : params.camb.NCAMB-1);

	      ThisfO[ThisOutput] = (ThisVel[i2]-ThisVel[i1]) / (CAMBScalefac[i2]-CAMBScalefac[i1]) * CAMBScalefac[ThisOutput] / ThisVel[ThisOutput];
	      ThisfO2[ThisOutput] = ( (3./7.*pow(ThisVel[i2],2.0)*pow(Omega(CAMBRedshifts[i2]),-0.007)) - 
				      (3./7.*pow(ThisVel[i1],2.0)*pow(Omega(CAMBRedshifts[i1]),-0.007)) )
		/ (CAMBScalefac[i2]-CAMBScalefac[i1]) * CAMBScalefac[ThisOutput] / 
		(3./7.*pow(ThisVel[ThisOutput],2.0)*pow(Omega(CAMBRedshifts[ThisOutput]),-0.007));
	    }
	  else
	    {
	      ThisfO[ThisOutput] = - ThisfO[ThisOutput-2] * (CAMBScalefac[ThisOutput]-CAMBScalefac[ThisOutput-1])/(CAMBScalefac[ThisOutput-1]-CAMBScalefac[ThisOutput-2]) 
		+ ThisfO[ThisOutput-1] *  (CAMBScalefac[ThisOutput]-CAMBScalefac[ThisOutput-2])/(CAMBScalefac[ThisOutput-1]-CAMBScalefac[ThisOutput-2]) ;
	      ThisfO2[ThisOutput] = - ThisfO2[ThisOutput-2] * (CAMBScalefac[ThisOutput]-CAMBScalefac[ThisOutput-1])/(CAMBScalefac[ThisOutput-1]-CAMBScalefac[ThisOutput-2]) 
		+ ThisfO2[ThisOutput-1] *  (CAMBScalefac[ThisOutput]-CAMBScalefac[ThisOutput-2])/(CAMBScalefac[ThisOutput-1]-CAMBScalefac[ThisOutput-2]) ;
	    }

#ifdef OUTPUT_GM
	  if (!ThisTask) 
	    fprintf(fd," %d %d   %g %g   %g %g %g   %g %g %g   %g %g %f  %g %g\n",
		    ThisRadius,ThisOutput,
		    CAMBRedshifts[ThisOutput],Smoothing.Radius[ThisRadius],
		    ThisMatter[ThisOutput],ThisVel[ThisOutput],GrowingMode(CAMBRedshifts[ThisOutput]),
		    3./7.*pow(ThisMatter[ThisOutput],2.0)*pow(Omega(CAMBRedshifts[ThisOutput]),-0.007),
		    3./7.*pow(ThisVel[ThisOutput],2.0)*pow(Omega(CAMBRedshifts[ThisOutput]),-0.007),
		    GrowingMode_2LPT(CAMBRedshifts[ThisOutput]),
		    ThisfO[ThisOutput], fomega(CAMBRedshifts[ThisOutput]), Omega(CAMBRedshifts[ThisOutput]),
		    ThisfO2[ThisOutput], fomega_2LPT(CAMBRedshifts[ThisOutput])
		    );
#endif

	}
      /* initialization of splines */
      gsl_spline_init(splineFomega[ThisRadius],  CAMBScalefac, ThisfO,  params.camb.NCAMB);
      gsl_spline_init(splineFomega2[ThisRadius], CAMBScalefac, ThisfO2, params.camb.NCAMB);

    }

#ifdef OUTPUT_GM
  if (!ThisTask) fclose(fd);
#endif

  gsl_spline_free(spline);
  gsl_interp_accel_free(accel);

  free(ThisfO2);
  free(ThisfO);
  free(ThisVel);
  free(ThisMatter);
  free(VarVel);
  free(VarMatter);

  return 0;
}


double IntegrandForSDMassVariance(double logk, void *param)
{
  double k,R;

  k=exp(logk);
  R=(ThisRadius<Smoothing.Nsmooth-1 ? Smoothing.Radius[ThisRadius] : params.InterPartDist/6.);
  return PowerFromCAMB(k) * exp(-k*k*R*R) * k*k*k / (2.*PI*PI);
}


double IntegrandForSDVelVariance(double logk, void *radius)
{
  double k,R;

  k=exp(logk);
  R=(ThisRadius<Smoothing.Nsmooth-1 ? Smoothing.Radius[ThisRadius] : params.InterPartDist/6.);
  return PowerFromCAMB(k) * exp(-k*k*R*R) * k / (2.*PI*PI);
}


double PowerFromCAMB(double k)
{
  /* This function gives the interpolated power spectrum from the CAMB outputs.
     The power spectrum is interpolated from */

  static int LastOutput=-99;

  /* initialize gsl spline the first time the function is called */
  if (spline==0x0 || ThisOutput!=LastOutput)
    {
      gsl_spline_init (spline, StoredLogK,
		       StoredLogTotalPowerSpectrum + ThisOutput * params.camb.Nkbins, 
		       params.camb.Nkbins);
    }

  LastOutput=ThisOutput;
  double result=exp(my_spline_eval (spline, log(k), accel));
  return result;
}


double MatterGrowingMode(double z)
{

  int mysmooth, interpolate;
  double w;

  /* computes the smoothing radius to be used (including the case of interpolation) */
  if (SDGM.ismooth>=0 && SDGM.ismooth<Smoothing.Nsmooth)
    {
      mysmooth=SDGM.ismooth;
      interpolate=0;
    }
  else
    {
      if (SDGM.radius>Smoothing.Radius[0])
	{
	  mysmooth=0;
	  interpolate=0;
	}
      else if (SDGM.radius<Smoothing.Radius[Smoothing.Nsmooth-1])
	{
	  mysmooth=Smoothing.Nsmooth-1;
	  interpolate=0;
	}
      else
	{
	  for (mysmooth=1; mysmooth<Smoothing.Nsmooth && SDGM.radius<=Smoothing.Radius[mysmooth]; mysmooth++)
	    ;
	  interpolate=1;
	}
    }

  if (interpolate)
    {
      w=(SDGM.radius-Smoothing.Radius[mysmooth])/(Smoothing.Radius[mysmooth-1]-Smoothing.Radius[mysmooth]);
      return ((1.-w) * my_spline_eval(splineGrowMatter[mysmooth  ], 1./(1.+z), accelGrow)
	      +   w  * my_spline_eval(splineGrowMatter[mysmooth-1], 1./(1.+z), accelGrow));
    }
  else
    return my_spline_eval(splineGrowMatter[mysmooth], 1./(1.+z), accelGrow);
}


double VelGrowingMode(double z)
{

  int mysmooth, interpolate;
  double w;

  /* computes the smoothing radius to be used (including the case of interpolation) */
  if (SDGM.ismooth>=0 && SDGM.ismooth<Smoothing.Nsmooth)
    {
      mysmooth=SDGM.ismooth;
      interpolate=0;
    }
  else
    {
      if (SDGM.radius>Smoothing.Radius[0])
	{
	  mysmooth=0;
	  interpolate=0;
	}
      else if (SDGM.radius<Smoothing.Radius[Smoothing.Nsmooth-1])
	{
	  mysmooth=Smoothing.Nsmooth-1;
	  interpolate=0;
	}
      else
	{
	  for (mysmooth=1; mysmooth<Smoothing.Nsmooth && SDGM.radius<=Smoothing.Radius[mysmooth]; mysmooth++)
	    ;
	  interpolate=1;
	}
    }

  if (interpolate)
    {
      w=(SDGM.radius-Smoothing.Radius[mysmooth])/(Smoothing.Radius[mysmooth-1]-Smoothing.Radius[mysmooth]);
      return ((1.-w) * my_spline_eval(splineGrowVel[mysmooth  ], 1./(1.+z), accelGrow)
	      +   w  * my_spline_eval(splineGrowVel[mysmooth-1], 1./(1.+z), accelGrow));
    }
  else
    return my_spline_eval(splineGrowVel[mysmooth], 1./(1.+z), accelGrow);
}


double InverseMatterGrowingMode(double D)
{

  int mysmooth, interpolate;
  double w;

  /* computes the smoothing radius to be used (including the case of interpolation) */
  if (SDGM.ismooth>=0 && SDGM.ismooth<Smoothing.Nsmooth)
    {
      mysmooth=SDGM.ismooth;
      interpolate=0;
    }
  else
    {
      if (SDGM.radius>Smoothing.Radius[0])
	{
	  mysmooth=0;
	  interpolate=0;
	}
      else if (SDGM.radius<Smoothing.Radius[Smoothing.Nsmooth-1])
	{
	  mysmooth=Smoothing.Nsmooth-1;
	  interpolate=0;
	}
      else
	{
	  for (mysmooth=1; mysmooth<Smoothing.Nsmooth && SDGM.radius<=Smoothing.Radius[mysmooth]; mysmooth++)
	    ;
	  interpolate=1;
	}
    }

  if (interpolate)
    {
      w=(SDGM.radius-Smoothing.Radius[mysmooth])/(Smoothing.Radius[mysmooth-1]-Smoothing.Radius[mysmooth]);
      return 1./((1.-w) * my_spline_eval(splineInvGrowMatter[mysmooth  ], D, accelGrow)
		 +   w  * my_spline_eval(splineInvGrowMatter[mysmooth-1], D, accelGrow)) -1.;
    }
  else
    return 1./my_spline_eval(splineInvGrowMatter[mysmooth], D, accelGrow) -1 ;
}


double VelfOmega(double z)
{

  int mysmooth, interpolate;
  double w;

  /* computes the smoothing radius to be used (including the case of interpolation) */
  if (SDGM.ismooth>=0 && SDGM.ismooth<Smoothing.Nsmooth)
    {
      mysmooth=SDGM.ismooth;
      interpolate=0;
    }
  else
    {
      if (SDGM.radius>Smoothing.Radius[0])
	{
	  mysmooth=0;
	  interpolate=0;
	}
      else if (SDGM.radius<Smoothing.Radius[Smoothing.Nsmooth-1])
	{
	  mysmooth=Smoothing.Nsmooth-1;
	  interpolate=0;
	}
      else
	{
	  for (mysmooth=1; mysmooth<Smoothing.Nsmooth && SDGM.radius<=Smoothing.Radius[mysmooth]; mysmooth++)
	    ;
	  interpolate=1;
	}
    }

  if (interpolate)
    {
      w=(SDGM.radius-Smoothing.Radius[mysmooth])/(Smoothing.Radius[mysmooth-1]-Smoothing.Radius[mysmooth]);
      return ((1.-w) * my_spline_eval(splineFomega[mysmooth  ], 1./(1.+z), accelGrow)
	      +   w  * my_spline_eval(splineFomega[mysmooth-1], 1./(1.+z), accelGrow));
    }
  else
    return my_spline_eval(splineFomega[mysmooth], 1./(1.+z), accelGrow);
}


double VelfOmega2(double z)
{

  int mysmooth, interpolate;
  double w;

  /* computes the smoothing radius to be used (including the case of interpolation) */
  if (SDGM.ismooth>=0 && SDGM.ismooth<Smoothing.Nsmooth)
    {
      mysmooth=SDGM.ismooth;
      interpolate=0;
    }
  else
    {
      if (SDGM.radius>Smoothing.Radius[0])
	{
	  mysmooth=0;
	  interpolate=0;
	}
      else if (SDGM.radius<Smoothing.Radius[Smoothing.Nsmooth-1])
	{
	  mysmooth=Smoothing.Nsmooth-1;
	  interpolate=0;
	}
      else
	{
	  for (mysmooth=1; mysmooth<Smoothing.Nsmooth && SDGM.radius<=Smoothing.Radius[mysmooth]; mysmooth++)
	    ;
	  interpolate=1;
	}
    }

  if (interpolate)
    {
      w=(SDGM.radius-Smoothing.Radius[mysmooth])/(Smoothing.Radius[mysmooth-1]-Smoothing.Radius[mysmooth]);
      return ((1.-w) * my_spline_eval(splineFomega2[mysmooth  ], 1./(1.+z), accelGrow)
	      +   w  * my_spline_eval(splineFomega2[mysmooth-1], 1./(1.+z), accelGrow));
    }
  else
    return my_spline_eval(splineFomega2[mysmooth], 1./(1.+z), accelGrow);
}
#endif

