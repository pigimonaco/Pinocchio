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
#include "fragment.h"

#define FTOZ(A) (A>0.0 ? A-1 : A);

#define FORTRAN

int compute_mf(int iout)
{
  /* computes the mass function of the groups and the analytic mass function */

  int i,ibin,idummy;
  char filename[LBLENGTH],lab1[10],lab2[10];
  double amass,x,dm,a,a1,a2,a3,m,D,mx,r,massvar,sigma,ni,fdummy;
  FILE *file;

  D=GrowingMode(outputs.z[iout],params.k_for_GM);

  /* sets counters to zero */
  for (i=0; i<mf.NBIN; i++)
    {
      mf.ninbin_local[i]=0;
      mf.massinbin_local[i]=0.0;
      mf.ninbin[i]=0;
      mf.massinbin[i]=0.0;
    }

  /* constructs the histograms; this is done for all slices */
  for (i=FILAMENT+1; i<=ngroups; i++)
    if (groups[i].point>0 && groups[i].good && groups[i].Mass >= params.MinHaloMass)
      {
	amass=groups[i].Mass*params.ParticleMass;
        ibin=(int)((log10(amass)-mf.mmin)/DELTAM);
        if (ibin < 0 || ibin >= mf.NBIN)
	  continue;
        mf.ninbin_local[ibin]++;
        mf.massinbin_local[ibin]+=amass;
      }

  /* sums counters over all tasks */
  MPI_Reduce(mf.ninbin_local, mf.ninbin, mf.NBIN, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(mf.massinbin_local, mf.massinbin, mf.NBIN, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* Writes result */
  if (!ThisTask)
    {
      /* if this is not the first slice, read the values that have already been stored */
      if (ThisSlice)
	{
	  sprintf(filename,"pinocchio.%6.4f.%s.mf.save",outputs.z[iout],params.RunFlag);
	  printf("[%s] Reading data from file  %s\n",fdate(),filename);
	  file=fopen(filename,"r");
	  for (i=0; i<mf.NBIN; i++)
	    {
	      (void)fread(&idummy,sizeof(int),1,file);
	      (void)fread(&fdummy,sizeof(double),1,file);
	      mf.ninbin[i]+=idummy;
	      mf.massinbin[i]+=fdummy;
	    }
	  fclose(file);
	  /* if this is the last slice remove the .save file */
	  if (ThisSlice==NSlices-1)
	    if (remove(filename)!=0)
	      printf("WARNING: Task %d could not remove file %s\n",ThisTask,filename);
	}

      /* if this is not the last slice, write the stored values */
      if (ThisSlice<NSlices-1)
	{
	  sprintf(filename,"pinocchio.%6.4f.%s.mf.save",outputs.z[iout],params.RunFlag);
	  printf("[%s] Storing data for slice %d into file %s\n",fdate(),ThisSlice,filename);
	  file=fopen(filename,"w");
	  for (i=0; i<mf.NBIN; i++)
	    {
	      fwrite(&mf.ninbin[i],sizeof(int),1,file);
	      fwrite(&mf.massinbin[i],sizeof(double),1,file);
	    }
	  fclose(file);
	}

      /* At the final slice writes the complete mass function */
      if (ThisSlice==NSlices-1)
	{
	  sprintf(filename,"pinocchio.%6.4f.%s.mf.out",outputs.z[iout],params.RunFlag);
	  printf("[%s] Writing mass function into file  %s\n",fdate(),filename);
	  file=fopen(filename,"w");

	  fprintf(file,"# Mass function for redshift %f\n",outputs.z[iout]);
	  if (params.OutputInH100)
	    {
	      strcpy(lab1,"/h");
	      strcpy(lab2,"h^4");
	    }
	  else
	    {
	      strcpy(lab1,"");
	      strcpy(lab2,"");
	    }

	  fprintf(file,"# 1) mass (Msun%s)\n",lab1);
	  fprintf(file,"# 2) n(m) (Mpc^-3 Msun^-1 %s)\n",lab2);
	  fprintf(file,"# 3) upper +1-sigma limit for n(m) (Mpc^-3 Msun^-1 %s)\n",lab2);
	  fprintf(file,"# 4) lower -1-sigma limit for n(m) (Mpc^-3 Msun^-1 %s)\n",lab2);
	  fprintf(file,"# 5) number of halos in the bin\n");
	  switch(params.AnalyticMassFunction)
	    {
	    case 0:
	      fprintf(file,"# 6) analytical n(m) from Press & Schechter 1974\n");
	      break;
	    case 1:
	      fprintf(file,"# 6) analytical n(m) from Sheth & Tormen 2001\n");
	      break;
	    case 2:
	      fprintf(file,"# 6) analytical n(m) from Jenkins et al. 2001\n");
	      break;
	    case 3:
	      fprintf(file,"# 6) analytical n(m) from Warren et al. 2006\n");
	      break;
	    case 4:
	      fprintf(file,"# 6) analytical n(m) from Reed et al. 2007\n");
	      break;
	    case 5:
	      fprintf(file,"# 6) analytical n(m) from Crocce et al. 2010\n");
	      break;
	    case 6:
	      fprintf(file,"# 6) analytical n(m) from Tinker et al. 2010\n");
	      break;
	    case 7:
	      fprintf(file,"# 6) analytical n(m) from Courtin et al. 2010\n");
	      break;
	    case 8:
	      fprintf(file,"# 6) analytical n(m) from Angulo et al. 2012\n");
	      break;
	    case 9:
	      fprintf(file,"# 6) analytical n(m) from Watson et al. 2013\n");
	      break;
	    case 10:
	      fprintf(file,"# 6) analytical n(m) from Crocce et al. 2010, universal\n");
	      break;
	    default:
	      break;
	    }
	  fprintf(file,"#\n");

	  for (i=0; i<mf.NBIN; i++)
	    {
	      x=mf.mmin+(i+0.5)*DELTAM;
	      m=pow(10.,x);
	      dm = params.ParticleMass * 
		(double)( (int) ( pow(10.,mf.mmin+(i+1)*DELTAM)/params.ParticleMass ) -
			  (int) ( pow(10.,mf.mmin+ i   *DELTAM)/params.ParticleMass ) );
	      if (dm>0.0)
		{
		  a=(double)mf.ninbin[i]/mf.vol/dm;
		  a1=((double)mf.ninbin[i]+sqrt((double)mf.ninbin[i]))/mf.vol/dm;
		  a2=((double)mf.ninbin[i]-sqrt((double)mf.ninbin[i]))/mf.vol/dm;
		}
	      else
		a=a1=a2=0.0;

	      if (mf.ninbin[i] > 1)
		mx=mf.massinbin[i]/(double)mf.ninbin[i];
	      else
		mx=m;
	      a3=AnalyticMassFunction(mx,outputs.z[iout]);

	      r=SizeForMass(mx);
	      massvar=MassVariance(r)*D*D;
	      sigma=sqrt(massvar);
	      ni=1.686/sigma;

	      fprintf(file,
		      " %15.8g %15.8g %15.8g %15.8g   %10d  %15.8g    %15.8g\n",
		      mx*mf.hfactor,
		      a/mf.hfactor4,
		      a1/mf.hfactor4,
		      a2/mf.hfactor4,
		      mf.ninbin[i],
		      a3/mf.hfactor4,
		      ni);
	    }
	  fclose(file);
	}
    }

  /* Bye! */
  return 0;
}


/* 
   GENERAL COMMUNICATION SCHEME

   The code wants to write a catalog splitted in params.NumFiles
   Each file will be written by this number of tasks
      NTasksPerFile=NTasks/params.NumFiles;

   Task ThisTask will write on file 
      ThisFile=ThisTask/NTasksPerFile;

   Among the tasks that access the file, one task will collect information 
   and write it on the file. This is:
      collector=ThisFile*NTasksPerFile;

   1: initialization of whatever

   2: each task constructs its own catalog (whatever it is)

   3: collector task opens the file and writes the header

   4: collector task writes the catalog on the file

   5: loop on all other tasks that write on the same file

     5a: the task sends its catalog to collector task
     5b: collector task receives the catalog

     5c: collector task writes the catalog on the file

   6: collector task closes the file

*/


int write_catalog(int iout)
{
  /* Writes the group catalogues */

  int igood,i,ngood,nhalos,j;
  double hfactor,GGrid[3],SGrid[3];
  char filename[LBLENGTH],labh[3];
  int NTasksPerFile,collector,itask,next,ThisFile;
  catalog_data *mycat;
  FILE *file;
  MPI_Status status;
#ifdef FORTRAN
  int idummy;
#endif

  /* ordering of coordinates to accomodate for rotation caused by fft ordering */
#ifdef ROTATE_BOX
  static int rot[3]={1,2,0};
#else
  static int rot[3]={0,1,2};
#endif

  if (params.DoNotWriteCatalogs)
    {
      if (!ThisTask)
	printf("Halo catalog at z=%f will not be written\n",outputs.z[iout]);
      return 0;
    }

  GGrid[0]=(double)MyGrids[0].GSglobal_x;
  GGrid[1]=(double)MyGrids[0].GSglobal_y;
  GGrid[2]=(double)MyGrids[0].GSglobal_z;
  SGrid[0]=(double)subbox.stabl_x;
  SGrid[1]=(double)subbox.stabl_y;
  SGrid[2]=(double)subbox.stabl_z;

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  /* output in H100 or in Htrue */
  if (params.OutputInH100)
    hfactor=params.Hubble100;
  else
    hfactor=1.0;

  /* each processor builds the catalogue */
  nhalos=0;
  for (i=FILAMENT+1, ngood=0; i<=ngroups; i++)
     if (groups[i].point > 0 && groups[i].good && 
	 groups[i].Mass >= params.MinHaloMass) 
       ngood++;

  /* space to store catalogs */
  mycat = (catalog_data *)wheretoplace_mycat;
  if (ngood * sizeof(catalog_data) > subbox.Npart/10*sizeof(histories_data))
    {
      printf("ERROR on task %d: surprisingly, memory reserved to mycat is insufficient in write_catalog\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  if (ngood)
    {
      for (i=FILAMENT+1, igood=0; i<=ngroups; i++)
	if (groups[i].point > 0 && groups[i].good && 
	    groups[i].Mass >= params.MinHaloMass) 
	  {
	    set_obj(i,outputs.F[iout],&obj1);
	    set_obj_vel(i,outputs.F[iout],&obj1);
	    mycat[igood].name=groups[i].name;
	    mycat[igood].n=groups[i].Mass;
	    mycat[igood].M=groups[i].Mass*params.ParticleMass*hfactor;
	    for (j=0; j<3; j++)
	      {
		/* q and x are in sub-box coordinates, they are transformed
		   to global box coordinates (forcing PBCs) */
		mycat[igood].q[j] = groups[i].Pos[rot[j]] + SGrid[rot[j]];
		if (mycat[igood].q[j]<0)
		  mycat[igood].q[j]+=GGrid[j];
		if (mycat[igood].q[j]>=GGrid[j])
		  mycat[igood].q[j]-=GGrid[j];
		mycat[igood].q[j] *= params.InterPartDist*hfactor;
		/* displacement is done up to ORDER_FOR_CATALOG */
		mycat[igood].x[j] = q2x(rot[j],&obj1,ORDER_FOR_CATALOG) + SGrid[rot[j]];
		if (mycat[igood].x[j]<0)
		  mycat[igood].x[j]+=GGrid[j];
		if (mycat[igood].x[j]>=GGrid[j])
		  mycat[igood].x[j]-=GGrid[j];
		mycat[igood].x[j]*=params.InterPartDist*hfactor;
		mycat[igood].v[j] = vel(rot[j],&obj1);
	      }
	    igood++;
	  }
    }

  /* The collector task opens the file and writes its catalogue */
  if (ThisTask==collector)
    {
      if (params.NumFiles>1)
	sprintf(filename,"pinocchio.%6.4f.%s.catalog.out.%d",
		outputs.z[iout],params.RunFlag,ThisFile);
      else
	sprintf(filename,"pinocchio.%6.4f.%s.catalog.out",
		outputs.z[iout],params.RunFlag);

      if (!ThisTask)
	{
	  if (ThisSlice)
	    printf("[%s] Updating file %s\n",fdate(),filename);
	  else
	    printf("[%s] Opening file %s\n",fdate(),filename);
	}

      if (ThisSlice)
	file=fopen(filename,"a");
      else
	file=fopen(filename,"w");

      if (params.CatalogInAscii)
	{
	  if (!ThisSlice)
	    {
	      fprintf(file,"# Group catalog for redshift %f and minimal mass of %d particle%s\n",
		      outputs.z[iout],params.MinHaloMass,(params.MinHaloMass==1?"":"s"));

	      if (params.OutputInH100)
		strcpy(labh,"/h");
	      else
		strcpy(labh,"");

	      fprintf(file,"#    1) group ID\n");
	      fprintf(file,"#    2) group mass (Msun%s)\n",labh);
	      fprintf(file,"# 3- 5) initial position (Mpc%s)\n",labh);
	      fprintf(file,"# 6- 8) final position (Mpc%s)\n",labh);
	      fprintf(file,"# 9-11) velocity (km/s)\n");
	      fprintf(file,"#   12) number of particles\n");
	      fprintf(file,"#\n");
	    }
	  if (ngood)
	    for (igood=0; igood<ngood; igood++)
	      fprintf(file," %12Lu %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
		      mycat[igood].name,
		      mycat[igood].M,
		      mycat[igood].q[0],mycat[igood].q[1],mycat[igood].q[2],
		      mycat[igood].x[0],mycat[igood].x[1],mycat[igood].x[2],
		      mycat[igood].v[0],mycat[igood].v[1],mycat[igood].v[2],
		      mycat[igood].n);
	}
      else
	{
	  if (!ThisSlice)
	    {
#ifdef FORTRAN
	      idummy=2*sizeof(int);
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	      fwrite(&NTasksPerFile,sizeof(int),1,file);
	      fwrite(&NSlices,sizeof(int),1,file);
#ifdef FORTRAN
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	    }
#ifdef FORTRAN
	  idummy=sizeof(int);
	  fwrite(&idummy,sizeof(int),1,file);
#endif
	  fwrite(&ngood,sizeof(int),1,file);
#ifdef FORTRAN
	  fwrite(&idummy,sizeof(int),1,file);
#endif

#ifdef FORTRAN
	  idummy=sizeof(catalog_data);
#endif
	  if (ngood)
	    for (igood=0; igood<ngood; igood++)
	      {
#ifdef FORTRAN
		fwrite(&idummy,sizeof(int),1,file);
#endif
		fwrite(mycat+igood,sizeof(catalog_data),1,file);
#ifdef FORTRAN
		fwrite(&idummy,sizeof(int),1,file);
#endif
	      }
	}

      nhalos+=ngood;
    }

  /* loop on other tasks that must give data to collector

     itask sends the number of good groups, collector receives it
     then itask sends the groups to collector (and deallocates mycat)
     while collector (allocates mycat and) receives the groups;
     finally collector writes mycat on the file 
  */
  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;

      if (ThisTask==collector)
	{
	  ngood=0;
	  MPI_Recv(&ngood, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);

	  if (!params.CatalogInAscii)
	    {
#ifdef FORTRAN
	      idummy=sizeof(int);
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	      fwrite(&ngood,sizeof(int),1,file);
#ifdef FORTRAN
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	    }

	  if (ngood)
	    {
	      MPI_Recv(mycat, ngood*sizeof(catalog_data), MPI_CHAR, itask, 0, MPI_COMM_WORLD, &status);
	    
	      if (params.CatalogInAscii)
		{
		  for (igood=0; igood<ngood; igood++)
		    fprintf(file," %12Lu %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
			    mycat[igood].name,
			    mycat[igood].M,
			    mycat[igood].q[0],mycat[igood].q[1],mycat[igood].q[2],
			    mycat[igood].x[0],mycat[igood].x[1],mycat[igood].x[2],
			    mycat[igood].v[0],mycat[igood].v[1],mycat[igood].v[2],
			    mycat[igood].n);
		}
	      else
		{
#ifdef FORTRAN
		  idummy=sizeof(catalog_data);
#endif
		  for (igood=0; igood<ngood; igood++)
		    {
#ifdef FORTRAN
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      fwrite(mycat+igood,sizeof(catalog_data),1,file);
#ifdef FORTRAN
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		    }
		}
	      
	      nhalos+=ngood;
	    }
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&ngood, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  if (ngood)
	    {
	      MPI_Send(mycat, ngood*sizeof(catalog_data), MPI_CHAR, collector, 0, MPI_COMM_WORLD);
	    }
	}

    }
  /* end of loop on tasks */
  if (ThisTask==collector) 
    {
      fclose(file);
      printf("[%s] Task %d has written %d halos on file %s\n",fdate(),ThisTask,nhalos,filename);
      fflush(stdout);
    }


  return 0;
}

#ifdef PLC
int write_PLC()
{
  /* Writes the plc catalogues */

  int i,nhalos,nstored;
  double hfactor;
  char filename[LBLENGTH],labh[3];
  int NTasksPerFile,collector,itask,next,ThisFile;
  FILE *file;
  MPI_Status status;
  static int FirstCall=1;
  struct towrite
  {
    unsigned long long int name;
    double red,x,y,z,vx,vy,vz,Mass,theta,phi,v_los,obsz;
  } writethis;
#ifdef FORTRAN
  int idummy;
#endif


  if (params.DoNotWriteCatalogs)
    {
      if (!ThisTask && FirstCall)
	printf("PLC catalog will not be written\n");
      return 0;
    }

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  /* output in H100 or in Htrue */
  if (params.OutputInH100)
    hfactor=params.Hubble100;
  else
    hfactor=1.0;
  nhalos=0;

  /* The collector task opens the file and writes its catalogue */
  if (ThisTask==collector)
    {
      if (params.NumFiles>1)
	sprintf(filename,"pinocchio.%s.plc.out.%d",params.RunFlag,ThisFile);
      else
	sprintf(filename,"pinocchio.%s.plc.out",params.RunFlag);

      if (!ThisTask)
	{
	  if (FirstCall)
	    printf("[%s] Opening file %s\n",fdate(),filename);
	  else
	    printf("[%s] Updating file %s\n",fdate(),filename);
	}

      if (FirstCall)
	file=fopen(filename,"w");
      else
	file=fopen(filename,"a");

      if (params.CatalogInAscii)
	{
	  if (FirstCall)
	    {
	      fprintf(file,"# Group catalog on the Past Light Cone for a minimal mass of %d particle%s\n",
		      params.MinHaloMass,(params.MinHaloMass==1?"":"s"));

	      if (params.OutputInH100)
		strcpy(labh,"/h");
	      else
		strcpy(labh,"");

	      fprintf(file,"#    1) group ID\n");
	      fprintf(file,"#    2) true redshift\n");
	      fprintf(file,"#  3-5) comoving position (Mpc%s)\n",labh);
	      fprintf(file,"#  6-8) velocity (km/s)\n");
	      fprintf(file,"#    9) group mass (Msun%s)\n",labh);
	      fprintf(file,"#   10) theta (degree)\n");
	      fprintf(file,"#   11) phi (degree)\n");
	      fprintf(file,"#   12) peculiar velocity along the line-of-sight (km/s)\n");
	      fprintf(file,"#   13) observed redshift\n");
	      fprintf(file,"#\n");
	    }

	  for (i=0; i<plc.Nstored; i++)
	    {
	      writethis.name=plcgroups[i].name;
	      writethis.red=plcgroups[i].z;
	      writethis.x=plcgroups[i].x[0]*hfactor;
	      writethis.y=plcgroups[i].x[1]*hfactor;
	      writethis.z=plcgroups[i].x[2]*hfactor;
	      writethis.vx=plcgroups[i].v[0];
	      writethis.vy=plcgroups[i].v[1];
	      writethis.vz=plcgroups[i].v[2];
	      writethis.Mass=plcgroups[i].Mass*params.ParticleMass*hfactor;
	      writethis.theta=plcgroups[i].theta;
	      writethis.phi=plcgroups[i].phi;
              writethis.v_los=(plcgroups[i].x[0]*plcgroups[i].v[0]+
			       plcgroups[i].x[1]*plcgroups[i].v[1]+
			       plcgroups[i].x[2]*plcgroups[i].v[2])/plcgroups[i].rhor;
              writethis.obsz=plcgroups[i].z+writethis.v_los/SPEEDOFLIGHT*(1.0+plcgroups[i].z);

	      fprintf(file," %12Lu %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %15.8e %16.6f %16.6f %16.6f %16.6f\n",
		      writethis.name,
		      writethis.red,
		      writethis.x,writethis.y,writethis.z,
		      writethis.vx,writethis.vy,writethis.vz,
		      writethis.Mass,
		      writethis.theta,
		      writethis.phi,
		      writethis.v_los,
		      writethis.obsz);
	      nhalos++;

	    }
	}
      else
	{
	  for (i=0; i<plc.Nstored; i++)
	    {
	      writethis.name=plcgroups[i].name;
	      writethis.red=plcgroups[i].z;
	      writethis.x=plcgroups[i].x[0]*hfactor;
	      writethis.y=plcgroups[i].x[1]*hfactor;
	      writethis.z=plcgroups[i].x[2]*hfactor;
	      writethis.vx=plcgroups[i].v[0];
	      writethis.vy=plcgroups[i].v[1];
	      writethis.vz=plcgroups[i].v[2];
	      writethis.Mass=plcgroups[i].Mass*params.ParticleMass*hfactor;
	      writethis.theta=plcgroups[i].theta;
	      writethis.phi=plcgroups[i].phi;
              writethis.v_los=(plcgroups[i].x[0]*plcgroups[i].v[0]+
			       plcgroups[i].x[1]*plcgroups[i].v[1]+
			       plcgroups[i].x[2]*plcgroups[i].v[2])/plcgroups[i].rhor;
              writethis.obsz=plcgroups[i].z+writethis.v_los/SPEEDOFLIGHT*(1.0+plcgroups[i].z);
#ifdef FORTRAN
	      idummy=sizeof(struct towrite);
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	      fwrite(&writethis,sizeof(struct towrite),1,file);
#ifdef FORTRAN
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	      nhalos++;
	    }
	}
    }

  /* loop on other tasks that collect data on collector

     itask sends the number of good groups, collector receives it
     then itask sends the groups to collector (and deallocates mycat)
     while collector (allocates mycat and) receives the groups;
     finally collector writes mycat on the file 
  */
  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;

      if (ThisTask==collector)
	{
	  nstored=0;
	  MPI_Recv(&nstored, 1, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);

	  if (nstored)
	    {

	      MPI_Recv(plcgroups, nstored*sizeof(plcgroup_data), MPI_CHAR, itask, 0, MPI_COMM_WORLD, &status);

	      if (params.CatalogInAscii)
		{
		  for (i=0; i<nstored; i++)
		    {
		      writethis.name=plcgroups[i].name;
		      writethis.red=plcgroups[i].z;
		      writethis.x=plcgroups[i].x[0]*hfactor;
		      writethis.y=plcgroups[i].x[1]*hfactor;
		      writethis.z=plcgroups[i].x[2]*hfactor;
		      writethis.vx=plcgroups[i].v[0];
		      writethis.vy=plcgroups[i].v[1];
		      writethis.vz=plcgroups[i].v[2];
		      writethis.Mass=plcgroups[i].Mass*params.ParticleMass*hfactor;
		      writethis.theta=plcgroups[i].theta;
		      writethis.phi=plcgroups[i].phi;
		      writethis.v_los=(plcgroups[i].x[0]*plcgroups[i].v[0]+
				       plcgroups[i].x[1]*plcgroups[i].v[1]+
				       plcgroups[i].x[2]*plcgroups[i].v[2])/plcgroups[i].rhor;
		      writethis.obsz=plcgroups[i].z+writethis.v_los/SPEEDOFLIGHT*(1.0+plcgroups[i].z);

		      fprintf(file," %12Lu %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %15.8e %16.6f %16.6f %16.6f %16.6f\n",
			      writethis.name,
			      writethis.red,
			      writethis.x,writethis.y,writethis.z,
			      writethis.vx,writethis.vy,writethis.vz,
			      writethis.Mass,
			      writethis.theta,
			      writethis.phi,
			      writethis.v_los,
			      writethis.obsz);
		      nhalos++;
		    }
		}
	      else
		{
		  for (i=0; i<nstored; i++)
		    {
		      writethis.name=plcgroups[i].name;
		      writethis.red=plcgroups[i].z;
		      writethis.x=plcgroups[i].x[0]*hfactor;
		      writethis.y=plcgroups[i].x[1]*hfactor;
		      writethis.z=plcgroups[i].x[2]*hfactor;
		      writethis.vx=plcgroups[i].v[0];
		      writethis.vy=plcgroups[i].v[1];
		      writethis.vz=plcgroups[i].v[2];
		      writethis.Mass=plcgroups[i].Mass*params.ParticleMass*hfactor;
		      writethis.theta=plcgroups[i].theta;
		      writethis.phi=plcgroups[i].phi;
		      writethis.v_los=(plcgroups[i].x[0]*plcgroups[i].v[0]+
				       plcgroups[i].x[1]*plcgroups[i].v[1]+
				       plcgroups[i].x[2]*plcgroups[i].v[2])/plcgroups[i].rhor;
		      writethis.obsz=plcgroups[i].z+writethis.v_los/SPEEDOFLIGHT*(1.0+plcgroups[i].z);

#ifdef FORTRAN
		      idummy=sizeof(struct towrite);
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      fwrite(&writethis,sizeof(struct towrite),1,file);
#ifdef FORTRAN
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      nhalos++;
		    }
		}
	    }
	}
      else if (ThisTask==itask)
	{
	  MPI_Send(&plc.Nstored, 1, MPI_INT, collector, 0, MPI_COMM_WORLD);
	  if (plc.Nstored)
	    {
	      MPI_Send(plcgroups, plc.Nstored*sizeof(plcgroup_data), MPI_CHAR, collector, 0, MPI_COMM_WORLD);
	    }
	}
      
    }
  /* end of loop on tasks */
  if (ThisTask==collector) 
    {
      fclose(file);
      printf("[%s] Task %d has written %d halos on file %s\n",fdate(),ThisTask,nhalos,filename);
      fflush(stdout);
    }

  FirstCall=0;

  return 0;
}
#endif


int write_histories(void)
{
  /* writes merger histories at the end of the run */

  int i, commint[2], ntrees, nbranch_all, nbranch_tree, ntrees_global, nbranch_global,
    ibranch, thisbranch, thistree, itree;
  char filename[LBLENGTH];
  int NTasksPerFile,collector,itask,next,ThisFile;
  histories_data *mycat;
  FILE *file;
  MPI_Status status;
#ifdef FORTRAN
  int idummy;
#endif

  if (params.DoNotWriteHistories)
    {
      if (!ThisTask)
	printf("Merger histories will not be written\n");
      return 0;
    }

  NTasksPerFile=NTasks/params.NumFiles;
  ThisFile=ThisTask/NTasksPerFile;
  collector=ThisFile*NTasksPerFile;

  for (i=0; i<=ngroups; i++)
    groups[i].trackT=0;
  for (i=0; i<=ngroups; i++)
    groups[i].trackC=0;


  /* each processor builds its catalogue */

  /* this loop counts the total number of branches in the local catalog */
  nbranch_all=ntrees=0;
  for (i=FILAMENT+1, nbranch_all=0;  i<=ngroups; i++)
    if (groups[groups[i].halo_app].good &&                       /* if the halo belongs to a good group */
	groups[groups[i].halo_app].Mass >= params.MinHaloMass)   /* and if its mass is large enough */
      nbranch_all++;

  mycat = (histories_data *)wheretoplace_mycat;
  if (nbranch_all * sizeof(histories_data) > subbox.Npart/10*sizeof(histories_data))
    {
      printf("ERROR on task %d: surprisingly, memory reserved to mycat is insufficient in write_histories\n",ThisTask);
      fflush(stdout);
      return 1;
    }

  /* assign the catalog by looping on all the trees */
  if (nbranch_all)
    {
      ntrees=0;      /* number of trees (of main halos at the final redshift) */
      thisbranch=0;  /* this is the counter of the number of stored halos (branches) */

      /* loop on main halos */
      for (i=FILAMENT+1, ntrees=0;  i<=ngroups; i++)
	if (groups[i].halo_app==i &&               /* if this is a main halo */
	    groups[i].good &&                      /* and if it is a good one */
	    groups[i].Mass >= params.MinHaloMass)  /* and if its mass is large enough */
	  {
	    ++ntrees;

	    /* counts the number of branches in the tree */
	    nbranch_tree=0;
	    next=i;
	    do
	      {
		++nbranch_tree;
		next=groups[next].ll;  /* next branch in the tree */
	      }
	    while (next != i);

	    ibranch=0;     /* this is the counter of branches within the tree, starting from zero
			      but zero is assigned to the main halo, whose nick will be nbranch_tree */
	    next=i;        /* this is the pointer to the halo in the main catalog */
	    /* NB: we start from the main halo */

	    do
	      {
		/* track vector keeps track of the relation between main catalog and numbering within the single tree */
		groups[thisbranch].trackT = next;
		mycat[thisbranch].nick = groups[next].trackC = (ibranch ? ibranch : nbranch_tree);
		mycat[thisbranch].ll   = ++ibranch;
		mycat[thisbranch].mass = groups[next].Mass;
		mycat[thisbranch].name = groups[next].name;
		mycat[thisbranch].mam  = groups[next].mass_at_merger;
		mycat[thisbranch].zap  = FTOZ(groups[next].t_appear);
		mycat[thisbranch].zpe  = FTOZ(groups[next].t_peak);
		mycat[thisbranch].zme  = FTOZ(groups[next].t_merge);
		++thisbranch;

		next=groups[next].ll;  /* next branch in the tree */
	      }  /* this loop closes when we get to the original halo, that is stored as the last one */
	    while (next != i);

	    /* now it fixes all the merged_with fields */
	    for (ibranch=thisbranch-nbranch_tree; ibranch<thisbranch; ibranch++)
	      if (groups[groups[ibranch].trackT].merged_with>FILAMENT)
		mycat[ibranch].mw = groups[groups[groups[ibranch].trackT].merged_with].trackC;
	      else
		mycat[ibranch].mw = -1;
	  }
    }

  /* computes the total number of trees and branches in the file */
  /* NB: this could be done by a reduce, but redefining the communicator is boring */
  if (ThisTask==collector)
    {
      ntrees_global=ntrees;
      nbranch_global=nbranch_all;
    }
  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;
      if (ThisTask==collector)
	{
	  MPI_Recv(commint, 2, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  ntrees_global+=commint[0];
	  nbranch_global+=commint[1];
	}
      else if (ThisTask==itask)
	{
	  commint[0]=ntrees;
	  commint[1]=nbranch_all;
	  MPI_Send(commint, 2, MPI_INT, collector, 0, MPI_COMM_WORLD);
	}
    }

  /* collector task opens the file and writes the header */
  if (ThisTask==collector)
    {
      if (params.NumFiles>1)
	sprintf(filename,"pinocchio.%s.histories.out.%d",params.RunFlag,ThisFile);
      else
	sprintf(filename,"pinocchio.%s.histories.out",params.RunFlag);

      if (!ThisTask)
	{
	  if (ThisSlice)
	    printf("[%s] Updating file %s\n",fdate(),filename);
	  else
	    printf("[%s] Opening file %s\n",fdate(),filename);
	}

      if (ThisSlice)
	file=fopen(filename,"a");
      else
	file=fopen(filename,"w");
      thistree=0;

      if (params.CatalogInAscii)
	{
	  if (!ThisSlice)
	    {
	      fprintf(file,"# Merger histories for halos with minimal mass of %d particle%s\n",
		      params.MinHaloMass,(params.MinHaloMass==1?"":"s"));

	      fprintf(file,"#  1) group ID\n");
	      fprintf(file,"#  2) index within the tree\n");
	      fprintf(file,"#  3) linking list\n");
	      fprintf(file,"#  4) merged with\n");
	      fprintf(file,"#  5) mass of halo at merger (particles)\n");
	      fprintf(file,"#  6) mass of main halo it merges with, at merger (particles)\n");
	      fprintf(file,"#  7) merger redshift\n");
	      fprintf(file,"#  8) redshift of peak collapse\n");
	      fprintf(file,"#  9) redshift at which the halo overtakes the minimal mass\n");
	      fprintf(file,"#\n");
	      fprintf(file,"# Nslices: \n");
	      fprintf(file," %d\n",NSlices);
	    }	
	  fprintf(file,"# Ntrees & Nbranches: \n");
	  fprintf(file," %d  %d\n",ntrees_global, nbranch_global);
	}
      else
	{
	  if (!ThisSlice)
	    {
#ifdef FORTRAN
	      idummy=sizeof(int);
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	      fwrite(&NSlices,sizeof(int),1,file);
#ifdef FORTRAN
	      fwrite(&idummy,sizeof(int),1,file);
#endif
	    }
#ifdef FORTRAN
	  idummy=2*sizeof(int);
	  fwrite(&idummy,sizeof(int),1,file);
#endif
	  fwrite(&ntrees_global,sizeof(int),1,file);
	  fwrite(&nbranch_global,sizeof(int),1,file);
#ifdef FORTRAN
	  fwrite(&idummy,sizeof(int),1,file);
#endif
	}
    }

  /* collector writes its catalog */
  if (ThisTask==collector)
    {
      if (nbranch_all)
	{
	  thisbranch=0;
	  if (params.CatalogInAscii)
	    {
	      for (itree=0; itree<ntrees; itree++)
		{
		  nbranch_tree=mycat[thisbranch].nick;
		  fprintf(file,"#Tree %d, Nbranches=%d\n", thistree++, nbranch_tree);
		  for (ibranch=0; ibranch<nbranch_tree; ibranch++)
		    {
		      fprintf(file," %12Ld %6d %6d %6d %9d %9d %9.4f %9.4f %9.4f\n",
			      mycat[thisbranch].name,
			      mycat[thisbranch].nick,
			      mycat[thisbranch].ll,
			      mycat[thisbranch].mw,
			      mycat[thisbranch].mass,
			      mycat[thisbranch].mam,
			      mycat[thisbranch].zme,
			      mycat[thisbranch].zpe,
			      mycat[thisbranch].zap);
		      ++thisbranch;
		    }
		}
	    }
	  else
	    {
	      for (itree=0; itree<ntrees; itree++)
		{
		  nbranch_tree=mycat[thisbranch].nick;
#ifdef FORTRAN
		  idummy=2*sizeof(int);
		  fwrite(&idummy,sizeof(int),1,file);
#endif
		  fwrite(&thistree,sizeof(int),1,file);
		  fwrite(&nbranch_tree,sizeof(int),1,file);
		  ++thistree;
#ifdef FORTRAN
		  fwrite(&idummy,sizeof(int),1,file);
#endif
		  for (ibranch=0; ibranch<nbranch_tree; ibranch++)
		    {
#ifdef FORTRAN
		      idummy=sizeof(histories_data);
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      fwrite(mycat+thisbranch++,sizeof(histories_data),1,file);
#ifdef FORTRAN
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		    }
		}
	    }
	}
    }

  /* loop on other tasks that must give data to collector

     itask sends the number of good groups, collector receives it
     then itask sends the groups to collector (and deallocates mycat)
     while collector (allocates mycat and) receives the groups;
     finally collector writes mycat on the file 
  */

  for (next=1; next<NTasksPerFile; next++)
    {
      itask=collector+next;

      if (ThisTask==collector)
	{
	  MPI_Recv(commint, 2, MPI_INT, itask, 0, MPI_COMM_WORLD, &status);
	  nbranch_all=commint[0];
	  ntrees=commint[1];

	  if (nbranch_all)
	    {
	      MPI_Recv(mycat, nbranch_all*sizeof(histories_data), MPI_CHAR, itask, 0, MPI_COMM_WORLD, &status);

	      /* writes the catalog that has just received */
	      thisbranch=0;
	      if (params.CatalogInAscii)
		{
		  for (itree=0; itree<ntrees; itree++)
		    {
		      nbranch_tree=mycat[thisbranch].nick;
		      fprintf(file,"#Tree %d, Nbranches=%d\n", thistree++, nbranch_tree);
		      for (ibranch=0; ibranch<nbranch_tree; ibranch++)
			{
			  fprintf(file," %12Ld %6d %6d %6d %9d %9d %9.4f %9.4f %9.4f\n",
				  mycat[thisbranch].name,
				  mycat[thisbranch].nick,
				  mycat[thisbranch].ll,
				  mycat[thisbranch].mw,
				  mycat[thisbranch].mass,
				  mycat[thisbranch].mam,
				  mycat[thisbranch].zme,
				  mycat[thisbranch].zpe,
				  mycat[thisbranch].zap);
			  ++thisbranch;
			}
		    }
		}
	      else
		{
		  for (itree=0; itree<ntrees; itree++)
		    {
		      nbranch_tree=mycat[thisbranch].nick;
#ifdef FORTRAN
		      idummy=2*sizeof(int);
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      fwrite(&thistree,sizeof(int),1,file);
		      fwrite(&nbranch_tree,sizeof(int),1,file);
		      ++thistree;
#ifdef FORTRAN
		      fwrite(&idummy,sizeof(int),1,file);
#endif
		      for (ibranch=0; ibranch<nbranch_tree; ibranch++)
			{
#ifdef FORTRAN
			  idummy=sizeof(histories_data);
			  fwrite(&idummy,sizeof(int),1,file);
#endif
			  fwrite(mycat+thisbranch++,sizeof(histories_data),1,file);
#ifdef FORTRAN
			  fwrite(&idummy,sizeof(int),1,file);
#endif
			}
		    }
		}
	    }
	}
      else if (ThisTask==itask)
	{
	  /* itask sends its data to collector */
	  commint[0]=nbranch_all;
	  commint[1]=ntrees;
	  MPI_Send(commint, 2, MPI_INT, collector, 0, MPI_COMM_WORLD);

	  if (nbranch_all)
	    {
	      MPI_Send(mycat, nbranch_all*sizeof(histories_data), MPI_CHAR, collector, 0, MPI_COMM_WORLD);
	    }
	}

    }

  /* end of loop on tasks */
  if (ThisTask==collector) 
    {
      fclose(file);
      printf("[%s] Task %d has written %d trees and %d branches on file %s\n",fdate(),ThisTask,ntrees_global,nbranch_global,filename);
      fflush(stdout);
    }

  return 0;
}
