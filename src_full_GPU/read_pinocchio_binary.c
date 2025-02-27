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


/* 
   
   This code reads a catalog written by pinocchio in binary format 
   and checks its integrity
   If required, it writes the catalog in ascii format

 */

/* to compile the code: */
/*    gcc -O -o read_pinocchio_binary read_pinocchio_binary.c  */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int main(int argc, char **argv)
{
  typedef struct
  {
    unsigned long long int id;
    double M;
    double q[3], x[3], v[3];
    int N, pad;
  } cat;
  cat catdata;

  typedef struct
  {
    unsigned long long int id;
    double z;
    double x[3], v[3];
    double M, th, ph, vl, zo;
  } plc;
  plc plcdata;

  typedef struct
  {
    unsigned long long int name;
    int nick, ll, mw, mass, mam;
    double zme, zpe, zap;
  }  hist;
  hist histdata;


  unsigned long long int id,n,ntrees_tot,nbranch_tot;
  int nproc,iproc,ngood,igood,i,dummy,ninhead,NSlices,ThisSlice,ntrees,nbranch,
    nbranch_tree,itree,ibranch,thistree,u,nblocks,print=0,First=1;
  char runname[100],filename[100],redshift[100];
  FILE *fd;
  FILE *fout;
  char outputf[100];


  printf("This is an example code to read pinocchio catalogs in binary format\n");
  printf("Pinocchio generates catalogs at fixed redshift in files\n");
  printf("like pinocchio.0.0000.SecondEuclidMocks0001.catalog.out,\n");
  printf("where SecondEuclidMocks0001 is the run name\n");
  printf("and files like pinocchio.SecondEuclidMocks0001.plc.out\n");
  printf("for the catalog on the past-light-cone\n");
  printf("or pinocchio.SecondEuclidMocks0001.histories.out\n");
  printf("for the halo merger trees\n");
  printf("\n");
  printf("Usage: read_pinocchio_binary runname redshift [ascii]\n");
  printf("     runname is the name of the run for all the files\n");
  printf("     redshift must be equal to the redshift given in the file name\n");
  printf("     to read the past-light-cone enter the string 'plc' in place of redshift\n");
  printf("     to read the histories file enter string 'hist' in place of redshift\n");
  printf("\n");
  printf("To print an ascii version of the file, add the string ascii as a further argument\n");
  printf("but BEWARE, the number of halos in the catalog can be large and ascii is not\n");
  printf("a recommended format to store catalogs\n\n");

  if (argc<3)
    return 0;

  strcpy(runname,argv[1]);
  strcpy(redshift,argv[2]);
  if (argc>=4 && !strcmp(argv[3],"ascii"))
    print=1;

  if (!strcmp(redshift,"plc"))
    /* CATALOG ON THE PAST LIGHT CONE */
    {
      sprintf(filename,"pinocchio.%s.plc.out",runname);
      printf("Opening file %s\n",filename);
      printf("I will write out the first ten halos\n");

      fd=fopen(filename,"r");
      if (fd==0x0)
	{
	  printf("Error: plc file %s not found\n",filename);
	  return 1;
	}

      if (print)
	{
	  sprintf(outputf,"pinocchio.%s.plc.out.ascii",runname);
	  printf("I will write the ascii catalog in the file %s\n",outputf);
	  fout=fopen(outputf,"w");
	  fprintf(fout,"# Group catalog on the Past Light Cone\n");
	  fprintf(fout,"#    1) group ID\n");
	  fprintf(fout,"#    2) true redshift\n");
	  fprintf(fout,"#  3-5) comoving position\n");
	  fprintf(fout,"#  6-8) velocity\n");
	  fprintf(fout,"#    9) group mass\n");
	  fprintf(fout,"#   10) theta\n");
	  fprintf(fout,"#   11) phi\n");
	  fprintf(fout,"#   12) peculiar velocity along the line-of-sight\n");
	  fprintf(fout,"#   13) observed redshift\n");
	  fprintf(fout,"#\n");
	}

      n=0;
      while (1)
	{
	  if (!fread(&dummy,sizeof(int),1,fd))
	    break;
	  u=fread(&plcdata,dummy,1,fd);
	  u=fread(&dummy,sizeof(int),1,fd);

	  if (n<10)
	    printf(" %12Ld %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %15.8e %16.6f %16.6f %16.6f %16.6f\n",
		   plcdata.id,plcdata.z,
		   plcdata.x[0], plcdata.x[1], plcdata.x[2],
		   plcdata.v[0], plcdata.v[1], plcdata.v[2],
		   plcdata.M, plcdata.th, plcdata.ph, plcdata.vl, plcdata.zo);
	  n++;

	  if (print)
	    fprintf(fout," %12Lu %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %15.8e %16.6f %16.6f %16.6f %16.6f\n",
		    plcdata.id,
		    plcdata.z,
		    plcdata.x[0],plcdata.x[1],plcdata.x[2],
		    plcdata.v[0],plcdata.v[1],plcdata.v[2],
		    plcdata.M,
		    plcdata.th,
		    plcdata.ph,
		    plcdata.vl,
		    plcdata.zo);

	}
      fclose(fd);
      if (print)
      fclose(fout);


      printf("I found %Ld halos in the file\n",n);
      printf("This is the last halo in the catalog:\n");
      printf(" %12Ld %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %16.6f %15.8e %16.6f %16.6f %16.6f %16.6f\n",
      	     plcdata.id,plcdata.z,
      	     plcdata.x[0], plcdata.x[1], plcdata.x[2],
      	     plcdata.v[0], plcdata.v[1], plcdata.v[2],
      	     plcdata.M, plcdata.th, plcdata.ph, plcdata.vl, plcdata.zo);
    }
  else if (!strcmp(redshift,"hist"))
    /* MERGER HISTORIES */
    {
      sprintf(filename,"pinocchio.%s.histories.out",runname);
      printf("Opening file %s\n",filename);
      printf("I will write out the first tree\n");

      fd=fopen(filename,"r");
      if (fd==0x0)
	{
	  printf("Error: history file %s not found\n",filename);
	  return 1;
	}

      if (print)
	{
	  sprintf(outputf,"pinocchio.%s.histories.out.ascii",runname);
	  printf("I will write the ascii catalog in the file %s\n",outputf);
	  fout=fopen(outputf,"w");

	  fprintf(fout,"# Merger histories\n");
	  fprintf(fout,"#  1) group ID\n");
	  fprintf(fout,"#  2) index within the tree\n");
	  fprintf(fout,"#  3) linking list\n");
	  fprintf(fout,"#  4) merged with\n");
	  fprintf(fout,"#  5) mass of halo at merger (particles)\n");
	  fprintf(fout,"#  6) mass of main halo it merges with, at merger (particles)\n");
	  fprintf(fout,"#  7) merger redshift\n");
	  fprintf(fout,"#  8) redshift of peak collapse\n");
	  fprintf(fout,"#  9) redshift at which the halo overtakes the minimal mass\n");
	  fprintf(fout,"#\n");
	}

      u=fread(&dummy,sizeof(int),1,fd);
      u=fread(&NSlices,sizeof(int),1,fd);
      u=fread(&dummy,sizeof(int),1,fd);

      if (print)
	{
	  fprintf(fout,"# Nslices: \n");
	  fprintf(fout," %d\n",NSlices);
	}

      printf("The run used %d slices\n",NSlices);

      ntrees_tot=nbranch_tot=0;
      for (ThisSlice=0; ThisSlice<NSlices; ThisSlice++)
	{
	  u=fread(&dummy,sizeof(int),1,fd);
	  u=fread(&ntrees,sizeof(int),1,fd);
	  u=fread(&nbranch,sizeof(int),1,fd);
	  u=fread(&dummy,sizeof(int),1,fd);
	  ntrees_tot+=ntrees;
	  nbranch_tot+=nbranch;

	  if (print)
	    {
	      fprintf(fout,"# Ntrees & Nbranches: \n");
	      fprintf(fout," %d  %d\n",ntrees, nbranch);
	    }


	  for (itree=0; itree<ntrees; itree++)
	    {

	      u=fread(&dummy,sizeof(int),1,fd);
	      u=fread(&thistree,sizeof(int),1,fd);
	      u=fread(&nbranch_tree,sizeof(int),1,fd);
	      u=fread(&dummy,sizeof(int),1,fd);

	      if (print)
		fprintf(fout,"#Tree %d, Nbranches=%d\n", thistree, nbranch_tree);

	      for (ibranch=0; ibranch<nbranch_tree; ibranch++)
		{
		  u=fread(&dummy,sizeof(int),1,fd);
		  u=fread(&histdata,sizeof(hist),1,fd);
		  u=fread(&dummy,sizeof(int),1,fd);

		  if (First)
		    printf(" %12Ld %6d %6d %6d %9d %9d %9.4f %9.4f %9.4f\n",
			   histdata.name,
			   histdata.nick,
			   histdata.ll,
			   histdata.mw,
			   histdata.mass,
			   histdata.mam,
			   histdata.zme,
			   histdata.zpe,
			   histdata.zap);

		  if (print)
		    fprintf(fout," %12Ld %6d %6d %6d %9d %9d %9.4f %9.4f %9.4f\n",
			   histdata.name,
			   histdata.nick,
			   histdata.ll,
			   histdata.mw,
			   histdata.mass,
			   histdata.mam,
			   histdata.zme,
			   histdata.zpe,
			   histdata.zap);

		}
	      //if (First)
	      First=0;
	    }
	}

      printf("I found %Ld trees and %Ld branches in the file\n",ntrees_tot,nbranch_tot);
      printf("This is the last branch in the catalog:\n");
      printf(" %12Ld %6d %6d %6d %9d %9d %9.4f %9.4f %9.4f\n",
	     histdata.name,
	     histdata.nick,
	     histdata.ll,
	     histdata.mw,
	     histdata.mass,
	     histdata.mam,
	     histdata.zme,
	     histdata.zpe,
	     histdata.zap);

      if (print)
	fclose(fout);
    }
  else
    /* CATALOG AT FIXED REDSHIFT */
    {
      /* reads the catalog */
      sprintf(filename,"pinocchio.%s.%s.catalog.out",redshift,runname);
      printf("Opening file %s\n",filename);
      printf("I will write out the first ten halos\n");

      fd=fopen(filename,"r");
      if (fd==0x0)
	{
	  printf("Error: catalog file %s not found\n",filename);
	  return 1;
	}

      if (print)
	{
	  sprintf(outputf,"pinocchio.%s.%s.catalog.out.ascii",redshift,runname);
	  printf("I will write the ascii catalog in the file %s\n",outputf);
	  fout=fopen(outputf,"w");
	}

      n=0;
      u=fread(&ninhead,sizeof(int),1,fd);
      ninhead/=sizeof(int);
      u=fread(&nproc,sizeof(int),1,fd);
      if (ninhead==2)
	u=fread(&NSlices,sizeof(int),1,fd);
      u=fread(&dummy,sizeof(int),1,fd);
      printf("The file has been written by %d tasks\n",nproc);
      if (ninhead==2)
	printf("The box has been fragmented in %d slices\n",NSlices);
      else
	printf("This old version catalog does not give the number of slices used\n");

      if (print)
	{
	  fprintf(fout,"# Group catalog for redshift %s\n",redshift);
	  fprintf(fout,"#    1) group ID\n");
	  fprintf(fout,"#    2) group mass\n");
	  fprintf(fout,"# 3- 5) initial position\n");
	  fprintf(fout,"# 6- 8) final position\n");
	  fprintf(fout,"# 9-11) velocity\n");
	  fprintf(fout,"#   12) number of particles\n");
	  fprintf(fout,"#\n");
	}

      nblocks=0;
      while (!feof(fd))
	{
	  if (!fread(&dummy,sizeof(int),1,fd))
	    break;
	  ++nblocks;
	  u=fread(&ngood,sizeof(int),1,fd);
	  u=fread(&dummy,sizeof(int),1,fd);
	  n+=ngood;
	  printf("  found %d halos in block N. %d\n",ngood,nblocks);
	  for (igood=0; igood<ngood; igood++)
	    {
	      u=fread(&dummy,sizeof(int),1,fd);
	      u=fread(&catdata,dummy,1,fd);
	      u=fread(&dummy,sizeof(int),1,fd);

	      if (nblocks==1 && igood<10)
		printf(" %12Ld %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
		       catdata.id, catdata.M,
		       catdata.q[0], catdata.q[1], catdata.q[2],
		       catdata.x[0], catdata.x[1], catdata.x[2],
		       catdata.v[0], catdata.v[1], catdata.v[2],
		       catdata.N);

	      if (print)
		fprintf(fout," %12Lu %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
			catdata.id,
			catdata.M,
			catdata.q[0],catdata.q[1],catdata.q[2],
			catdata.x[0],catdata.x[1],catdata.x[2],
			catdata.v[0],catdata.v[1],catdata.v[2],
			catdata.N);
	    }
	}
      fclose(fd);
      if (print)
	fclose(fout);

      printf("I found %Ld halos in the file\n",n);
      printf("This is the last halo in the catalog:\n");
      printf(" %12Ld %13.6e %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %12d\n",
	     catdata.id, catdata.M,
	     catdata.q[0], catdata.q[1], catdata.q[2],
	     catdata.x[0], catdata.x[1], catdata.x[2],
	     catdata.v[0], catdata.v[1], catdata.v[2],
	     catdata.N);

      if (ninhead==1)
	{
	  printf("   nproc=%d,  nblocks=%d,  reconstructed nslices=%d\n",nproc,nblocks,nblocks/nproc);
	  if (nblocks%nproc)
	    printf("WARNING: the numbers of blocks and tasks are not consistent!\n");
	}

      else
	{
	  printf("   nproc=%d,  nblocks=%d,  read nslices=%d\n",nproc,nblocks,NSlices);
	  if (nblocks != nproc*NSlices)
	    printf("WARNING: the number of blocks is not as expected: nproc * nslices = %d\n",
		   nproc*NSlices);
	}
    }

  return 0;
    
}


