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
   The similarity of this code with the read_parameter_file() function in the Gadget code,
   by V. Springel (http://www.gadgetcode.org), is not a coincidence.  
   I have adopted the same code structure.
*/


#include "pinocchio.h"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define LOGICAL 4
#define INT3 5
#define DOUBLE3 6
#define INT_SKIP 99
#define MAXTAGS 100

int read_parameter_file()
{
  FILE *fd;
  char buf[SBLENGTH], buf1[SBLENGTH], buf2[SBLENGTH], buf3[SBLENGTH], buf4[SBLENGTH];
  int i, j, nt, number_of_fields;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][SBLENGTH];
  double z;

  if (!ThisTask)
    {

      nt=0;

      /* list of requested parameters */
      strcpy(tag[nt], "RunFlag");
      addr[nt] = params.RunFlag;
      id[nt++] = STRING;

      strcpy(tag[nt], "OutputList");
      addr[nt] = params.OutputList;
      id[nt++] = STRING;

      strcpy(tag[nt], "Omega0");
      addr[nt] = &params.Omega0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaLambda");
      addr[nt] = &params.OmegaLambda;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PrimordialIndex");
      addr[nt] = &params.PrimordialIndex;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "Sigma8");
      addr[nt] = &params.Sigma8;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "Hubble100");
      addr[nt] = &params.Hubble100;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "OmegaBaryon");
      addr[nt] = &params.OmegaBaryon;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DEw0");
      addr[nt] = &params.DEw0;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "DEwa");
      addr[nt] = &params.DEwa;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "StartingzForPLC");
      addr[nt] = &params.StartingzForPLC;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "LastzForPLC");
      addr[nt] = &params.LastzForPLC;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PLCAperture");
      addr[nt] = &params.PLCAperture;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "PLCProvideConeData");
      addr[nt] = &params.PLCProvideConeData;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "PLCCenter");
      addr[nt] = &(params.PLCCenter);
      id[nt++] = DOUBLE3;

      strcpy(tag[nt], "PLCAxis");
      addr[nt] = &(params.PLCAxis);
      id[nt++] = DOUBLE3;

      strcpy(tag[nt], "BoxSize");
      addr[nt] = &params.BoxSize;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "BoxInH100");
      addr[nt] = &params.BoxInH100;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "RandomSeed");
      addr[nt] = &params.RandomSeed;
      id[nt++] = INT;

      strcpy(tag[nt], "GridSize");
      addr[nt] = &(params.GridSize[0]);
      id[nt++] = INT;

      strcpy(tag[nt], "BoundaryLayerFactor");
      addr[nt] = &params.BoundaryLayerFactor;
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "MinHaloMass");
      addr[nt] = &params.MinHaloMass;
      id[nt++] = INT;

      strcpy(tag[nt], "WriteFmax");
      addr[nt] = &params.WriteFmax;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "WriteVmax");
      addr[nt] = &params.WriteVmax;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "WriteRmax");
      addr[nt] = &params.WriteRmax;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "NumFiles");
      addr[nt] = &params.NumFiles;
      id[nt++] = INT_SKIP;

      strcpy(tag[nt], "MaxMem");
      addr[nt] = &params.MaxMem;
      id[nt++] = INT;

      strcpy(tag[nt], "AnalyticMassFunction");
      addr[nt] = &params.AnalyticMassFunction;
      id[nt++] = INT;

      strcpy(tag[nt], "CatalogInAscii");
      addr[nt] = &params.CatalogInAscii;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "DoNotWriteCatalogs");
      addr[nt] = &params.DoNotWriteCatalogs;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "DoNotWriteHistories");
      addr[nt] = &params.DoNotWriteHistories;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "WriteSnapshot");
      addr[nt] = &params.WriteSnapshot;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "WriteTimelessSnapshot");
      addr[nt] = &params.WriteTimelessSnapshot;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "OutputInH100");
      addr[nt] = &params.OutputInH100;
      id[nt++] = LOGICAL;

      strcpy(tag[nt], "InputSpectrum_UnitLength_in_cm");
      addr[nt] = &(params.InputSpectrum_UnitLength_in_cm);
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "FileWithInputSpectrum");
      addr[nt] = params.FileWithInputSpectrum;
      id[nt++] = STRING;

      strcpy(tag[nt], "WDM_PartMass_in_kev");
      addr[nt] = &(params.WDM_PartMass_in_kev);
      id[nt++] = DOUBLE;

      strcpy(tag[nt], "TabulatedEoSfile");
      addr[nt] = params.TabulatedEoSfile;
      id[nt++] = STRING;

#ifdef READ_PK_TABLE
      strcpy(tag[nt], "CAMBMatterFileTag");
      addr[nt] = params.camb.MatterFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CAMBTransferFileTag");
      addr[nt] = params.camb.TransferFile;
      id[nt++] = STRING;

      strcpy(tag[nt], "CAMBRunName");
      addr[nt] = params.camb.RunName;
      id[nt++] = STRING;

      strcpy(tag[nt], "CAMBRedsfhitsFile");
      addr[nt] = params.camb.RedshiftsFile;
      id[nt++] = STRING;
#endif

      strcpy(tag[nt], "MaxMemPerParticle");
      addr[nt] = &(params.MaxMemPerParticle);
      id[nt++] = DOUBLE;


      for (j=0; j<nt; j++)     /* All logical tags are FALSE by default */
	if (id[j]==LOGICAL)
	  *((int *) addr[j]) = 0;

      printf("Reading parameters from file %s\n",params.ParameterFile);
      fflush(stdout);

      if((fd = fopen(params.ParameterFile, "r")))
	{
	  while(!feof(fd))
	    {
	      *buf = 0;
	      (void)fgets(buf, SBLENGTH, fd);
	      number_of_fields = sscanf(buf, "%s %s %s %s", buf1, buf2, buf3, buf4);
	      if(number_of_fields < 1)
		continue;

	      if(buf1[0] == '#' || buf1[0]=='%')
		continue;

	      /* searches for a match with the tag list */
	      for(i = 0, j = -1; i < nt; i++)
		if(strcmp(buf1, tag[i]) == 0)
		  {
		    j = i;
		    tag[i][0] = 0;
		    break;
		  }

	      if(j >= 0)
		{
		  switch (id[j])
		    {
		    case DOUBLE:
		      if (number_of_fields<2)
			j=-10;
		      else
			*((double *) addr[j]) = atof(buf2);
		      break;

		    case STRING:
		      if (number_of_fields<2)
			j=-10;
		      else
			strcpy(addr[j], buf2);
		      break;

		    case INT:
		      if (number_of_fields<2)
			j=-10;
		      else
			*((int *) addr[j]) = atoi(buf2);
		      break;

		    case INT_SKIP:
		      if (number_of_fields<2)
			j=-10;
		      else
			*((int *) addr[j]) = atoi(buf2);
		      break;

		    case INT3:
		      if (number_of_fields<4)
			j=-10;
		      else
			{
			  *((int *)addr[j]  ) = atoi(buf2);
			  *((int *)addr[j]+1) = atoi(buf3);
			  *((int *)addr[j]+2) = atoi(buf4);
			}
		      break;

		    case DOUBLE3:
		      if (number_of_fields<4)
			j=-10;
		      else
		        {
			  *((double *) addr[j])   = atof(buf2);
			  *((double *) addr[j]+1) = atof(buf3);
			  *((double *) addr[j]+2) = atof(buf4);
		        }
		      break;

		    case LOGICAL:
		      *((int *) addr[j]) = 1;
		      break;

		    default:
		      break;

		    }

		  if (j<0)
		    {
		      printf("ERROR on task 0: file %s, Tag %s has no argument.\n",
			     params.ParameterFile, buf1);
		      fflush(stdout);
		      return 1;
		    }
		}
	      else
		{
		  printf("ERROR on task 0: file %s, Tag %s is unknown.\n",
			 params.ParameterFile, buf1);
		  fflush(stdout);
		  return 1;
		}
	    }
	  fclose(fd);

	}
      else
	{
	  printf("ERROR: Parameter file %s not found.\n", params.ParameterFile);
	  return 2;
	}

      /* Checks that all parameters are found */
      for(i = 0; i < nt; i++)
	{
	  if(*tag[i])
	    {
	      if (id[i]==LOGICAL || id[i]==INT_SKIP || id[i]==DOUBLE3) 
		*((int *) addr[i]) = 0;
	      else
		{
		  printf("ERROR on task 0: I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], params.ParameterFile);
		  fflush(stdout);
		  return 1;
		}

	    }
	}

      /* Now opens and reads the list of wanted outputs */
      outputs.n=0;
      if((fd = fopen(params.OutputList, "r")))
	while(!feof(fd))
	  {
	    *buf = 0;
	    (void)fgets(buf, SBLENGTH, fd);
	    if (buf[0] != '#' && buf[0] != '%' && sscanf(buf,"%lf",&z) != EOF)
	      outputs.z[outputs.n++]=z;
	    if (outputs.n == MAXOUTPUTS)
	      {
		printf("ERROR on task 0: too many output redshifts in file %s\n",params.OutputList);
		printf("                 increase MAXOUTPUTS in pinocchio.h\n");
		fflush(stdout);
		return 1;
	      }

	  }
      else
	{
	  printf("Error on task 0: file %s not found\n",params.OutputList);
	  fflush(stdout);
	  return 1;
	}
      fclose(fd);
      outputs.zlast=outputs.z[outputs.n-1];

      if (outputs.n<=0)
	{
	  printf("ERROR on task 0: no valid lines found in file %s\n",params.OutputList);
	  fflush(stdout);
	  return 1;
	}

      params.GridSize[1]=params.GridSize[2]=params.GridSize[0];

    }

  /* processor 0 broadcasts the parameters to all other processors */
  MPI_Bcast(&params, sizeof(param_data), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&outputs.n, sizeof(output_data), MPI_BYTE, 0, MPI_COMM_WORLD);

  return 0;
}
