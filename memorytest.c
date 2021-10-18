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

//#define VERBOSE

void abort_code(void);
int my_initialization(void);

int main(int argc, char **argv)
{

  int i,j,k, i1,j1,k1, i2,j2,k2,NS2, NN,NN1,N1,N2,N3,ssafe;
  double tdis,size,this,BytesPerParticle,FmaxBPP,FragG_BPP,FragP_BPP,
    TotalNP,TotalNP_pertask,ratio,smallest,cc,MemPerTask;

#ifdef SCALE_DEPENDENT_GROWTH
  SDGM.flag=-1;
#endif

  int NSlices=0, Nplanes;

  /* Initialize MPI */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);

  if (NTasks>1)
    {
      printf("Please run this code on 1 task only\n");
      MPI_Finalize();
      return 0;
    }


  if (argc<3)
    {
      printf("Usage: memorytest paramfile NTasks\n");
      MPI_Finalize();
      return 0;
    }


  strcpy(params.ParameterFile,argv[1]);
  if (my_initialization())
    abort_code();

  NTasks = atoi(argv[2]);

  printf("Number of MPI tasks used for the run: %d\n",NTasks);

  Nplanes=MyGrids[0].GSglobal_x/NTasks;
  if (Nplanes*NTasks<MyGrids[0].GSglobal_x)
    ++Nplanes;
  printf("Nplanes=%d\n",Nplanes);


  /* typical displacement at zlast */
  tdis = GrowingMode(outputs.zlast) * sqrt( DisplVariance(params.InterPartDist) );

  /* mass of the largest halo expected in the box */
  params.Largest=1.e18;
  cc=1./pow(params.BoxSize_htrue,3.0);

  double aa=AnalyticMassFunction(params.Largest,outputs.zlast);
  while (aa*params.Largest<cc)
    {
      params.Largest*=0.99; 
      aa=AnalyticMassFunction(params.Largest,outputs.zlast);
    }

  /* boundary layer */
  size=SizeForMass(params.Largest);
  subbox.SafetyBorder = params.BoundaryLayerFactor * size;
  subbox.safe = (int)(subbox.SafetyBorder/params.InterPartDist)+1;

  if (!ThisTask)
    {
      printf("\n");
      printf("Determination of the boundary layer\n");
      printf("   growing mode at z=%f: %f\n",outputs.zlast, GrowingMode( outputs.zlast ));
      printf("   largest halo expected in this box at z=%f: %e Msun\n",
	     outputs.zlast, params.Largest);
      printf("   its Lagrangian size: %f Mpc\n",size);
      printf("   typical displacement: %f \n",tdis);
      printf("   the boundary layer will be %f, a factor of %f with respect to the typical displacement\n",
	     subbox.SafetyBorder, subbox.SafetyBorder/tdis);
    }

  /* finds the optimal number of sub-boxes to use for the fragmentation */
  ssafe=2.*subbox.safe;
  FmaxBPP = (double)sizeof(product_data) + 10.0*(double)sizeof(double) + 
    (double)sizeof(int) * (double)NTasks / (double)MyGrids[0].GSglobal_z;
  TotalNP = (double)MyGrids[0].GSglobal_x * (double)MyGrids[0].GSglobal_y * (double)MyGrids[0].GSglobal_z;
  TotalNP_pertask = TotalNP/(double)NTasks;

  FragP_BPP=(double)sizeof(product_data);
  FragG_BPP=3.0*(double)sizeof(int)+(double)(sizeof(group_data) +sizeof(histories_data))/10.0;
#ifdef PLC
  FragG_BPP+=(double)sizeof(plcgroup_data)/10.;
#endif

  printf("Number of particles per task: %f millions\n", TotalNP_pertask/1.e6);
  printf("Total number of particles: %f millions\n", TotalNP/1.e6);

  smallest=1.e10;
  NSlices=0;

  do
    {
      ++NSlices;

      BytesPerParticle=1.e10;
      printf("Iteration with NSlices=%d:\n",NSlices);
      for (k=1; k<=NTasks; k++)
	for (j=1; j<=NTasks/k; j++)
	  for (i=1; i<=NTasks/k/j; i++)
	    /* the three indices must be exact divisors of the three grid lengths */
	    if (i*j*k==NTasks)
	      {
		/* number of particles in the sub-box */
		N1 = find_length(MyGrids[0].GSglobal_x,i,0);
		N2 = find_length(MyGrids[0].GSglobal_y,j,0);
		N3 = find_length(MyGrids[0].GSglobal_z,k*NSlices,0);
		if (N1<ssafe || N2<ssafe || N3<ssafe)
		  continue;
		NN = (N1 + (i==1? 0 : ssafe))
		  *  (N2 + (j==1? 0 : ssafe)) 
		  *  (N3 + (k*NSlices==1? 0 : ssafe));

		ratio = (double)NN/TotalNP_pertask;
		if (NSlices>1)
		  this=(double)sizeof(product_data) + ratio * (FragP_BPP + FragG_BPP);
		else		    
		  this=( (double)sizeof(product_data) > ratio * FragG_BPP ?
			 (double)sizeof(product_data) : ratio * FragG_BPP) +
		    ratio * FragP_BPP;
		if (this<FmaxBPP)
		  this=FmaxBPP;
		
		if (this<smallest)
		  {
		    smallest=this;
		    i2=i; j2=j; k2=k; NS2=NSlices;
		  }

#ifdef VERBOSE
		printf ("   case %d %d %d: N1=%d, N2=%d, N3=%d, NN=(%d * %d * %d)=%d - %f  (overhead=%f)\n",
			i,j,k,N1,N2,N3,
			(N1 + (i==1? 0 : ssafe)),
			(N2 + (j==1? 0 : ssafe)),
			(N3 + (k==1? 0 : ssafe)),NN,this,(double)NN/(double)(N1*N2*N3));
#endif

		if (this < BytesPerParticle)
		  {
		    BytesPerParticle=this;
		    NN1=NN;
		    i1=i;
		    j1=j;
		    k1=k;
		  }
	      }

      printf("subboxes: %d %d %d\n",i1,j1,k1);
      printf("bytes per particle: %f\n",BytesPerParticle);

      if (BytesPerParticle>1000.)
	break;
    }
  while (BytesPerParticle>params.MaxMemPerParticle);


  if (BytesPerParticle>1000.)
    {
      printf("ERROR: no possible division of sub-boxes found up to Nslices=%d\n", 
	     NSlices);
      printf("lowest possible value of memory per particle is %f ",smallest);
      printf("found on a subdivision %d-%d-%d on %d slices\n",i2,j2,k2,NS2);
      printf("please decrease BoundaryLayerFactor or increase MaxMemPerParticle\n");

      MPI_Finalize();
      return NSlices;

    }

  subbox.nbox_x=i1;
  subbox.nbox_y=j1;
  subbox.nbox_z_thisslice=k1;
  subbox.nbox_z_allslices=k1*NSlices;

  subbox.safe_x = (subbox.nbox_x>1 ? subbox.safe : 0);
  subbox.safe_y = (subbox.nbox_y>1 ? subbox.safe : 0);
  subbox.safe_z = (subbox.nbox_z_allslices>1 ? subbox.safe : 0);

  subbox.pbc_x = (subbox.nbox_x==1);
  subbox.pbc_y = (subbox.nbox_y==1);
  subbox.pbc_z = (subbox.nbox_z_allslices==1);

  /* this will be mybox for the first slice */
  NN=subbox.nbox_y*subbox.nbox_z_thisslice;
  subbox.mybox_x=ThisTask/NN;
  NN1=ThisTask-subbox.mybox_x*NN;
  subbox.mybox_y=NN1/subbox.nbox_z_thisslice;
  subbox.mybox_z=NN1-subbox.mybox_y*subbox.nbox_z_thisslice;

  subbox.Lgrid_x = find_length(MyGrids[0].GSglobal_x,subbox.nbox_x,subbox.mybox_x);
  subbox.Lgrid_y = find_length(MyGrids[0].GSglobal_y,subbox.nbox_y,subbox.mybox_y);
  subbox.Lgrid_z = find_length(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.Lgwbl_x = subbox.Lgrid_x + 2*subbox.safe_x; 
  subbox.Lgwbl_y = subbox.Lgrid_y + 2*subbox.safe_y;
  subbox.Lgwbl_z = subbox.Lgrid_z + 2*subbox.safe_z;

  subbox.Npart = subbox.Lgwbl_x * subbox.Lgwbl_y * subbox.Lgwbl_z;

  subbox.start_x = find_start(MyGrids[0].GSglobal_x,subbox.nbox_x,subbox.mybox_x);
  subbox.start_y = find_start(MyGrids[0].GSglobal_y,subbox.nbox_y,subbox.mybox_y);
  subbox.start_z = find_start(MyGrids[0].GSglobal_z,subbox.nbox_z_allslices,subbox.mybox_z);

  subbox.stabl_x = subbox.start_x - subbox.safe_x;
  subbox.stabl_y = subbox.start_y - subbox.safe_y;
  subbox.stabl_z = subbox.start_z - subbox.safe_z;

  subbox.overhead=(double)subbox.Npart/(double)(subbox.Lgrid_x * subbox.Lgrid_y * subbox.Lgrid_z);

  double MemPerPlane = BytesPerParticle * (double)MyGrids[0].GSglobal_x * (double)MyGrids[0].GSglobal_x / 1024. / 1024. / 1024.;
  MemPerTask  = BytesPerParticle * TotalNP_pertask / 1024. / 1024. / 1024.;
  int NPlanesPerTask = (int)(params.MaxMem/1024.0/MemPerPlane);

  printf("\n");
  printf("FRAGMENTATION:\n");
  if (NSlices>1)
    printf("The box will be fragmented in %d slices\n",NSlices);
  else
    printf("The box will be fragmented in one go\n");
  printf("Number of sub-boxes per dimension: %d %d %d\n",subbox.nbox_x,subbox.nbox_y,subbox.nbox_z_allslices);
  printf("Boundary layer (true Mpc):         %f\n",subbox.SafetyBorder);
  printf("Boundary layer (gridpoints):       %d\n",subbox.safe);
  printf("Cores will work on a grid:         %d %d %d\n",subbox.Lgwbl_x,subbox.Lgwbl_y,subbox.Lgwbl_z);
  printf("Number of particles per core:      %d\n",subbox.Npart);
  printf("The resolved box will be:          %d %d %d\n",subbox.Lgrid_x,subbox.Lgrid_y,subbox.Lgrid_z);
  printf("Periodic boundary conditions:      %d %d %d\n",subbox.pbc_x,subbox.pbc_y,subbox.pbc_z);
  printf("Required bytes per fft particle:   %f\n",BytesPerParticle);
  printf("The overhead for fragmentation is: %f\n",subbox.overhead);
  printf("Required memory per task:          %4.0fMb - Maxmem=%dMb\n", MemPerTask*1024.,params.MaxMem);
  printf("\n");
  printf("\n");

  if (MemPerTask > params.MaxMem/1024.0)
    printf("FATAL ERROR: your requirements overshoot the available memory per MPI task, the run will fail\n");
  printf("Required memory per plane: %fGb\n", MemPerPlane);
  printf("You will be able to load at most %d planes per MPI task (%d in total over %d required)\n",
	 NPlanesPerTask,NPlanesPerTask * NTasks,MyGrids[0].GSglobal_x);
  if (NPlanesPerTask * NTasks < MyGrids[0].GSglobal_x)
    printf("FATAL ERROR: you will not be able to load all the planes on a task, the run will fail (%d)\n",NPlanesPerTask * NTasks);
  printf("Total required memory: %fGb\n",
	 BytesPerParticle * TotalNP / 1024. / 1024. / 1024.);


  MPI_Finalize();
  return NSlices;
}

void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}



int my_initialization(void)
{
  workspace = gsl_integration_workspace_alloc(NWINT);
  if (set_parameters())
    return 1;
#ifdef SCALE_DEPENDENT_GROWTH
  /* reads the P(k) tables from CAMB */
  if (read_power_table_from_CAMB())
    return 1;
#endif
  if (initialize_cosmology())
    return 1;
  if (set_grids())
    return 1;
  /* now it re-initializes the variance with a top-hat filter */
  WindowFunctionType=2;
  if (initialize_MassVariance())
    return 1;
#ifdef SCALE_DEPENDENT_GROWTH
  /* initialize scale-dependent growing modes and fomegas */
  if (initialize_ScaleDependentGrowth())
    return 1;
#endif
  return 0;
}
