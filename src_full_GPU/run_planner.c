/* ######HEADER###### */

#include "pinocchio.h"
#include "def_splines.h"

#define MARGIN 1.0
#define MAXCOUNT 10

void abort_code(void);

/* 
   Order-of-magnitude estimate of the needed overhead: 
   sigma < 6 : overhead = (sigma/6)**2   * 3.8 * 1.1 (margin)
   sigma > 6 : overhead = (sigma/6)**0.6 * 3.8 * 1.1 (margin)   
*/

int main(int argc, char **argv, char **envp)
{

  /* Initialize MPI */
  int got_level;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &got_level);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTasks);

#ifdef _OPENMP
  /* initialization of OpemMP */
#pragma omp parallel
  {
#pragma omp master
    internal.nthreads_omp = omp_get_num_threads();
  }
#endif


  if (NTasks>1)
    {
      printf("Please run this code in a scalar mode, on 1 task only\n");
      MPI_Finalize();
      return 0;
    }

  if (argc<4)
    {
      printf("Usage: run_planner paramfile RamPerNode (Gb) TasksPerNode [Nnodes]\n");
      MPI_Finalize();
      return 0;
    }

  int RamPerNode = atoi(argv[2]);
  int TasksPerNode = atoi(argv[3]);
  int RamPerTask = (int)( (1024. * (double)RamPerNode) / (double)TasksPerNode);
  double overhead;
  int ForceNnodes;
  if (argc>=5)
    {
      ForceNnodes = atoi(argv[4]);
    }
  else
    ForceNnodes = 0;

  printf("run_planner planning a run on nodes with %d Gb or RAM each, with %d tasks per node\n",RamPerNode,TasksPerNode);
  if (ForceNnodes)
    printf("The number of nodes will be forced to %d\n",ForceNnodes);


  printf("\n");
  printf("********************************************** \n");
  printf("     This is the standard PINOCCHIO output     \n");
  printf("********************************************** \n");

  greetings();

  workspace = gsl_integration_workspace_alloc(NWINT);

  memset(&params, 0, sizeof(param_data));
  strcpy(params.ParameterFile,argv[1]);

  if (set_parameters())
    exit(1);

  if (initialize_cosmology())
    exit(1);

  /* now it re-initializes the variance with a SKS filter */
  WindowFunctionType=1;
  if (initialize_MassVariance())
    return 1;

  printf("\n");
  printf("********************************************** \n");
  printf("     End of standard PINOCCHIO output \n");
  printf("********************************************** \n");
  printf("\n");

  /* sets the main grid variables */
  Ngrids=1;

  MyGrids=(grid_data*)malloc(Ngrids * sizeof(grid_data));

  MyGrids[0].GSglobal[_x_] = params.GridSize[_x_];
  MyGrids[0].GSglobal[_y_] = params.GridSize[_y_];
  MyGrids[0].GSglobal[_z_] = params.GridSize[_z_];

  MyGrids[0].Ntotal = (unsigned long long)MyGrids[0].GSglobal[_x_] * 
    (unsigned long long)MyGrids[0].GSglobal[_y_] * 
    (unsigned long long)MyGrids[0].GSglobal[_z_];

  MyGrids[0].BoxSize = params.BoxSize_htrue;

  memset(&memory,0,sizeof(memory_data));

  double sigma = sqrt(MassVariance(params.BoxSize_htrue/(double)params.GridSize[_x_]/PI/NYQUIST));
#ifndef CLASSIC_FRAGMENTATION
  /* needed overhead is estimated based on the mass variance on the grid */
  overhead = (sigma < 6 ? pow(sigma/6.0,2.0) * 3.8 * MARGIN : pow(sigma/6.0,0.6) * 3.8 * MARGIN);
#endif


#ifdef TWO_LPT
#ifdef THREE_LPT
  double baseline = 155;
#else
  double baseline = 130;
#endif
#else
  double baseline = 95;
#endif
#ifdef SNAPSHOT
  baseline += 4;
#endif


  /* setting MaxMemPerParticle */
  params.MaxMemPerParticle = baseline;
  double old=0.0; int Nnodes=0;

#ifndef CLASSIC_FRAGMENTATION
  printf("needed overhead: %f\n",overhead);
#endif
  int count=0, result;

  if (ForceNnodes)
    {
      /* in this case the number of nodes is fixed */
      params.MaxMemPerParticle = (int)(ForceNnodes / ((double)MyGrids[0].Ntotal / GBYTE / (double)RamPerNode) ) - 1;
      if (params.MaxMemPerParticle<baseline)
	{
	  printf("**************************************************************************\n");
	  printf("The number of nodes you provided is insufficient, I will compute it myself\n");
	  printf("**************************************************************************\n");
	  Nnodes=0;
	}
      else
	{
	  NTasks = ForceNnodes * TasksPerNode;

	  MyGrids[0].ParticlesPerTask = (int)((double)MyGrids[0].Ntotal / (double)NTasks);
	  /* proxies for pfft allocation */
	  MyGrids[0].total_local_size = MyGrids[0].ParticlesPerTask;
	  MyGrids[0].total_local_size_fft = MyGrids[0].total_local_size
	    * (MyGrids[0].GSglobal[_x_] + 1) / MyGrids[0].GSglobal[_x_];

	  printf("\n");
	  printf("********************************************** \n");
	  printf("                set_subboxes output            \n");
	  printf("********************************************** \n");
	  printf("\n");
	  result = set_subboxes();
	  printf("\n");
	  printf("********************************************** \n");
	  printf("\n");

	  if (result)
	    {
	      printf("**************************************************************************\n");
	      printf("The number of nodes you provided is insufficient, I will compute it myself\n");
	      printf("**************************************************************************\n");
	      Nnodes=0;
	    }
	  else
	    {
	      printf("I found a successful configuration on %d nodes\n",ForceNnodes);
	      Nnodes = ForceNnodes;
	    }
	}
    }

  if (!Nnodes)
    {
      /* in this case it finds the number of needed nodes */
      do
	{
#ifdef CLASSIC_FRAGMENTATION
	  if (old==params.MaxMemPerParticle)
	    params.MaxMemPerParticle+=10;
#endif

	  old = params.MaxMemPerParticle;
	  /* estimate of the number of needed nodes */
	  Nnodes = (int)((double)MyGrids[0].Ntotal * params.MaxMemPerParticle / GBYTE / (double)RamPerNode + 1);
	  NTasks = Nnodes * TasksPerNode;

	  MyGrids[0].ParticlesPerTask = (int)((double)MyGrids[0].Ntotal / (double)NTasks);

	  /* proxies for pfft allocation */
	  MyGrids[0].total_local_size = MyGrids[0].ParticlesPerTask;
	  MyGrids[0].total_local_size_fft = MyGrids[0].total_local_size
	    * (MyGrids[0].GSglobal[_x_] + 1) / MyGrids[0].GSglobal[_x_];

	  printf("\n");
	  printf("********************************************** \n");
	  printf("  set_subboxes output, IGNORE ERROR MESSAGES   \n");
	  printf("********************************************** \n");
	  printf("\n");
	  result = set_subboxes();
	  printf("\n");
	  printf("********************************************** \n");
	  printf("\n");

#ifdef CLASSIC_FRAGMENTATION
	  /* I add a 5% margin here */
	  overhead = (double)subbox.Npart / (double)MyGrids[0].ParticlesPerTask * 1.05;
#endif

	  params.MaxMemPerParticle = 
	    overhead * (double)(sizeof(product_data) + FRAGFIELDS * sizeof(int)) +
	    (double)(memory.prods + memory.groups) / (double)MyGrids[0].ParticlesPerTask;
	  if (params.MaxMemPerParticle<baseline)
	    params.MaxMemPerParticle = baseline;

	  printf("I tried MaxMemPerParticle = %d; Nnodes = %d; NTasks = %d; MyGrids[0].ParticlesPerTask = %d; new MaxMemPerParticle = %d\n",(int)old,Nnodes,NTasks,MyGrids[0].ParticlesPerTask,(int)params.MaxMemPerParticle);
	  count++;

	} while ((params.MaxMemPerParticle > old || result==1) && count<=MAXCOUNT);

      if (count>MAXCOUNT)
	{
	  printf("Sorry, I could not find an acceptable configuration for this run!\n");
#ifdef CLASSIC_FRAGMENTATION
	  printf("This is too much for CLASSIC_FRAGMENTATION, please uncomment this option in the Makefile\n");
#else
	  printf("Try to contact pierluigi.monaco@inaf for some help\n");
#endif
	  return 1;
	}
    }


  params.MaxMem            = RamPerTask;
  params.PredPeakFactor    = (sigma < 6 ? (sigma-6)/3.0 + 1.5 : 1.5);
  if (params.PredPeakFactor<0.6)
    params.PredPeakFactor = 0.6;
#ifdef CLASSIC_FRAGMENTATION
  params.BoundaryLayerFactor = 1;
#else
  params.BoundaryLayerFactor = 3;
#endif

  printf("\n********************************************** \n");
  printf("Density std dev on grid for this run: %f\n",sigma);
  printf("We assume a value of MaxMemPerParticle of %d\n",(int)params.MaxMemPerParticle);
  printf("The number of nodes needed for this run is: %d\n",Nnodes);
  printf("Number of MPI tasks used for the run: %d\n",NTasks);

  if ((int)(MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle / MBYTE + 1.0) > params.MaxMem)
    {
      printf("ERROR: MaxMem of %d Mb per task is insufficient to store %d bytes per particle\n",
	     params.MaxMem, (int)params.MaxMemPerParticle);
      printf("       please increase MaxMem to at least %d\n",
	     (int)(MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle / MBYTE + 1.0));

      return 1;
    }


  /* now it re-initializes the variance with a top-hat filter */
  WindowFunctionType=2;
  if (initialize_MassVariance())
    return 1;

  printf("\n");
  printf("A successful PINOCCHIO run is determined by these parameters:\n");
  printf("- MaxMem: is the memory available to the MPI task;\n");
  printf("- MaxMemPerParticle: is the memory per particle that can be allocated;\n");
  printf("- BoundaryLayerFactor: is the depth of the boundary layer.\n");
  printf("\n");
  printf("  Each task adds to the fragmentation sub-box a layer as deep\n");
  printf("  as BoundaryLayerFactor times the Lagrangian size of the largest\n");
  printf("  halo expected in the box. The augmented sub-box will be smaller than the whole box.\n");
  printf("  In the new fragmentation, halos at the border of the fragmentation sub-volume\n");
  printf("  are augmented by BoundaryLayerFactor times their Lagrangian size.\n");
  printf("  Suggested value: 1.0 for the classic fragmentation, 3.0 for the new fragmentation.\n");
  printf("\n");
  printf("- PredPeakFactor: the number of allocated peaks will be PredPeakFactor times 1/6 of\n");
  printf("  the particles in the sub-box (without boundary layer). If this is kept small \n");
  printf("  there might be not enough space to allocate peaks.\n");
  printf("  Suggested value: 0.6 at low resolution (Mpart~1e11 Msun/h), \n");
  printf("  0.8 at medium resolution (Mpart~1e9 Msun/h), >1.0 at higher resolution.\n");
  printf("  At the end of a run, the code suggests a minimum value for PredPeakFactor,\n");
  printf("  but keep a margin to it.\n");
  printf("\n");

  struct map{
    int lz;
    double oh,am,m1,m2,m3,m4,m5,m6,m7,m8;
  } mymap;

  mymap.am=(double)memory.all/MBYTE;
  mymap.oh=(double)subbox.Nalloc/(double)MyGrids[0].ParticlesPerTask;
  mymap.m1=(double)memory.prods/(double)MyGrids[0].ParticlesPerTask;
  mymap.m2=(double)(memory.fields+memory.fields_to_keep)/(double)MyGrids[0].ParticlesPerTask;
  mymap.m3=(double)memory.fft/(double)MyGrids[0].ParticlesPerTask;
  mymap.m4=(double)memory.fmax_total/(double)MyGrids[0].ParticlesPerTask;
  mymap.m5=(double)memory.frag_prods/(double)MyGrids[0].ParticlesPerTask;
  mymap.m6=(double)(memory.groups+memory.frag_arrays)/(double)MyGrids[0].ParticlesPerTask;
  mymap.m7=(double)memory.frag_total/(double)MyGrids[0].ParticlesPerTask;
  mymap.m8=(double)memory.all/(double)MyGrids[0].ParticlesPerTask;


  printf("\n");
  printf("Map of memory usage for Task 0:\n");
  printf("Task N.    mem(MB) overhead   products   fields     ffts     fmax  frag pr.  groups fragment  total bytes per particle\n");

  printf("%6d   %8.0f  %6.1f       %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f   %6.1f\n",
	 0, mymap.am, mymap.oh, mymap.m1, mymap.m2, mymap.m3, mymap.m4, mymap.m5, mymap.m6, mymap.m7, mymap.m8);

  printf("\n");
  printf("Complete memory map\n");
  printf("  memory.prods:           %12zu, %6.1f bpp\n",memory.prods,           (double)memory.prods/(double)MyGrids[0].total_local_size);
  printf("  memory.fields_to_keep   %12zu, %6.1f bpp\n",memory.fields_to_keep,  (double)memory.fields_to_keep/(double)MyGrids[0].total_local_size);
  printf("  memory.fields           %12zu, %6.1f bpp\n",memory.fields,          (double)memory.fields/(double)MyGrids[0].total_local_size);
  printf("  memory.first_allocated: %12zu, %6.1f bpp\n",memory.first_allocated, (double)memory.first_allocated/(double)MyGrids[0].total_local_size);
  printf("  memory.fft:             %12zu, %6.1f bpp\n",memory.fft,             (double)memory.fft/(double)MyGrids[0].total_local_size);
  printf("  memory.fmax_total:      %12zu, %6.1f bpp\n",memory.fmax_total,      (double)memory.fmax_total/(double)MyGrids[0].total_local_size);
  printf("  memory.frag_prods:      %12zu, %6.1f bpp\n",memory.frag_prods,      (double)memory.frag_prods/(double)MyGrids[0].total_local_size);
  printf("  memory.frag_arrays:     %12zu, %6.1f bpp\n",memory.frag_arrays,     (double)memory.frag_arrays/(double)MyGrids[0].total_local_size);
  printf("  memory.groups:          %12zu, %6.1f bpp\n",memory.groups,	        (double)memory.groups/(double)MyGrids[0].total_local_size);
  printf("  memory.frag_allocated:  %12zu, %6.1f bpp\n",memory.frag_allocated,  (double)memory.frag_allocated/(double)MyGrids[0].total_local_size);	       
  printf("  memory.frag_total:      %12zu, %6.1f bpp\n",memory.frag_total,      (double)memory.frag_total/(double)MyGrids[0].total_local_size);	       
  printf("  memory.all:             %12zu, %6.1f bpp\n",memory.all,	        (double)memory.all/(double)MyGrids[0].total_local_size);

  printf("\n");
  printf("Number of nodes needed for this run:    %d\n",Nnodes);
  if (ForceNnodes && ForceNnodes!=Nnodes)
    printf("WARNING: this is different from your request!\n");
  printf("Number of MPI tasks used for the run:   %d\n",NTasks);
  printf("Each MPI task will need at least %f Mb\n",MyGrids[0].ParticlesPerTask * params.MaxMemPerParticle / MBYTE + 1.0);
  double Nfloat = (double)MyGrids[0].Ntotal*(double)params.MaxMemPerParticle / (RamPerNode * GBYTE);
  printf("This run will occupy memory of %5.2f nodes, %4.2f percent of available memory\n",
	 Nfloat, 100.*Nfloat / (double)Nnodes);
  printf("Density standard deviation on the grid: %f\n",sigma);
  printf("Predicted overhead:                     %f\n",overhead);
  printf("You can copy and paste these into the parameter file:\n");
  printf("   MaxMem                %d\n",params.MaxMem);
  printf("   MaxMemPerParticle     %d\n",(int)params.MaxMemPerParticle);
  printf("   PredPeakFactor        %3.1f\n",params.PredPeakFactor);
  printf("   BoundaryLayerFactor   %3.1f\n",params.BoundaryLayerFactor);
  if (sigma>5.)
    printf("WARNING: this is a %shigh resolution run, it may fail at first attempt,\n  but in this case the code will suggest where the problem is.\n",(sigma>6?"very ":""));

#ifndef CLASSIC_FRAGMENTATION
  if (sigma>6.)
    printf("  Please consider using 2.5 or even 2.0 for BoundaryLayerFactor to increase the probability of success\n");
#endif
  printf("\n");

  (void)estimate_file_size();
  printf("\n");

  printf("**************************************************************************\n");
  printf(" I'm now checking parameters and directives, this may give error messages \n");
  printf("**************************************************************************\n");
  if (check_parameters_and_directives())
    {
      printf("**************************************************************************\n");
      printf("Please follow these instructions before running the code\n\n");
    }
  else
    printf("**************************************************************************\n");


  printf("This is the number of sub-boxes per dimension: %d %d %d\n",subbox.nbox[_x_],subbox.nbox[_y_],subbox.nbox[_z_]);
  printf("Their products MUST give the number of tasks.\n");
  printf("The more similar they are, the better (unless some of them is 1).\n");
  printf("If the run is large and the three numbers are not as similar as possible,\ntry to change the number of tasks per node or to ask for a specific number of nodes to achieve a better balance.\n\n");

  /* done */
  printf("run_planner done!\n");
  MPI_Finalize();

  return 0;
}


void abort_code(void)
{
  printf("Task %d aborting...\n",ThisTask);
  MPI_Abort(MPI_COMM_WORLD,1);
}
