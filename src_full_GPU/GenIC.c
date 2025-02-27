/*****************************************************************
 *                        PINOCCHI0  V4.0                        *
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

  The code contained in this file has been adapted from the 2LPTic code
  by Roman Scoccimarro, downloadable from:

  http://cosmo.nyu.edu/roman/2LPT/

  The 2LPTic code is a 2LPT extension of the N-GenIC code by V. Springel,
  downloadable from:

  http://www.gadgetcode.org/right.html#ICcode

  and distributed under GNU/GPL license.  The part of the 2LPTic code
  we use here is entirely contained in N-GenIC.

*/


#include "pinocchio.h"

#define GRID MyGrids[ThisGrid]

#define BOTTOM 0
#define TOP    1

typedef long long int map_index;
typedef int map_point[2];             // NOTE: this "limits" the plane to a side length of 2^31 points

unsigned int *OLDSEEDTABLE;


unsigned int *SEEDTABLE;
unsigned int *SEED_k0[2], *SEED_x0, *SEED_y0;
map_point     SEED_k0_BL[2];         // bottom-left corner of the seed regions for symmetries on k=0 plane
int           SEED_k0_nx[2];         // x-dimension of the seed regions for symmetries on k=0 plane
int           SEED_k0_ny[2];         // y-dimension of the seed regions for symmetries on k=0 plane


int  generate_seeds_plane(int, map_point, map_point, unsigned, int);
int  quadrant            (map_point, int);  
void dump_seed_plane     (int, map_point, map_point, int[3], int, int, int);



int GenIC_large(int ThisGrid)
{
  int           Nmesh, Nsample, VectorLength;
  int           C[3], local_n[3], local_start[3];
  double        Box, fac;

  VectorLength   = GRID.total_local_size_fft;
  Nmesh          = GRID.GSglobal[_x_];
  Nsample        = GRID.GSglobal[_x_];
  Box            = GRID.BoxSize;
  
  fac = pow(1./Box,1.5);


  dprintf(VDBG, ThisTask, "==== %d in k-space:: (%td, %td) - (%td, %td) - (%td, %td)\n",
	  ThisTask,
	  GRID.GSstart_k[0], GRID.GSstart_k[0] + GRID.GSlocal_k[0],
	  GRID.GSstart_k[1], GRID.GSstart_k[1] + GRID.GSlocal_k[1],
	  GRID.GSstart_k[2], GRID.GSstart_k[2] + GRID.GSlocal_k[2]);
  
  
  // -----------------------------------------------
  // generate seed plane
  // bottom-left and top-right corners for the region of k-plane
  // that this Task must deal with


  // prepare coordinates
  if ( !params.use_transposed_fft )
    // this is the non-transposed order x, y, z    
    for(int i = 0; i < 3; i++)
      C[i] = i;
  else
    {
      if ( internal.tasks_subdivision_dim > 1 )
  	C[0] = _y_, C[1] = _z_, C[2] = _x_;
      else
  	C[0] = _y_, C[1] = _x_, C[2] = _z_;
    }

  // set up starting points and lenghts
  for(int i = 0; i< 3; i++)
    {
      local_n    [i] = GRID.GSlocal_k[i];
      local_start[i] = GRID.GSstart_k[i];
    }


  // actually generate seeds
  {
    map_point     bottom_left, top_right;
    map_point     plane_bottom_left, plane_top_right;

    bottom_left[_x_] = GRID.GSstart_k[_x_];
    bottom_left[_y_] = GRID.GSstart_k[_y_];
    
    top_right[_x_]   = GRID.GSstart_k[_x_] + GRID.GSlocal_k[_x_];
    top_right[_y_]   = GRID.GSstart_k[_y_] + GRID.GSlocal_k[_y_];

    if( !internal.large_plane )
      {
	plane_bottom_left[_x_] = 0;
	plane_bottom_left[_y_] = 0;

	plane_top_right[_x_] = Nmesh;
	plane_top_right[_y_] = Nmesh;
      }
    else
      {
	plane_bottom_left[_x_] = bottom_left[_x_];
	plane_bottom_left[_y_] = bottom_left[_y_];
	plane_top_right[_x_]   = top_right[_x_];
	plane_top_right[_y_]   = top_right[_y_];
      }

    int ret = generate_seeds_plane(ThisGrid, plane_bottom_left, plane_top_right, Nmesh, (local_start[_z_] == 0));

    if( ret > 0 )
      return ret;
    
    if ( internal.dump_seedplane )
      {
	
	if ( params.use_transposed_fft && internal.tasks_subdivision_dim > 1 )
	  C[1] = _x_;
	dump_seed_plane(ThisGrid, bottom_left, top_right, C, local_n[_x_], GRID.GSstart_k[_z_], Nmesh);
      }
  }

  // -----------------------------------------------	

  

  
  // the vector is initialized to zero
  memset( kdensity[ThisGrid], 0, VectorLength * sizeof(double) );

  int           Nmesh_2, Nmesh_odd;
  gsl_rng      *k0_generator;

  Nmesh_2   = Nmesh / 2;
  Nmesh_odd = (Nmesh % 2);  


  // initialize the pseudo-random number chain
  
  gsl_rng_set(random_generator, params.RandomSeed);
  
  // initialize a second random generator, used on the plane k = 0
  
  k0_generator = gsl_rng_alloc(gsl_rng_ranlxd1);

  // --------------------------------------------------------------
  // this is the main loop on grid modes
  //--------------------------------------------------------------
  
  for( int i = 0; i < local_n[_x_]; i++ )
    {

      // ii is the x grid coordinate
      int ii = local_start[_x_] + i;
      if ( ii == Nmesh_2 )
	continue;

      double kvec[3];
      
      if ( ii < Nmesh_2 )
	kvec[0] = ii * 2 * PI / Box;
      else
	kvec[0] = -(Nmesh - ii) * 2 * PI / Box;

      double kmag2_i = kvec[0] * kvec[0];
      int iii = ii;
      
      
      for(int j = 0; j < local_n[_y_]; j++)	    
	{

	  // jj is the y grid coordinate
	  int jj = local_start[_y_] + j;
	  if(jj == Nmesh_2)
	    continue;
	  int jjj = jj;
	  
	  if(jj < Nmesh_2)
	    kvec[1] = jj * 2 * PI / Box;
	  else
	    kvec[1] = -(Nmesh - jj) * 2 * PI / Box;

	  double kmag2_ij = kmag2_i + kvec[1] * kvec[1];

	  unsigned int seed;

	  // initialize the random generator with the seed found in the
	  // seedtable at the (i,j) coordinate in grid coordinates
	  if(internal.large_plane)
	    {
	      gsl_rng_set(random_generator, SEEDTABLE[j * local_n[_x_] + i]);
	      seed = SEEDTABLE[j * local_n[_x_] + i];
	    }
	  else
	    {
	      gsl_rng_set(random_generator, SEEDTABLE[jj * Nmesh + ii]);
	      seed = SEEDTABLE[jj * Nmesh + ii];
	    }

	  // since the chain of random numbers is subsequent, to
	  // reproduce all the time the same sequence, this
	  // process must generate all the random numbers on the
	  // k column that are below its starting k coordinate
	  
	  double phase;
	  if(local_start[_z_] > 0 && local_start[_z_] < Nmesh_2)
	    for(int RR = 0; RR < local_start[_z_]; RR++)
	      {
	  	phase = gsl_rng_uniform(random_generator);
	  	do
	  	  phase = gsl_rng_uniform(random_generator);
	  	while(phase == 0);
	      }

	  // inner loop
	  for(int k = 0; (k < local_n[_z_]) && (k+local_start[_z_] < Nmesh_2); k++)
	    {
	      // kk is the z grid coordinate
	      int kk = local_start[_z_] + k;

	      // generate phase and amplitude
	      phase = gsl_rng_uniform(random_generator) * 2 * PI;
	      double ampl;
	      do
		ampl = gsl_rng_uniform(random_generator);
	      while(ampl == 0);	      	      

	      // blind points on the cube
	      if(ii == 0 && jj == 0 && kk == 0)
		continue;
	      if(kk == Nmesh_2)
		continue;
	      
	      if(kk < Nmesh_2)
		kvec[2] = kk * 2 * PI / Box;
	      else
		kvec[2] = -(Nmesh - kk) * 2 * PI / Box;

	      double kmag2_local = kmag2_ij + kvec[2] * kvec[2];
	      double kmag = sqrt(kmag2_local);

	      if(kmag * Box / (2 * PI) > NYQUIST * Nsample / 2)
		continue;		  	      

	      double p_of_k = PowerSpectrum(kmag);	      

	      double sign = 1.0;	      
	      int addr_j = j;	     

	      // special symmetries on the plane k = 0
	      if( kk == 0 )
		{

		  // some points are left empty
		  if( (ii == 0) && (jj == Nmesh_2 || jj == Nmesh_2 + Nmesh_odd) )
		    continue;
		  if( (ii == Nmesh_2) || (ii == Nmesh_2 + Nmesh_odd) )
		    continue;
	      
		  // the half-plane for ii > Nmesh_2 takes the seeds
		  // from points in the other half-plane
		  if( ( ii > Nmesh_2 ) ||
		      ( ii == 0 && jj > Nmesh/2 ) )
		    {
		      // symmetries on j
		      jjj = Nmesh -jj;
		      if(jjj == Nmesh)
			jjj = 0;

		      if(Nmesh_odd && jj == Nmesh_2+1)
			{
			  jjj = Nmesh_2+1;
			  addr_j = Nmesh_2;
			}

		      // symmetries on i
		      if(ii > Nmesh_2)
			iii = Nmesh - ii;

		      
		      //unsigned int seed;
		      if(internal.large_plane)
			// each task owns only a region of the seeds' plane
			{
			  if(jjj == 0)
			    // the line below would read
			    //   seed = SEED_y0[iii];
			    // however, since SEED_y0 starts from 0 we
			    // have to displace accordingly to i position
			    // reversed because of the requested symmetry
			    {
			      seed = SEED_y0[local_n[_x_] - i];
			    }
			  else if(iii == 0)
			    // the line below would read
			    //   seed = SEED_x0[jjj];
			    // however, since SEED_y0 starts from 0 we
			    // have to displace accordingly to j position
			    // reversed because of the requested symmetry
			    {
			      seed = SEED_x0[local_n[_y_] - j];
			    }
			  else
			    {
			      map_point V;
			      V[_x_] = ii;
			      V[_y_] = jj;
			      int Q = quadrant(V , Nmesh/2) / 2;
			      
			      seed = SEED_k0[Q][(jjj - SEED_k0_BL[Q][_y_]) * SEED_k0_nx[Q] + (iii - SEED_k0_BL[Q][_x_]) ];
			    }
			}
		      else
			{
			  // each task owns the whole seeds plane
			  seed = SEEDTABLE[jjj * Nmesh + iii];
			}

		      
		      sign = -1.0;

		      // re-initialize the spare random generator to the right seed
		      gsl_rng_set(k0_generator, seed);
		      phase = gsl_rng_uniform(k0_generator) * 2 * PI;
		      do
			ampl = gsl_rng_uniform(k0_generator);
		      while(ampl == 0);
		      
		    }
		}

#ifndef NO_RANDOM_MODULES
	      p_of_k *= -log(ampl);
#endif
	      double delta = fac * sqrt(p_of_k);

	      int addr;
	      
	      // calculate the storing address in local coordinates
	      if(!params.use_transposed_fft)
		addr    = 2*((i * local_n[_y_] + addr_j) * local_n[_z_] + k);
	      else
		{
		  // the addressing here is obviously awful and inefficient because
		  //   when the transposed fft are used the memory ordering is
		  //   row-major in z, y, x instead of x, y, z.
		  // 
		  // however the nested loops have been designed for x, y, z order
		  //   that is "natural" because of the z=0 symmetries.
		  //
		  // due to this fact, for the sake of a clearer coding, we prefer
		  //   not to write a different version of the loops, taking into
		  //   account that this part of the code is not expèected to be
		  //   dominant

		  if(internal.tasks_subdivision_dim > 1)
		    addr    = 2*((addr_j * local_n[_z_] + k) * local_n[_x_] + i);
		  else
		    addr    = 2*((addr_j * local_n[_x_] + i) * local_n[_z_] + k);
		}
	      
	      kdensity[ThisGrid][addr    ] =  delta * cos(phase);
	      kdensity[ThisGrid][addr + 1] =  sign * delta * sin(phase);

	    }

	}
    }

  MPI_Barrier(MPI_COMM_WORLD);
  
  gsl_rng_free(k0_generator);
  
  if(SEED_k0[0] != NULL)
    free(SEED_k0[0]);
  if(SEED_k0[1] != NULL)
    free(SEED_k0[1]);
  if(SEED_y0 != NULL)
    free(SEED_y0);
  if(SEED_x0 != NULL)
    free(SEED_x0);
  free(SEEDTABLE);


  /* This is needed to achieve proper normalization for pinocchio */
  
  fac=pow((double)Nmesh, 3.0);
  unsigned int myVectorLength = VectorLength - VectorLength % 8;
  int i;
  for (i = 0; i < myVectorLength; i += 8)
    {
      kdensity[ThisGrid][i  ] *= fac;
      kdensity[ThisGrid][i+1] *= fac;
      kdensity[ThisGrid][i+2] *= fac;
      kdensity[ThisGrid][i+3] *= fac;
      kdensity[ThisGrid][i+4] *= fac;
      kdensity[ThisGrid][i+5] *= fac;
      kdensity[ThisGrid][i+6] *= fac;
      kdensity[ThisGrid][i+7] *= fac;
    }
  for( ; i < VectorLength; i++)
    kdensity[ThisGrid][i] *= fac;
  
  // dump the kdensity cube on a ascii file
  // this part shall be removed in production,
  //   or made optional

  if(internal.dump_kdensity)
    dump_cvector(kdensity[ThisGrid],
		 params.use_transposed_fft*internal.tasks_subdivision_dim,
    		 MyGrids[ThisGrid].GSglobal[_x_],
    		 MyGrids[ThisGrid].GSlocal_k,
    		 MyGrids[ThisGrid].GSstart_k, "my.dkdensity", 0);
  

  return 0;
}


// =====================================================================================


#define max(x,y) ((x) > (y) ? (x) : (y))

#define BL 0  // bottom-left
#define TR 1  // top-right

void      symmetry                 (map_point *, map_point *, int);
map_index get_map                  (map_point);
int       cmp_map_points           (const void *, const void *);
int       generate_seeds_subregion (map_point, map_point, uint, uint **, uint, map_index *);
void      transpose_subregion      (int, uint, map_point *, map_point *);
int       get_plane_subregions     (uint, map_point, map_point, map_point [4][2], uint [4], uint[4]);
int       copy_seeds_subregion     (map_point, map_point, uint *, uint, uint);




int generate_seeds_plane(int ThisGrid, map_point bottom_left, map_point top_right, unsigned Nmesh, int k0symmetries)
{

  /* ---------------------------------------
   *
   *  STEP 1
   *
   *  generate old seed table if requested
    --------------------------------------- */

  
  if(internal.mimic_original_seedtable)
    {

      OLDSEEDTABLE = (unsigned int*) calloc( GRID.GSglobal[_x_] * GRID.GSglobal[_y_], sizeof(unsigned int));
      
      if( !internal.large_plane )
	SEEDTABLE = OLDSEEDTABLE;

      gsl_rng_set(random_generator, params.RandomSeed);

      for(int i = 0; i < Nmesh / 2; i++)
	{
	  // proceed inwards from the 4 corners
	  
	  // 0, 0
	  for(int j = 0; j < i ; j++)
	    OLDSEEDTABLE[j * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	  for(int j = 0; j < i + 1; j++)
	    OLDSEEDTABLE[i * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);	  
	  
	  // Nmesh, 0
	  for(int j = 0; j < i; j++)
	    OLDSEEDTABLE[j * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	  for(int j = 0; j < i + 1; j++)
	    OLDSEEDTABLE[i * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	  
	  // 0, Nmesh
	  for(int j = 0; j < i; j++)
	    OLDSEEDTABLE[(Nmesh - 1 - j) * Nmesh + i] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	  for(int j = 0; j < i + 1; j++)
	    OLDSEEDTABLE[(Nmesh - 1 - i) * Nmesh + j] = 0x7fffffff * gsl_rng_uniform(random_generator);	  
	  
	  
	  // Nmesh, Nmesh
	  for(int j = 0; j < i; j++)
	    OLDSEEDTABLE[(Nmesh - 1 - j) * Nmesh + (Nmesh - 1 - i)] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	  for(int j = 0; j < i + 1; j++)
	    OLDSEEDTABLE[(Nmesh - 1 - i) * Nmesh + (Nmesh - 1 - j)] = 0x7fffffff * gsl_rng_uniform(random_generator);
	  
	}
	
    }  

  if( (!internal.large_plane) && internal.mimic_original_seedtable )
    return 0;

  
  /* ---------------------------------------
   *
   *  STEP 2
   *
   *  continue if either a large plane is
   *  requested or the new seed table is to
   *  be generated
   --------------------------------------- */


  int           N_subregions;
  map_index     error;  
  map_point     subregions[4][2] = { {{0, 0}, {0, 0}}, {{0, 0}, {0, 0}} };
  uint          subregions_offsets[4] = {0 ,0, 0, 0}, subregions_starts[4] = {0, 0, 0, 0};

  
  // individuate the sub-regions of this task's domain to map them onto spiral-plane
  
  N_subregions = get_plane_subregions(Nmesh, bottom_left, top_right, subregions, subregions_offsets, subregions_starts);
  dprintf(VDBG, ThisTask, "Task %d: subregions done\n", ThisTask); fflush(stdout);

  
  // show some infos
  
  if(internal.verbose_level > VDBG)
    for(int T = 0; T < NTasks; T++)
      {
	if( ThisTask == T )
	  for(int i = 0; i < N_subregions; i++)
	    {
	      dprintf(VDBG, ThisTask, "\ttask: %3d\n"
		      "\t\tReg %d: bottom left : %d, %d\n"""
		      "\t\tReg %d: top right : %d, %d\n"
		      "\t\tReg %d: st/off: %u, %u\n ", ThisTask,
		      i, subregions[i][BL][0], subregions[i][BL][1], i, subregions[i][TR][0], subregions[i][TR][1],
		      i, subregions_starts[i], subregions_offsets[i]);
	    }
	
	MPI_Barrier(MPI_COMM_WORLD);
      }


  /* ---------------------------------------
   *
   *  STEP 3
   *
   *  either copy old seeds or generate seeds 
   *  accordingly to spiral order in each
   *  subregion
   * --------------------------------------- */

  {
    uint      ystride, *start_pointer;  
    map_index Npoints;


      // calculate the total length of the whole region along y
    ystride = subregions[N_subregions-1][TR][_y_] - subregions[0][BL][_y_];
    
    // calculate how many points will be generated
    Npoints = ystride * (subregions[N_subregions-1][TR][_x_] - subregions[0][BL][_x_]);
    
    // allocate memory for all the seeds
    
    if( ( SEEDTABLE = calloc( sizeof(uint), (size_t)Npoints) ) == NULL )
      {
	dprintf(VXERR, ThisTask, "Unable to allocate %lld bytes\n", Npoints * sizeof(uint));   
	return 1;
      }

    if( internal.mimic_original_seedtable )
      // copy same region from the old seed plane in the right place
      
      for(int i = 0; i < N_subregions; i++)
	{
	  start_pointer = &(SEEDTABLE[subregions_starts[i]]);
	  dprintf(VDBG, ThisTask, "\t\t ccc task %d copying seed in region %d "
		  ":: (%d, %d) -> (%d, %d) starting at %u, offset %u\n",
		  ThisTask, i, subregions[i][BOTTOM][_x_], subregions[i][BOTTOM][_y_],
		  subregions[i][TOP][_x_], subregions[i][TOP][_y_], subregions_starts[i], subregions_offsets[i]);
	  
	  copy_seeds_subregion(subregions[i][BOTTOM], subregions[i][TOP], start_pointer, Nmesh, subregions_offsets[i]);		
	}
    else
      for( int i = 0; i < N_subregions; i++ )
	{
	  
	  map_point transposed_subregion[2];
	  
	  transposed_subregion[BOTTOM][_x_] = subregions[i][BOTTOM][_x_];
	  transposed_subregion[BOTTOM][_y_] = subregions[i][BOTTOM][_y_];
	  transposed_subregion[TOP][_x_]    = subregions[i][TOP][_x_];
	  transposed_subregion[TOP][_y_]    = subregions[i][TOP][_y_];
	  
	  // transpose subregion in spiral's coordinates
	  transpose_subregion(0, Nmesh, &transposed_subregion[0], &transposed_subregion[1]);
	  
	  start_pointer = &SEEDTABLE[subregions_starts[i]];
	  
	  // generate the seeds and put them in memory at the right place
	  generate_seeds_subregion(transposed_subregion[0], transposed_subregion[1],
				   params.RandomSeed, &start_pointer, subregions_offsets[i], &error);
	}
    
  }
  
  dprintf(VDBG, ThisTask, "Task %d: seeds generated\n", ThisTask); fflush(stdout);



/* ---------------------------------------
   *
   *  STEP 4
   *
   *  if a large plane is used AND this
   *  task's domain include k=0 symmetries
   *  populate buffers to fullfill those
   *  symmetries
   * --------------------------------------- */

  
  if(internal.large_plane && k0symmetries)
    {
  
      // now check for additional symmetries on plane k = 0:
      //
      // - quadrant 1 and 3 just pick seeds from regions in
      //   quadrants 2 and 0 respectively, symmetrical with
      //   respect to the centre (Nmesh/2, Nmesh/2)
      // - row x=0 in quadrant 1 picks seed from row x=0
      //   in quadrant 0 (symmetrically to y=Nmesh/2)
      // - column y=0 in quadrant 2 picks seeds from column y=0
      //   in quadrant 0 (symmetrically to x=Nmesh/2)
      //
      //                    ^
      //              Q 2   | Q 3
      //                    |
      //           ------------------>
      //                    |
      //              Q 0   | Q 1
      //
      // 

      map_point    SBL, STR; // Symmetric Bottom-Left, Top-Right
      unsigned int Np;

  
      for(int i = 0; i < N_subregions; i++)
	{
	  int Q             = quadrant(subregions[i][0], Nmesh/2);
	  int SEED_k0_index = Q / 2;
	  
	  if(Q > 0)
	    {
	      switch(Q)
		{
	      
		case(1):  // row y = 0 and/or plane symmetry

		  if(subregions[i][BOTTOM][_y_] == 0)
		    // this subregion in quadrant 1 contains the stripe y = 0
		    {
		      // set-up corners for seed sub-region generation
		      SBL[_y_] = 0, SBL[_x_] = Nmesh - subregions[i][TOP][_x_];// + 1;
		      STR[_y_] = 1, STR[_x_] = Nmesh - subregions[i][BOTTOM][_x_] + 1;

		      // find horizontal extension in this quadrant
		      Np = STR[_x_] - SBL[_x_];
		  
		      // allocate memory for seeds in this stripe
		      SEED_y0 = (unsigned int*)malloc(sizeof(unsigned int) * Np);
		      
		      if(!internal.mimic_original_seedtable)
			{
			  // transpose subregion in spiral's coordinates
			  transpose_subregion(0, Nmesh, &SBL, &STR);

			  generate_seeds_subregion(SBL, STR, params.RandomSeed, &SEED_y0, 0, 0);
			}
		      else
			copy_seeds_subregion(SBL, STR, SEED_y0, Nmesh, 0);

		    }

		  // NOTE ---> the previous lines dealed with the y=0 border of Q1;
		  //           now continue with the operations for plane symmetry.
		  //           Since they are the the same than for the plane symmetry
		  //           in case(3), we simply avoid the usual break and let it
		  //           continue below with the code for case(3)
		  
		case(3): // plane symmetry

		  // set-up corners for seed sub-region generation
		  SBL[_x_] = subregions[i][BOTTOM][_x_], SBL[_y_] = subregions[i][BOTTOM][_y_];
		  STR[_x_] = subregions[i][TOP   ][_x_], STR[_y_] = subregions[i][TOP   ][_y_];

		  // take the symmetrical region
		  symmetry(&SBL, &STR, Nmesh);
		  SEED_k0_BL[SEED_k0_index][_x_] = STR[_x_], SEED_k0_BL[SEED_k0_index][_y_] = STR[_y_];
		  SEED_k0_nx[SEED_k0_index]      = SBL[_x_] - STR[_x_];
		  SEED_k0_ny[SEED_k0_index]      = SBL[_y_] - STR[_y_];
		  
		  // transpose into spiral-plane region;
		  // note that taking the symmetrical with respect to the centre
		  // the bottom-left and top-right corners are swapped
		  if(!internal.mimic_original_seedtable)
		    transpose_subregion(0, Nmesh, &STR, &SBL);

		  // calculate how many points are in region
		  Np = (SBL[_x_] - STR[_x_]) * (SBL[_y_] - STR[_y_]);

		  // allocate memory
		  SEED_k0[SEED_k0_index] = (unsigned int*)malloc(sizeof(unsigned int) * Np);

		  if(internal.mimic_original_seedtable)
		    copy_seeds_subregion(STR, SBL, SEED_k0[SEED_k0_index], Nmesh, 0);
		  else
		    // generate appropriate seeds
		    generate_seeds_subregion(STR, SBL, params.RandomSeed, &SEED_k0[SEED_k0_index], 0, 0);

		  break;

		case(2): // column x = 0
		  if(subregions[i][BOTTOM][_x_] == 0)
		    {
		      
		      // set-up corners for seed sub-region generation
		      SBL[_x_] = 0, SBL[_y_] = Nmesh - subregions[i][TOP][_y_];
		      STR[_x_] = 1, STR[_y_] = Nmesh - subregions[i][BOTTOM][_y_] + 1;

		      // find vertical extension in this quadrant
		      Np = STR[_y_] - SBL[_y_];
		  
		      // allocate memory for seeds in this stripe
		      SEED_x0 = (unsigned int*)malloc(sizeof(unsigned int) * Np);

		      if(!internal.mimic_original_seedtable)
			{
			  // transpose subregion in spiral's coordinates
			  transpose_subregion(0, Nmesh, &SBL, &STR);

			  generate_seeds_subregion(SBL, STR, params.RandomSeed, &SEED_x0, 0, 0);
			}
		      else
			copy_seeds_subregion(SBL, STR, SEED_x0, Nmesh, 0);
		      
		    }     
		  break;

		  
		}  // ends switch(Q)
	      
	    }  // ends if (Q > 0)
	  
	}  // ends subregions loop
      
    }  // ends k0symmetries


  dprintf(VDBG, ThisTask, "Task %d: symmetries generated\n", ThisTask); fflush(stdout);
  
  return 0;
}




void symmetry(map_point *bottom_left, map_point *top_right, int N)
{
  for(int i = 0; i < 2; i++)
    {
      // both the "+ 1 " terms below account for the upper limits being open
      
      (*bottom_left)[i] = abs(N - (*bottom_left)[i]) + 1;  
      (*top_right)[i] = abs(N - (*top_right)[i] + 1);
      
    }

  return;
}




int quadrant(map_point P, int N_2)
{
  int Q;
  Q = (P[0] >= N_2);
  Q += (P[1] >= N_2) * 2;
  return Q;
}




map_index get_map(map_point P)
// returns the ordinal number of a given point along the spiral
{
  map_index l, d;
  uint mx, my, c;
  
  mx = abs(P[_x_]);
  my = abs(P[_y_]);
  l = 2 * max(mx, my);

  c = (P[_y_] > P[_x_]) + (P[_x_] > 0)*(P[_x_] == P[_y_]);
  
  d = (c) ? l * 3 + P[_x_] + P[_y_] : l - P[_x_] - P[_y_];

  return (l-1)*(l-1) + d;
}


map_index *map;


int cmp_map_points(const void *A, const void *B)
// Comparison function to be called by qsort.
// Performs indexed-sorting: "map" (that is the
// indexed array) must be in its scope
{
  map_index a = map[*(uint*)A];
  map_index b = map[*(uint*)B];

  return (a - b);
}




int generate_seeds_subregion(map_point bottom_left, map_point top_right, uint seed, uint **seed_plane, uint offset, map_index *value)
// this function manages the generation of seeds for a single subregion
// - bottom_left, top_right are the corners of the subregion
// - seed is the initial seed for the random number generator
// - seed_plane is a pointer to the pointer that will host the seeds
// - offset is the offset to be accounted while filling the memory pointed by seed_plane
// - xstride is the xlength of the whole memory matrix pointed by seed_plane, it's used to fill it correctly
// - value is a returning value in case of error. Basically it handles the case in which some allocation fails, and
//   it will return the amount of bytes that were requested.
// we assume that the seed plane is stored in memory startign from bottom_left  

{
  map_index  N, rnd_counter, Ncnt;
  map_point  p;
  uint       xlength, ylength;
  uint      *idxs;

  // calculate the region's dimensions
  xlength = top_right[_x_] - bottom_left[_x_];
  ylength = top_right[_y_] - bottom_left[_y_];
  N = xlength * ylength;

  // Allocate actual memory for the random seeds.
  // The difference with map is that we do not need to store 8-bytes
  // integers, while we needed it to store in the map array the
  // spiral ordered points
  if(*seed_plane == NULL)
    {
      if((*seed_plane = calloc( (size_t)N, sizeof(uint) )) == NULL)
	{
	  if(value != NULL)
	    *value = N * sizeof(uint);
	  return 3;
	}
      offset = 0;
    }
  
  // allocate temporary memory for the random seeds
  if((map = calloc( (size_t)N, sizeof(map_index) )) == NULL)
    {
      if(value != NULL)
	*value = N * sizeof(map_index);
      return 1;
    }

  // allocate temporary memory for the indexes that will be used to
  // sort the points along the spiral
  if((idxs = calloc( (size_t)N, sizeof(uint) )) == NULL)
    {
      if(value != NULL)
	*value = N * sizeof(uint);      
      return 2;
    }

    
  // initialize indexes
  // for the following loop check whether:
  // - leave unroll to the compiler
  // - or explicitly use simd instruction
  for(int i = 0; i < N; i++)
    idxs[i] = i;

  // get the map value from spiral for each point in the region
  for(int j = 0; j < ylength; j++)
    {
      p[_y_] = bottom_left[_y_]+j;

      int j_offset = j*xlength;
      for(int i = 0; i < xlength; i++)
	{
	  p[_x_] = bottom_left[_x_]+i;
	  map[j_offset + i] = get_map(p);
	}
    }

  // index-sort the map
  qsort(idxs, (size_t)N, sizeof(uint), cmp_map_points);

  // generate random numbers following the spiral traversal order

  // set the random number generator
  gsl_rng *R;
  
  R = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(R, seed);
  
  // actually generates pseudo-random number
  rnd_counter = 0;
  Ncnt = 0;
  for(int i = 0; i < N; i++)
    {
      if(rnd_counter < map[idxs[i]]-1)
	for(; rnd_counter < map[idxs[i]]-1; rnd_counter++)
	  gsl_rng_get(R);
      
      map[idxs[i]] = (map_index)gsl_rng_get(R);
      Ncnt++;
      rnd_counter++;
    }
  
  gsl_rng_free(R);  // releases generator
  
  free(idxs);       // releases indexes

  //copy the seeds from map into seed_plane, as uint instead of as long long
  for(uint k = 0, j = 0; j < N; j++)
    {
      if((k > 0) && (j % xlength == 0))
	k += offset;
      (*seed_plane)[k++] = (uint)map[j];
    }
  
  free(map);
  map = 0x0;
  return 0;
}




int copy_seeds_subregion(map_point bottom_left, map_point top_right, uint *seed_plane, uint Nmesh, uint offset)
{
  int jdim, idim;
  
  if((idim = top_right[_x_] - bottom_left[_x_]) == 0)
    idim = 1;
  if((jdim = top_right[_y_] - bottom_left[_y_]) == 0)
    jdim = 1;

  for(int k = 0, j = bottom_left[_y_]; j < top_right[_y_]; j++)
    {
      for(int i = bottom_left[_x_] ; i < top_right[_x_]; i++)
  	seed_plane[k++] = OLDSEEDTABLE[j*Nmesh + i];
    
      k += offset;
    }

  return 0;
}



void transpose_subregion(int mode, uint N, map_point *bottom_left, map_point *top_right)
// note: by construction single regions are supposed to lie within quadrants
{
  uint N_2 = N/2;

  map_point copy_1 = {(*bottom_left)[0], (*bottom_left)[1]};
  map_point copy_2 = {(*top_right)[0], (*top_right)[1]};


  
  if(mode == 0)
    // from plane to spiral
    {     
      if((*bottom_left)[_x_] >= N_2)
	{
	  (*bottom_left)[_x_] = (*bottom_left)[_x_] - N;
	  (*top_right)[_x_]   = (*top_right)[_x_] - N;
	}
      
      if((*bottom_left)[_y_] >= N_2)
	{
	  (*bottom_left)[_y_] = (*bottom_left)[_y_] - N;
	  (*top_right)[_y_]   = (*top_right)[_y_] - N;
	}
    }
  else
    // from spiral back to plane
    {
      if((*bottom_left)[_x_] < 0)
	{
	  (*bottom_left)[_x_] = N + (*bottom_left)[_x_];
	  (*top_right)[_x_]   = N + (*top_right)[_x_];
	}
      
      if((*bottom_left)[_y_] < 0)
	{
	  (*bottom_left)[_y_] = N + (*bottom_left)[_y_];
	  (*top_right)[_y_]   = N + (*top_right)[_y_];
	}
    }

  dprintf(VDBG, ThisTask, "\t\t Task %d: region [%d, %d] -> "
	  "[%d, %d[ transposes into [%d, %d] -> [%d, %d[\n",
	  ThisTask, copy_1[0], copy_1[1], copy_2[0], copy_2[1],
	  (*bottom_left)[0], (*bottom_left)[1], (*top_right)[0], (*top_right)[1]);
  return;
}




int get_plane_subregions(uint Nmesh, map_point bottom_left, map_point top_right, map_point Regions[4][2], uint Offsets[4], uint Starts[4])
{
#define NCOORDS 2

  int  c[NCOORDS] = {0, 1};
  uint Nmesh_2 = Nmesh/2;
  int  Nregions = 1;

  // copy the bottom-left and the top-right corner in the first subregion
  for(int i = 0; i < NCOORDS; i++)
    {
      Regions[0][BL][i] = bottom_left[i];
      if(bottom_left[i] > Nmesh)
	return -1;
    }
  
  for(int i = 0; i < NCOORDS; i++)
    {
      Regions[0][TR][i] = top_right[i];
      if(top_right[i] > Nmesh)
	return -1;
    }

  for(int i = 0; i < 4; i++)
    Offsets[i] = 0, Starts[i] = 0;  
  
  
  // sub-divide for each coordinate if the extension
  // along that coordinate encompassed Nmesh/2
  //

  int xlen = top_right[_x_] - bottom_left[_x_];

  for(int k = 0; k < NCOORDS; k++)
    {
      int m = 1;
      
      // cycle over all existing sub-regions
      for(int R = Nregions-1; R >= 0; R--)
	
	// if a region overlaps N/2 along the current coordinate, we need to sub-divide
	if(Regions[R][BL][c[k]] < Nmesh_2 &&
	   Regions[R][TR][c[k]] > Nmesh_2)
	  {
	    int j = Nregions + R;

	    for(int i = 0; i < NCOORDS; i++)
	      Regions[j][BL][c[i]] = Regions[R][BL][c[i]];
	    Regions[j][BL][c[k]] = Nmesh_2;
	    
	    for(int i = 0; i < NCOORDS; i++)
	      Regions[j][TR][c[i]] = Regions[R][TR][c[i]];
	    Regions[R][TR][c[k]] = Nmesh_2;

	    if(k == 0)
	      {
		Offsets[j] = Nmesh_2 - Regions[R][BL][c[k]];
		Offsets[R] = Regions[j][TR][c[k]] - Nmesh_2;
	      }
	    else
	      Offsets[j] = Offsets[R];

	    Starts[j] = (Regions[j][BL][c[_y_]] - bottom_left[_y_]) * xlen +
	      (Regions[j][BL][c[_x_]] - bottom_left[_x_]);

	    m = 2;
	  }

      Nregions *= m;
    }
  

  return Nregions;  
}


void dump_seed_plane(int ThisGrid, map_point bottom_left, map_point top_right, int C[3], int n_x, int z0, int Nmesh)
// dump the seed plane
// NOTE:: the transposition is taken into account through the rotation array C
//        Then, the seedplane is always written in the "right" orientation, i.e. with x and y axis in their
//        "natural" orientation
{
  char  buffer[400];

  sprintf(buffer, "seedplane_partial_%d", ThisTask);
  dprintf(VDIAG, ThisTask, "\t\tTask %d writing seedplane: %d, %d -> %d, %d\n", ThisTask, 
	  bottom_left[_x_], bottom_left[_y_], top_right[_x_], top_right[_y_]);
  
  
  for(int tt = 0; tt < NTasks; tt++)
    {
      if( ThisTask == tt)
	{
	  if( z0 == 0)
	    // write data only if this task owes the z = 0 plane
	    // this avoids to write the same portion of plane multiple times
	    {
	      FILE *file = fopen(buffer, "w");
	      int k, offset;
	      
	      if(internal.large_plane == 0)
		{
		  k = bottom_left[_y_] * Nmesh + bottom_left[_x_];
		  offset = Nmesh - n_x;
		}
	      else
		k = offset = 0;


	      for(int jj = bottom_left[_y_]; jj < top_right[_y_]; jj++)
		{
		  if( internal.dump_seedplane == 1 )
		    // write the plane in its "natural order"
		    {
		      unsigned int line = jj*Nmesh;
		      for(int ii = bottom_left[_x_]; ii < top_right[_x_]; ii++)
			fprintf(file, "%d %d %d %u\n", ii, jj, line + ii, SEEDTABLE[k++]);
		    }
		  else
		    // write the plane in its transposed order
		    {
		      for(int ii = bottom_left[_x_]; ii < top_right[_x_]; ii++)
			fprintf(file, "%d %d %d %u\n", jj, ii, ii*Nmesh + jj, SEEDTABLE[k++]);
		    }
		  k += offset;
		}
	      
	      fclose(file);
	    }
	  else 
	    dprintf(VDBG, ThisTask, "\t\t\tTask %d skips because z0 == %d\n", ThisTask, z0);
	}
      MPI_Barrier(MPI_COMM_WORLD);
    }
  
  if(ThisTask == 0)
    {
      int err;
      char filename[200];
      sprintf(filename, "seed_plane.Ng_%d_Nt_%d.dat", Nmesh, NTasks);
      sprintf(buffer, "cat seedplane_partial_* | sort -k3 -n > %s", filename);

      err = system(buffer);
      if(err != 0)
	dprintf(VXX, 0, "something went wrong while writing %s file\n", filename);
      
      err = system("rm -f seedplane_partial_*");
      if(err != 0)
	dprintf(VXX, 0, "something went wrong while removing temporary files seedplane_partial_*\n");
    }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  return;
}


//double VarianceOnGrid(int ThisGrid, double Redshift) //, double ThisRadius)
//{
//
//  int i,j,k,ii,jj,kk;
//  double fundamental = 2. * PI / MyGrids[ThisGrid].BoxSize;
//  double nyquist = NYQUIST * MyGrids[ThisGrid].GSglobal[_x_]/2 * fundamental;
//  double d3k=pow(fundamental,3.);
//  double D, kmag;
//  double Variance=0.0;
//
//  for (i=0; i<=MyGrids[ThisGrid].GSglobal[_x_]; i++)
//    {
//      ii = i-MyGrids[ThisGrid].GSglobal[_x_]/2;
//      for (j=0; j<=MyGrids[ThisGrid].GSglobal[_y_]; j++)
//	{
//	  jj = j-MyGrids[ThisGrid].GSglobal[_y_]/2;
//	  for (k=0; k<=MyGrids[ThisGrid].GSglobal[_z_]; k++)
//	    {
//	      kk = k-MyGrids[ThisGrid].GSglobal[_z_]/2;
//	     
//	      
//	      kmag=sqrt((double)(ii*ii+jj*jj+kk*kk))*fundamental;
//	      if (kmag>0.0 && kmag <= nyquist)
//		{
//		  //w = WindowFunction(kmag * ThisRadius);
//		  D = GrowingMode(Redshift,kmag);
//		  Variance += PowerSpectrum(kmag) * D*D * d3k;  // * w*w
//		}
//	    }
//	}
//    }
//
//  Variance /= pow(2.*PI,3.);
//
//  return Variance;
//}
//
