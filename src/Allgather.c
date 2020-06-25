#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <mpi.h>

#define N_default 1000000

int main( int argc, char **argv )
{
  MPI_Status Status;
  int        Ntasks, Me;
  
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &Ntasks );
  MPI_Comm_rank( MPI_COMM_WORLD, &Me );

  int  N = N_default;
  int *Table;
  
  if ( argc > 1 )
    N = atoi( *(argv+1) );

  if ( (Table = (int*)calloc( N, sizeof(int) )) == NULL )
    {
      printf("Not enough memory to allocate %d ints\n", N );
      MPI_Finalize();
      return 1;
    }

  int mystart, mysize, myend;
  int *recvcounts;
  int *recvdispls;

  int bunch     = N / Ntasks;
  int remainder = N % Ntasks;
  
  recvcounts = (int*)calloc( Ntasks, sizeof(int) );
  recvdispls = (int*)calloc( Ntasks, sizeof(int) );

  for( int i = 0; i < Ntasks; i++ )
    {
      recvcounts[i] = bunch + (i < remainder);
      recvdispls[i] = bunch*i + ((remainder > i)? i : remainder );
    }

  myend     = recvdispls[ Me ] + recvcounts[ Me ];
  
  for ( int i = 0; i < Ntasks; i++ )
    {
      if ( Me == i )
	printf( "%d: my start is at %d, my size is %d, my end will be %d\n", Me, recvdispls[ Me ], recvcounts[ Me ], myend-1 );
      MPI_Barrier( MPI_COMM_WORLD );
    }

  for ( int i = recvdispls[ Me ]; i < myend; i++ )
    Table[i] = Me;

  MPI_Allgatherv( MPI_IN_PLACE, mysize, MPI_INT, Table, recvcounts, recvdispls, MPI_INT, MPI_COMM_WORLD );
  

  // check that everything is correct
  unsigned err = 0;
  for( int i = 0, j = 0; i < N; i++ )
    {
      j += (i == recvdispls[j+1]);
      err += (Table[i] != j);
    }

  unsigned err_tot;
  MPI_Reduce( &err, &err_tot, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

  if ( (Me == 0) && (err_tot > 0) )
    printf("something incorrect\n");

      
  
  free( recvdispls );
  free( recvcounts );
  free( Table );
  
  MPI_Finalize( );
  return 0;
}
