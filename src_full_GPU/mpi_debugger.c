#ifdef MPI_ATTACH_DEBUGGER

#include <unistd.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MASTER 0

void mpi_attach_debugger(MPI_Comm comm)
{
  // get MPI rank
  int rank = -1;
  MPI_Comm_rank(comm, &rank);

  // get MPI communicator size                                                                                                                                                           
  int Nranks = -1;
  MPI_Comm_size(comm, &Nranks);
  // get the hostname

  char hostname[MPI_MAX_PROCESSOR_NAME];
  int resultlen = -1;
  MPI_Get_processor_name(hostname, &resultlen);

  // get the pid of the current process
  const pid_t pid = getpid();

  for (int task=0 ; task<Nranks ; task++)
    {
      if (task == rank)
	{
          printf("\n\t Task: %d - pid: %i - hostname: %s",
                 rank, pid, hostname);
          fflush(stdout);
        }
      MPI_Barrier(comm);
    }

  // master rank collects all pids
  pid_t *all_pids = NULL;
  if (rank == MASTER)
    all_pids = (pid_t *)malloc(Nranks * sizeof(pid_t));

  MPI_Gather(&pid,     sizeof(pid_t), MPI_BYTE,
             all_pids, sizeof(pid_t), MPI_BYTE,
             MASTER, comm);

  if (rank == MASTER)
    {
      FILE *fd = fopen("pid_list_for_debugger.txt", "w");
      for (int task=0 ; task<Nranks ; task++)
        fprintf(fd, "%i \n", all_pids[task]);
      fclose(fd);

      free(all_pids);

      printf("\n\n\n\t pid_list_for_debugger.txt written \n\n");
      fflush(stdout);
    }

  MPI_Barrier(comm);

  volatile int continue_run = 0;
  while (continue_run == 0) /* continue_run needsto be set to 1 ("set var
                               continue_run = 1") by debugger to continue */
    {
      sleep(1);
    }

  return;
}

#endif // MPI_ATTACH_DEBUGGER
