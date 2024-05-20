#!/bin/bash

# This script is meant to set the parameters
# to compile and perform code profiling and
# tracing using Score-p tool

# set the compiler
COMPILER=mpicc

# set the machine name on the src/Makefile (SYSTYPE)
SYSTEM=LeonardoBoost

# executable basename
BASE=Pinocchio
# executable with/without OMP or Debug support
EXEC=( "${BASE}" "${BASE}OMP" )

# Select the PAPI metrics (empty means do not use PAPI counters)
# PAPI_METRIC=PAPI_DP_OPS,PAPI_L3_TCM
PAPI_METRIC=

# Number of MPI tasks (array)
NTASKS=( 16 )

# MPI mapping list (e.g. socket, numa, core)
MPI_MAP_BY=( numa )

# Number of OMP threads per MPI task (array)
OMP_THR=( 2 )

# Number of nodes (set automatically by Slurm)
if [ ${#SLURM_NNODES} -eq 0 ]
then
    NODES=1
else
    NODES=${SLURM_NNODES}
fi
