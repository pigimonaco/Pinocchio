#!/bin/bash

ml gcc/12.2 binutils/2.41 nvhpc/24.3 openmpi/4.1.6--nvhpc--24.3

# edit the libraries' path
LIB=/leonardo_scratch/fast/CNHPC_1498509/lib/nvhpc-23.11
FTTW_LIB=${LIB}/fftw/fftw-3.3.10/lib
GSL_LIB=${LIB}/gsl/gsl-2.7.1/lib
PFFT_LIB=${LIB}/pfft/pfft/lib
PMT_LIB=/leonardo_scratch/fast/CNHPC_1498509/lib/pmt/local/lib64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FTTW_LIB}:${GSL_LIB}:${PFFT_LIB}:${PMT_LIB}

export OMP_TARGET_OFFLOAD=mandatory
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=ACTIVE