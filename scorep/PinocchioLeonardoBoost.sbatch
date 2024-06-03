#!/bin/bash

#SBATCH --job-name=PinocchioScorep
#SBATCH --partition=boost_usr_prod
#SBATCH --account=IscrC_SCGPCT
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=PinocchioScorep-%j.out
#SBATCH --error=PinocchioScorep-%j.err
#SBATCH --time=3:00:00

module load scorep/8.1--openmpi--4.1.4--nvhpc--23.1-cuda-11.8
module load nvhpc/23.11
# gsl, fftw, pfft have been compiled locally using nvhpc
# src/Makefile has been set up accordingly

# edit the libraries' path
LIB=/leonardo/home/userexternal/dgoz0000/lib/nvhpc-23.11
FTTW_LIB=${LIB}/fftw/fftw-3.3.10/lib
GSL_LIB=${LIB}/gsl/gsl-2.7.1/lib
PFFT_LIB=${LIB}/pfft/pfft/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${FTTW_LIB}:${GSL_LIB}:${PFFT_LIB}:/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/binutils-2.41-gzwp3ai6qpc44qi4ajukobnetdxqd6za/lib/

############################################## TO EDIT #############################################################
####################################################################################################################

#------- set the working dir
export WORKDIR=/leonardo/home/userexternal/dgoz0000/FAST_IscrC_SCGPCT/Pinocchio

#------ switch on/off profiling
PROFILING=1

#------ switch on/off tracing
# Tracing is meant to be run after the profiling
TRACING=0

#------ set the scorep memory if tracing is enable
export SCOREP_TOTAL_MEMORY=4097kB # set this value according to estimated memory requirements during profiling analysis

# set the compiler
COMPILER=mpicc

# set the machine name on the src/Makefile accordingly
SYSTEM=LeonardoBoost

# executable basename
BASE=Pinocchio

# executable with/without OMP or Debug support
# EXEC=( "${BASE}" "${BASE}OMP" )
EXEC=( "${BASE}OMP" )

# Select the PAPI metrics (empty means do not use PAPI counters)
# PAPI_METRIC=PAPI_DP_OPS,PAPI_L3_TCM
PAPI_METRIC=

# Number of MPI tasks (array)
NTASKS=( 8 )

# MPI mapping list (e.g. socket, numa, core)
MPI_MAP_BY=( numa )

# Number of OMP threads per MPI task (array)
OMP_THR=( 2 4 )

# Parameter file
PARAMFILE=${WORKDIR}/example/parameter_file_profiling

# Outputs (name of file with required output redshifts)
OUTPUTS=${WORKDIR}/example/outputs_scorep

# Profiling/tracing directory (output directory)
OUT_DIR=${WORKDIR}/profiling_256_par_256_box


####################################################################################################################
####################################################################################################################
#
#
#
#
#
export SCOREP_PATH=${WORKDIR}/scorep

if [ ! -d ${SCOREP_PATH} ]
then
    printf "\n\t ${SCOREP_PATH} does not exist... aborting... \n\n"
    exit 1
fi

# load the run configuration
source ${SCOREP_PATH}/run.sh

printf "\n\t Running on ${NODES} nodes \n"
printf "\n\t Running Pinocchio on hosts: ${SLURM_NODELIST}"

# get the actual topology
source ${SCOREP_PATH}/topology.sh ${NODES}

if [ ${PROFILING} -eq 1 ]
then
    printf "\n\t Performing the profiling of Pinocchio \n\n"

    # compile the code with score-p support
    source ${SCOREP_PATH}/compile.sh ${COMPILER}

    # perform the profiling
    source ${SCOREP_PATH}/profiling.sh
fi

if [ ${TRACING} -eq 1 ]
then
    printf "\n\t Performing the tracing of Pinocchio \n\n"
    
    source ${SCOREP_PATH}/tracing.sh
fi

printf "\n\t END OF GAME!!! \n\n"

exit 0

