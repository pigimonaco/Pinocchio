#!/bin/bash

#SBATCH --job-name=PinocchioScorep
#SBATCH --partition=pleiadi
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --output=PinocchioScorep-%j.out
#SBATCH --error=PinocchioScorep-%j.err
#SBATCH --time=08:00:00

export MODULE_VERSION=5.0.1
source /opt/cluster/spack/share/spack/setup-env.sh

module purge
module load default-gcc-11.2.0
module load scorep-7.0-gcc-11.2.0-4eap4t5
module load papi-6.0.0.1-gcc-11.2.0-hnhu57s

OTF2=/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/otf2-2.3-wsvfcmjwktwfg3fgkyauk3arsqivj5b2/lib
CUB4W=/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/cubew-4.6-5cfpxt7b3lma5faqiwyep5phk75vqikg/lib/
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${OTF2}:${CUB4W}

###### set the working dir
export WORKDIR=/u/dgoz/beegfs/Pinocchio

export SCRIPT_PATH=${WORKDIR}/script
###### edit the parameters
source ${SCRIPT_PATH}/run.sh

###### switch profiling
PROFILING=0

# Tracing is meant to be run after the profiling
###### swith tracing
TRACING=1
export SCOREP_TOTAL_MEMORY=7MB # set this value according to estimated memory requirements during profiling analysis
####################################################################################################################

printf "\n\t Running on ${SLURM_NNODES} nodes \n"
printf "\n\t Running Pinocchio on hosts: ${SLURM_NODELIST}"

# get the actual topology
source ${SCRIPT_PATH}/topology.sh ${NODES}

if [ ${PROFILING} -eq 1 ]
then
    printf "\n\t Performing the profiling of Pinocchio \n\n"

    # compile the code with score-p support
    source ${SCRIPT_PATH}/compile.sh ${COMPILER}

    # perform the profiling
    source ${SCRIPT_PATH}/profiling.sh
fi

if [ ${TRACING} -eq 1 ]
then
    printf "\n\t Performing the tracing of Pinocchio \n\n"
    
    source ${SCRIPT_PATH}/tracing.sh
fi

printf "\n\t END OF GAME!!! \n\n"

exit 0
