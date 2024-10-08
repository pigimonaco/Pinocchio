module load fftw/3.3.10--openmpi--4.1.4--gcc--11.3.0
module load cuda/12.1
module load nvhpc/23.11
module load gsl/2.7.1--gcc--11.3.0


export LD_LIBRARY_PATH="/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/nvhpc-23.11-tgvw3c2exrfgdvn5qdw3rybzd3dbkkim/Linux_x86_64/23.11/math_libs/12.3/lib64:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/nvhpc-23.11-tgvw3c2exrfgdvn5qdw3rybzd3dbkkim/Linux_x86_64/23.11/comm_libs/12.3/nvshmem/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/nvhpc-23.11-tgvw3c2exrfgdvn5qdw3rybzd3dbkkim/Linux_x86_64/23.11/comm_libs/12.3/nccl/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/leonardo/pub/userexternal/glacopo0/PFFT/lib:$LD_LIBRARY_PATH"