#  *****************************************************************
#  *                        PINOCCHIO  V4.1                        *
#  *  (PINpointing Orbit-Crossing Collapsed HIerarchical Objects)  *
#  *****************************************************************
#  
#  This code was written by
#  Pierluigi Monaco
#  Copyright (C) 2016
#  
#  web page: http://adlibitum.oats.inaf.it/monaco/pinocchio.html
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#


# Directory settings #

RUNDIR = ./

# Executable name #

EXEC = pinocchio.x

# Compilation options #

DEBUG = YES
OMP   = YES

all: pinocchio_cufftmp

########## Preprocessor Definitions ########## 

# Displacement LPT order # 

OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# Dynamics of triaxial collapse #

OPTIONS += -DELL_CLASSIC
# OPTIONS += -DELL_SNG

# Building groups and fragmentation #

# OPTIONS += -DCLASSIC_FRAGMENTATION

# Interpolation of collapse time from a TABLE #

# OPTIONS += -DTABULATED_CT
# OPTIONS += -DTRILINEAR
# OPTIONS += -DBILINEAR_SPLINE
# OPTIONS += -DALL_SPLINE
# OPTIONS += -DASCII
# OPTIONS += -DHISTO

# Output #

# OPTIONS += -DPLC
# OPTIONS += -DSNAPSHOT
# OPTIONS += -DADD_RMAX_TO_SNAPSHOT
# OPTIONS += -DLIGHT_OUTPUT

# Modified Gravity models #

# OPTIONS += -DSCALE_DEPENDENT
# OPTIONS += -DMOD_GRAV_FR
# OPTIONS += -DFR0=1.e-4

# Other options # 

# OPTIONS += -DRECOMPUTE_DISPLACEMENTS
# OPTIONS += -DUSE_GPERFTOOLS
# OPTIONS += -DLONGIDS
# OPTIONS += -DWHITENOISE
# OPTIONS += -DNO_RANDOM_MODULES
# OPTIONS += -DDOUBLE_PRECISION_PRODUCTS
# OPTIONS += -DVETTORIALIZZA
OPTIONS += -DNORADIATION
# OPTIONS += -DPRINTJUNK

# OPTIONS += -DGPU_INTERPOLATION
# OPTIONS += -DCUSTOM_INTERPOLATION

# offloading options
OPTIONS += -DCUFFTMP

########## Preprocessor Definitions end ########## 

# OpenMP Configuration #

ifeq ($(OMP), YES)
OPTIONS += -DUSE_FFT_THREADS
OMP_FLAG = -fopenmp
OMP_FFTW = -lfftw3_omp
else
OMP_FLAG =
OMP_FFTW =
endif

########## System-specific Configuration ############

SYSTYPE="LEONARDO"

###################### MARCONI ######################

ifeq ($(SYSTYPE),"Marconi100")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O2 -Wno-unused-result  $(OMP_FLAG)
FFTW_LIBR   = -L$(HOME)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(HOME)/include
PFFT_LIBR   = -lpfft
PFFT_INCL   = 
MPI_LIBR    = -L$(SMPI_ROOT)/lib -lmpi_ibm
MPI_INCL    = -I$(SMPI_ROOT)/include
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

###################### LEONARDO ######################

ifeq ($(SYSTYPE),"LEONARDO")
CC          =  mpicc
CXX         =  mpic++
CDEBUG      = -g -O0 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
#COPTIMIZED  = -mp=gpu -gpu=cc80  -gpu=managed -Minfo=mp $(OMP_FLAG)                                                                                                                                                                          
# COPTIMIZED  = -O3 -Wno-unused-result -foffload=-lm -fno-fast-math -fno-associative-math $(OMP_FLAG)                                                                                                                                         
FFTW_LIBR   = -L/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/fftw-3.3.10-aqomamymhcna6zjrd6nlrync7oag4fkg/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/fftw-3.3.10-aqomamymhcna6zjrd6nlrync7oag4fkg/include
PUBLIC      = /leonardo/pub/userexternal/glacopo0/
PFFT_LIBR   = -L$(PUBLIC)/PFFT/lib -lpfft
PFFT_INCL   = -I$(PUBLIC)/PFFT/include
CUDA_LIB    = -L/leonardo/prod/opt/compilers/cuda/12.1/none/lib64 -lcuda -lcudart -lnvidia-ml
CUDA_INCL   = -I/leonardo/prod/opt/compilers/cuda/12.1/none/include
# NVSHMEM_INC = -I/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-23.11-wnrvac5a7nx7cnbtaamktaaqljmgorva/Linux_x86_64/23.11/comm_libs/nvshmem/include                                                               
# NVSHMEM_LIB = -L/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-23.11-wnrvac5a7nx7cnbtaamktaaqljmgorva/Linux_x86_64/23.11/comm_libs/nvshmem/lib -lnvshmem_host -lnvshmem_device                                   
CUFFT_LIB   = -L/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-23.11-wnrvac5a7nx7cnbtaamktaaqljmgorva/Linux_x86_64/23.11/math_libs/lib64 -lcufftMp
CUFFT_INCL  = -I/leonardo/prod/spack/5.2/install/0.21/linux-rhel8-icelake/gcc-8.5.0/nvhpc-23.11-wnrvac5a7nx7cnbtaamktaaqljmgorva/Linux_x86_64/23.11/math_libs/include/cufftmp
# MPI_LIBR    = -L/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/openmpi-4.1.4-si334jjzldoxyl3jbbb3cf6legzshlqn/lib -lmpi                                                                                                
# MPI_INCL    = -I/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/openmpi-4.1.4-si334jjzldoxyl3jbbb3cf6legzshlqn/include                                                                                                  
GSL_LIBR    = -L/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/gsl-2.7.1-wsvpkwadh36phlaba3key5sd4cemsxds/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/gsl-2.7.1-wsvpkwadh36phlaba3key5sd4cemsxds/include
endif


ifeq (CUFFTMP,$(findstring CUFFTMP,$(OPTIONS)))
# included libraries, shared libraries                                                                                                                                                                                                        
INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL) $(CUDA_INCL) $(NVSHMEM_INC) $(CUFFT_INCL)
LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(CUDA_LIB) $(NVSHMEM_LIB) $(CUFFT_LIB)
else
# included libraries, shared libraries


INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)
LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR)
endif


###################### PLEIADI ######################

ifeq ($(SYSTYPE),"pleiadi")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/fftw-3.3.10-bap3ootbcmiypnmo7rouk57jbdhfjuty/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/fftw-3.3.10-bap3ootbcmiypnmo7rouk57jbdhfjuty/include
PFFT_LIBR   = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/pfft-1.0.8-alpha-qy6lqerctfrhjzaoamzjhbpga2ofejqf/lib -lpfft
PFFT_INCL   = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/pfft-1.0.8-alpha-qy6lqerctfrhjzaoamzjhbpga2ofejqf/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/gsl-2.7.1-uehbowahpvubxundtx3cz4a4a5sdwzun/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I/opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/gsl-2.7.1-uehbowahpvubxundtx3cz4a4a5sdwzun/include
endif

###################### AMONRA ######################

ifeq ($(SYSTYPE),"amonra")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW_INC)
PFFT_LIBR   = -L/u/pmonaco/lib -lpfft
PFFT_INCL   = -I/u/monaco/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

###################### MY_MACHINE ####################

ifeq ($(SYSTYPE),"my_machine")
CC          =  mpicc 
CDEBUG      =  -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  =  -g -O3 -Wno-unused-result -march=native $(OMP_FLAG)
FFTW_LIBR   = -L/usr/local/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/usr/local/include 
PFFT_LIBR   = -L/usr/local/lib -lpfft
PFFT_INCL   = 
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L/usr/local/lib -lgsl -lgslcblas
GSL_INCL    = -I/usr/local/include
endif

########## System-specific Configuration end ########## 

ifeq (CUFFTMP,$(findstring CUFFTMP,$(OPTIONS)))
# included libraries, shared libraries 
INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL) $(CUDA_INCL) $(NVSHMEM_INC) $(CUFFT_INCL)
LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) $(CUDA_LIB) $(NVSHMEM_LIB) $(CUFFT_LIB)
else
# included libraries, shared libraries 
INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL) 
LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) 
endif


# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Source files #

OBJECTS = fmax.o  variables.o cosmo.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o

ifeq (CUFFTMP,$(findstring CUFFTMP,$(OPTIONS)))
NVCC               = nvcc
#OBJECTS            += cuda_functions.o                                                                                                                                                                                                       
COMPTIMIZED_NVIDIA = -O0 -g -G -std=c++17 --generate-code arch=compute_80,code=sm_80
cuda_functions.o: cuda_functions.cu pinocchio.h
	$(NVCC) $(INC) $(COMPTIMIZED_NVIDIA) $(OPTIONS) -c $< $(LIB)
endif

# rules                                                                                                                                                                                                                                       
%.o: %.c Makefile pinocchio.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

ifeq (CUFFTMP,$(findstring CUFFTMP,$(OPTIONS)))
pinocchio_cufftmp: $(OBJECTS) pinocchio.o cuda_functions.o Makefile
	$(CXX) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o cuda_functions.o $(OBJECTS) $(LIB)
else
pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)
endif

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)

clean:
	\rm -f *.o *~ $(EXEC) run_planner
