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

EXEC ?= pinocchio.x

# Compilation options #

DEBUG ?= NO
OMP   ?= NO
GPU   ?= NO

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
OPTIONS += -DCUSTOM_INTERPOLATION

####################### CPU flag #####################
ifeq ($(DEBUG), YES)
OPTIONS += -DCPU_DEBUG
endif
######################################################
####################### GPU flag #####################
ifeq ($(GPU), YES)
OPTIONS += -DGPU_OMP # Switch-on GPU using OMP
GPU_FLAG = -gpu=ccnative -target=gpu -gpu=nomanaged
ifeq ($(DEBUG), YES)
OPTIONS += -DGPU_OMP_DEBUG
GPU_FLAG = -gpu=ccnative,debug,lineinfo -target=gpu -gpu=nomanaged -mp=gpu
endif
else
GPU_FLAG =
endif
#######################################################


########## Preprocessor Definitions end ########## 

# OpenMP Configuration #

ifeq ($(OMP), YES)
# OPTIONS += -DUSE_FFT_THREADS
OMP_FLAG = -mp=multicore
OMP_FFTW = -lfftw3_omp
ifeq ($(GPU), YES)
OMP_FLAG = -mp=multicore,gpu
endif
endif

ifeq ($(OMP), NO)
ifeq ($(GPU), YES)
OMP_FLAG = -mp=gpu
endif
OMP_FLAG =
OMP_FFTW =
endif

########## System-specific Configuration ############

SYSTYPE ?= LeonardoBoost

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
CC          = nvc
# CC          = mpicc
# CDEBUG      = -ggdb3 -Wall 
# COPTIMIZED  = -acc -gpu=ccnative -Minfo=accel
COPTIMIZED  =    $(OMP_FLAG)
# COPTIMIZED  = -mp=gpu -gpu=cc80 -target=gpu -gpu=managed -Minfo=all -v  
# COPTIMIZED  = -g -mp=gpu -gpu=ccnative,debug,lineinfo -target=gpu -Minfo=all -v
FFTW_LIBR   = -L/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/fftw-3.3.10-aqomamymhcna6zjrd6nlrync7oag4fkg/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/fftw-3.3.10-aqomamymhcna6zjrd6nlrync7oag4fkg/include
PFFT_LIBR   = -L$(HOME)/lib -lpfft
PFFT_INCL   = -I$(HOME)/include

MPI_LIBR    = -L/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/openmpi-4.1.4-si334jjzldoxyl3jbbb3cf6legzshlqn/lib -lmpi
MPI_INCL    = -I/leonardo/prod/spack/03/install/0.19/linux-rhel8-icelake/gcc-11.3.0/openmpi-4.1.4-si334jjzldoxyl3jbbb3cf6legzshlqn/include
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

###################### LeonardoBoost ######################

ifeq ($(SYSTYPE), LeonardoBoost)

CC          = nvc -Minfo=all -v -Mneginfo $(GPU_FLAG)
#CDEBUG      = -g -O0 -Wall -Wextra $(OMP_FLAG)
CDEBUG      = -g -O0 $(OMP_FLAG)
#COPTIMIZED  = -O3 -Wall -Wextra -Wno-unused-result $(OMP_FLAG)
COPTIMIZED  = -O3 -fast $(OMP_FLAG)

MPI_LIBR    = -lmpi
MPI_INCL    = 

LIBDIR      =  /leonardo_work/IscrC_SCGPCT/lib/nvhpc-23.11

FFTW        = $(LIBDIR)/fftw/fftw-3.3.10
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include

GSL         = $(LIBDIR)/gsl/gsl-2.7.1
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include

PFFT        = $(LIBDIR)/pfft/pfft
PFFT_LIBR   = -L$(PFFT)/lib -lpfft
PFFT_INCL   = -I$(PFFT)/include

endif

##########################################################################

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

# Include paths for libraries #

INC =  $(PFFT_INCL) $(FFTW_INCL) $(MPI_INCL) $(GSL_INCL)

# Library flags #

LIB =  -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR)

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Source files #

OBJECTS = fmax.o  variables.o  cubic_spline_interpolation.o cosmo.o initialization.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o collapse_times.o collapse_times_GPU.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o

# Main targets and rules

pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)	


# Generic rule for compiling source files #

%.o: %.c Makefile pinocchio.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

# Clean target #

clean:
	\rm -f *.o *~ $(EXEC) run_planner


