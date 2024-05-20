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

#RUNDIR = ./

########## System-specific Configuration ############

SYSTYPE ?= base

##########################################################################
# The following variables should be passed from command line to make
# Executable name
EXEC     ?= Pinocchio
# Compilation options
COMPILER ?= gcc
DEBUG    ?= NO
OMP      ?= NO
SCOREP   ?= NO
#########################################################################


########## Preprocessor Definitions ########## 
# Displacement LPT order # 

OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT

# Dynamics of triaxial collapse #

OPTIONS += -DELL_CLASSIC
# OPTIONS += -DELL_SNG

# Building groups and fragmentation #

OPTIONS += -DCLASSIC_FRAGMENTATION

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
# OPTIONS += -DNORADIATION
# OPTIONS += -DPRINTJUNK

########## Preprocessor Definitions end ########## 

# OpenMP Configuration #

ifeq ($(OMP), YES)
# OPTIONS += -DUSE_FFT_THREADS
OMP_FLAG = -fopenmp
OMP_FFTW = -lfftw3_omp
else
OMP_FLAG =
OMP_FFTW =
endif

# Scorep configuration #
ifeq ($(SCOREP), YES)
SCOREP_FLAG = scorep --user
OPTIONS += -D_SCOREP
else
SCOREP_FLAG = 
endif

CC = $(SCOREP_FLAG) $(COMPILER)


ifeq ($(SYSTYPE), pleiadi)
CDEBUG      = -g -Wall -Wextra $(OMP_FLAG)
COPTIMIZED  = -O3 -Wall -Wextra -Wno-unused-result $(OMP_FLAG)

MPI         = /opt/cluster/spack/opt/spack/linux-centos7-broadwell/gcc-11.2.0/openmpi-4.1.3-djxjqlmzbqwq76bhh3wvgxaefnoczleg
MPI_LIBR    = -L$(MPI)/lib -lmpi
MPI_INCL    = -I$(MPI)/include

LIBDIR      = /u/dgoz/lib

FFTW        = $(LIBDIR)/fftw-3.3.10_pleiadi
FFTW_LIBR   = -L$(FFTW)/lib -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW)/include

PFFT        = $(LIBDIR)/pfft-1.0.8-alpha_pleiadi
PFFT_LIBR   = -L$(PFFT)/lib -lpfft
PFFT_INCL   = -I$(PFFT)/include

GSL         = $(LIBDIR)/gsl-2.7.1_pleiadi
GSL_LIBR    = -L$(GSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL)/include

endif


###################### LeonardoBoost ######################

ifeq ($(SYSTYPE), LeonardoBoost)

CDEBUG      = -g -O0 -Wall -Wextra $(OMP_FLAG)
#COPTIMIZED  = -O3 -Wall -Wextra -Wno-unused-result $(OMP_FLAG)
COPTIMIZED  = -O3 $(OMP_FLAG)

MPI_LIBR    = -lmpi
MPI_INCL    =

LIBDIR      = /leonardo/home/userexternal/dgoz0000/lib/nvhpc-23.11

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

###########################################################################


###################### BASE ######################

ifeq ($(SYSTYPE), base)
CDEBUG      = -ggdb3 -Wall $(OMP_FLAG)
COPTIMIZED  = -O3 -Wno-unused-result $(OMP_FLAG)
FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3 $(OMP_FFTW)
FFTW_INCL   = -I$(FFTW_INC)
PFFT_LIBR   = -L/u/mlepinzan/base/lib -lpfft
PFFT_INCL   = -I/u/mlepinzan/base/include

MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

###################### MY_MACHINE ####################

ifeq ($(SYSTYPE), my_machine)
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

INC = $(MPI_INCL) $(GSL_INCL) $(FFTW_INCL) $(PFFT_INCL)

# Library flags #

LIB = $(MPI_LIBR) $(GSL_LIBR) $(FFTW_LIBR) $(PFFT_LIBR) -lm

# Compiler options : choose from CDEBUG or COPTIMIZED #

ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# Source files #

OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o cosmo.o pinocchio.o

# Main targets and rules

$(EXEC): $(OBJECTS) Makefile
	$(CC) $(INC) $(COPTS) -o $(EXEC) $(OBJECTS) $(LIB)

run_planner: $(OBJECTS) run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o $(OBJECTS) $(LIB)	

# Generic rule for compiling source files #

%.o: %.c Makefile pinocchio.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<

memorytest: $(OBJECTS) memorytest.o Makefile
	$(CC) $(INC) $(COPTS) -o memorytest memorytest.o $(OBJECTS) $(LIB)

# Clean target #

clean:
	\rm -f *.o *~ $(EXEC) run_planner


