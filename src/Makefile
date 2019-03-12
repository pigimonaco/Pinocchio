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


RUNDIR = ./

EXEC     = pinocchio.x
OPTIONS += -DASCII_WRITE
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT
OPTIONS += -DPLC
#OPTIONS += -DROTATE_BOX
#OPTIONS += -DWHITENOISE
#OPTIONS += -DNO_RANDOM_MODULES
#OPTIONS += -DNORADIATION
#OPTIONS += -DSCALE_DEPENDENT_GROWTH


#PROFILER="gperftools"

ifdef PROFILER
ifeq ($(PROFILER), "gperftools")
OPTIONS += -DUSE_GPERFTOOLS
endif
endif

BUILDTYPE="DEBUG"
#BUILDTYPE="PRODUCTION""

SYSTYPE="my_machine"


ifeq ($(SYSTYPE),"my_machine")
CC          =  mpicc

ifeq ($(BUILDTYPE), "DEBUG")
COPTS       = -ggdb3 -fPIC -march=native -Wall -fno-omit-frame-pointer -fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero 
EXEC_SUFFIX =.d
else
COPTS  = -O3 -Wno-unused-result -march=native -fopt-info-missed -fPIC
EXEC_SUFFIX =.x
#CPROFILE    = -fprofile-generate -ftest-coverage #--coverage -fPIC
endif

FFTW_LIBR   = -L/usr/lib -lfftw3_mpi -lfftw3
FFTW_INCL   = -I/usr/local/include
PFFT_LIBR   = -L/usr/lib -lpfft
PFFT_INCL   = -I/usr/include
MPI_LIBR    = 
MPI_INCL    =
GSL_LIBR    = -lgsl -lgslcblas -lm
GSL_INCL    = -I/usr/include
endif


# included libraries, shared libraries 
INC =  $(FFTW_INCL)  $(PFFT_INCL) $(GSL_INCL)
LIB =  $(FFTW_LIBR)  $(PFFT_LIBR)  $(GSL_LIBR)

ifdef PROFILER
ifeq ($(PROFILER), "gperftools")
LIB += -lprofiler
endif
endif

EXEC = $(EXEC_BASE)$(EXEC_SUFFIX)

# source files:
OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o \
	ReadParamfile.o ReadWhiteNoise.o write_fields.o allocations.o LPT.o fragment.o \
	build_groups.o distribute.o write_halos.o write_snapshot.o cosmo.o

ifeq (SCALE_DEPENDENT_GROWTH,$(findstring SCALE_DEPENDENT_GROWTH,$(OPTIONS)))
OBJECTS += Pk_from_CAMB.o
endif


# rules
%.o: %.c Makefile pinocchio.h fragment.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<


pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(EXEC) pinocchio.o $(OBJECTS) $(LIB)

memorytest: $(OBJECTS) memorytest.o Makefile
	$(CC) $(INC) $(COPTS) -o memorytest memorytest.o $(OBJECTS) $(LIB)

clean:
	\rm -f *.o *~ $(EXEC) memorytest


