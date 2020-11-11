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
EXEC=pinocchio.x

DEBUG=YES
DISTRIBUTION=NEW

OPTIONS += -DCLASSIC_FRAGMENTATION
#OPTIONS += -DDOUBLE_PRECISION_PRODUCTS
OPTIONS += -DLIGHT_OUTPUT

OPTIONS += -DTWO_LPT
#OPTIONS += -DTHREE_LPT
#OPTIONS += -DPLC
OPTIONS += -DNORADIATION
#OPTIONS += -DTABULATED_CT
#OPTIONS += -DELL_SNG
OPTIONS += -DELL_CLASSIC
#OPTIONS += -DSCALE_DEPENDENT
#OPTIONS += -DMOD_GRAV_FR
#OPTIONS += -DFR0=1.e-4
#OPTIONS += -DRECOMPUTE_DISPLACEMENTS

# output
#OPTIONS += -DROTATE_BOX
OPTIONS += -DTIMELESS_SNAPSHOT

# other options
#OPTIONS += -DWHITENOISE
#OPTIONS += -DNO_RANDOM_MODULES
#OPTIONS += -DDOUBLE_PRECISION_PRODUCTS


SYSTYPE="my_machine"
#SYSTYPE="Marconi100"

ifeq ($(SYSTYPE),"Marconi100")
CC          =  xlc
CDEBUG      = -g
COPTIMIZED  = -O2
FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3
FFTW_INCL   = -I$(FFTW_INC)
MPI_LIBR    = -L/cineca/prod/opt/compilers/spectrum_mpi/10.3.1/binary/lib -lmpi_ibm
MPI_INCL    = -I/cineca/prod/opt/compilers/spectrum_mpi/10.3.1/binary/include
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

ifeq ($(SYSTYPE),"amonra")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall
COPTIMIZED  = -O3 -Wno-unused-result
FFTW_LIBR   = -L$(FFTW_LIB) -lfftw3_mpi -lfftw3
FFTW_INCL   = -I$(FFTW_INC)
MPI_LIBR    = -lmpi
MPI_INCL    =
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif

ifeq ($(SYSTYPE),"my_machine")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall #-fno-stack-protector
COPTIMIZED  = -O3 -Wno-unused-result -fno-stack-protector
FFTW_LIBR   = -L$(HOME)/lib -lfftw3_mpi -lfftw3
FFTW_INCL   = -I/usr/local/include -I$(HOME)/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -lgsl -lgslcblas -lm
GSL_INCL    = -I/usr/include
endif


# included libraries, shared libraries 
INC =  $(FFTW_INCL)  $(MPI_INCL) $(GSL_INCL)
LIB =  $(FFTW_LIBR)  $(MPI_LIBR) $(GSL_LIBR)


# compiler options : choose from CDEBUG or COPTIMIZED
ifeq ($(DEBUG), YES)
COPTS = $(CDEBUG)
else
COPTS = $(COPTIMIZED)
endif

# source files:
OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-fftw.o GenIC.o \
	ReadParamfile.o ReadWhiteNoise.o write_fields.o allocations.o LPT.o fragment.o \
	build_groups.o write_halos.o write_snapshot.o cosmo.o 

ifeq ($(DISTRIBUTION), OLD)
OBJECTS += distribute.o
else
OBJECTS += distribute_new.o
endif

# rules
%.o: %.c Makefile pinocchio.h fragment.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<


pinocchio: $(OBJECTS) pinocchio.o Makefile 
	$(CC) $(INC) $(COPTS) -o $(RUNDIR)$(EXEC) pinocchio.o $(OBJECTS) $(LIB)

new: $(OBJECTS) pinocchio.o Makefile 
	$(CC) $(INC) $(COPTS) -o pinocchio_new.x pinocchio.o $(OBJECTS) $(LIB)

run_planner: allocations.o initialization.o ReadParamfile.o cosmo.o variables.o run_planner.o Makefile
	$(CC) $(INC) $(COPTS) -o run_planner run_planner.o allocations.o initialization.o ReadParamfile.o cosmo.o variables.o $(LIB)

clean:
	\rm -f *.o *~ $(EXEC) run_planner


