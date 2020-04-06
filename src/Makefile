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

EXEC = pinocchio.MG.x

# stantard options
OPTIONS += -DPLC
OPTIONS += -DNORADIATION
OPTIONS += -DTWO_LPT

# for the standard version
OPTIONS += -DTHREE_LPT
OPTIONS += -DELL_CLASSIC

# to compute collapse times numerically
#OPTIONS += -DTHREE_LPT
#OPTIONS += -DTABULATED_CT
#OPTIONS += -DELL_SNG

# for neutrino cosmology:
#OPTIONS += -DTHREE_LPT   # NB: this was not possible in Rizzo et al. (2017)
#OPTIONS += -DSCALE_DEPENDENT
#OPTIONS += -DREAD_PK_TABLE
#OPTIONS += -DONLY_MATTER_POWER
#OPTIONS += -DELL_CLASSIC

# for f(R) MG
#OPTIONS += -DSCALE_DEPENDENT
#OPTIONS += -DRECOMPUTE_DISPLACEMENTS
#OPTIONS += -DMOD_GRAV_FR
#OPTIONS += -DFR0=1.e-4
#OPTIONS += -DTABULATED_CT
#OPTIONS += -DELL_SNG

# output
#OPTIONS += -DROTATE_BOX
#OPTIONS += -DTIMELESS_SNAPSHOT

# other options
#OPTIONS += -DWHITENOISE
#OPTIONS += -DNO_RANDOM_MODULES
#OPTIONS += -DDOUBLE_PRECISION_PRODUCTS


SYSTYPE="my_machine"


ifeq ($(SYSTYPE),"my_machine")
CC          =  mpicc
CDEBUG      = -ggdb3 -Wall
COPTIMIZED  = -O3 -Wno-unused-result
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
COPTS = $(CDEBUG)
#COPTS = $(COPTIMIZED)

# source files:
OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-fftw.o GenIC.o \
	ReadParamfile.o ReadWhiteNoise.o write_fields.o allocations.o LPT.o fragment.o \
	build_groups.o distribute.o write_halos.o write_snapshot.o cosmo.o

# ifeq (SCALE_DEPENDENT_GROWTH,$(findstring SCALE_DEPENDENT_GROWTH,$(OPTIONS)))
# OBJECTS += Pk_from_CAMB.o
# endif


# rules
%.o: %.c Makefile pinocchio.h fragment.h def_splines.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<


pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) $(INC) $(COPTS) -o $(EXEC) pinocchio.o $(OBJECTS) $(LIB)

memorytest: $(OBJECTS) memorytest.o Makefile
	$(CC) $(INC) $(COPTS) -o memorytest memorytest.o $(OBJECTS) $(LIB)

clean:
	\rm -f *.o *~ $(EXEC) memorytest


