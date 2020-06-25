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

# Compiler:
CC = /data/moretti/local/bin/mpicc

CFLAGS  = -ggdb3 -O0 -Wall -std=gnu99
#CFLAGS += -O3 -Wno-unused-result
CFLAGS += -DTWO_LPT
CFLAGS += -DNORADIATION

# Flag to choose how to compute collapse times
# Classic = LPT based
# SNG = differential equations
# MG only works with SNG!!!
CFLAGS += -DTABULATED_CT
CFLAGS += -DELL_SNG
#CFLAGS += -DELL_CLASSIC

# Flag for MG scale dependent growth
CFLAGS += -DSCALE_DEPENDENT
CFLAGS += -DMOD_GRAV_FR
CFLAGS += -DFR0=1.e-4

#CFLAGS += -DKEEP_DENSITY
#CFLAGS += -DROTATE_BOX
#CFLAGS += -DWHITENOISE
#CFLAGS += -DNO_RANDOM_MODULES
#CFLAGS += -DNORADIATION
#CFLAGS += -DSCALE_DEPENDENT_GROWTH

# Poderosa:
CFLAGS +=  -I/data/moretti/local/include
LDFLAGS = -L/data/moretti/local/lib

LDLIBS =  -lfftw3_mpi -lfftw3 -lgsl -lgslcblas -lm

all: pinocchio

clean:
	rm -f *.o pinocchio

pinocchio: fmax.o variables.o initialization.o collapse_times.o \
		   fmax-fftw.o GenIC.o ReadParamfile.o ReadWhiteNoise.o \
		   write_fields.o allocations.o LPT.o fragment.o build_groups.o \
		   distribute.o write_halos.o write_snapshot.o cosmo_new.o



fmax.o:           pinocchio.h
variables.o:	  pinocchio.h
initialization.o: pinocchio.h
collapse_times.o: pinocchio.h
fmax-fftw.o:	  pinocchio.h
GenIC.o:		  pinocchio.h
ReadParamfile.o:  pinocchio.h
ReadWhiteNoise.o: pinocchio.h
write_fields.o:   pinocchio.h
allocations.o:    pinocchio.h
LPT.o:			  pinocchio.h
fragment.o:		  pinocchio.h
build_groups.o:   pinocchio.h
distribute.o:	  pinocchio.h
write_halos.o:	  pinocchio.h
write_snapshot.o: pinocchio.h
cosmo_new.o:	  pinocchio.h

