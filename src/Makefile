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


# ============================
# build options
# ============================
#

#                                          the basename for the executable
EXEC_BASE   = pinocchio

#                                          compile with -g and sanitizers or with machine-specific opts
#BUILDTYPE   = DEBUG                     
BUILDTYPE   = PRODUCTION

#USE_MPIP    = ON

#                                          activate openmp regions in the code
#USE_OPENMP  = ON

#                                          use fftw threads pthreads
#FFT_THREADS = PTHREADS
#FFT_THREADS = OPENMP

#                                          determines the machine on which to compile
#SYSTYPE     = my_machine
#SYSTYPE     = exabed
#SYSTYPE     = hotcat
SYSTYPE     = CASS

#                                          activate a tracer
#PROFILER    = gperftools                


# ============================
# code options
# ============================
#
OPTIONS += -DASCII_WRITE
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT
#OPTIONS += -DPLC
#OPTIONS += -DROTATE_BOX
#OPTIONS += -DWHITENOISE
#OPTIONS += -DNO_RANDOM_MODULES
#OPTIONS += -DNORADIATION
#OPTIONS += -DSCALE_DEPENDENT_GROWTH

OPTIONS += -DTEST_ONLY


# ===========================
# set variables, flags and libraries
# ===========================
#

#                                             tracer
#
ifeq ($(PROFILER), gperftools)
OPTIONS += -DUSE_GPERFTOOLS
endif


#                                             OpenMP
#
ifeq ($(USE_OPENMP), ON)
EXEC_SUFFIX_OMP  =.omp
OMP_FLAG         = -fopenmp
endif

#                                             MpiP
#
ifeq ($(USE_MPIP), ON)
EXEC_SUFFIX_MPIP  =.mpiP
EXTRA_LIBS       += -L/scratch/Software/lib -lmpiP -lbfd -liberty -lunwind
MPIP_CFLAG        = -g
endif




#                                             FFT threads
#
ifeq ($(FFT_THREADS),PTHREADS)
FFTW_THREAD_LIB  = -lfftw3_threads
THREAD_LIB       = -lm -lpthread
FPIC             = -fPIC
EXEC_SUFFIX_FFT  =.fftthreads
OPTIONS         += -DUSE_FFT_THREADS
else
ifeq  ($(FFT_THREADS),OPENMP)
FFTW_THREAD_LIB  = -lfftw3_omp
THREAD_LIB       = -lm
FPIC             = -fPIC
ifneq ($(USE_OPENMP), ON)
OMP_FLAG         = -fopenmp
endif
EXEC_SUFFIX_FFT   =.fftomp
OPTIONS         += -DUSE_FFT_THREADS
endif
endif


#                                             build type
#
ifeq ($(BUILDTYPE), DEBUG)
EXEC_SUFFIX_T =.d
else
EXEC_SUFFIX_T =.x
endif

EXEC = $(EXEC_BASE)$(EXEC_SUFFIX_MPIP)$(EXEC_SUFFIX_OMP)$(EXEC_SUFFIX_FFT)$(EXEC_SUFFIX_T)

# ============================
# compile options
# ============================
#

x86_DBG = -ggdb3 -fPIC -march=native -Wall -fno-omit-frame-pointer -fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero
ARM_DBG = -ggdb3 -fPIC -march=native -Wall -fno-omit-frame-pointer -fsanitize=float-divide-by-zero -fsanitize=bounds -fsanitize=signed-integer-overflow -fsanitize=vla-bound -fsanitize=undefined -fsanitize=integer-divide-by-zero

x86_OPT = -O3 -Wno-unused-result -march=native -fomit-frame-pointer -fstrict-aliasing ${MPIP_CFLAG}
ARM_OPT = -O3 -Wno-unused-result -march=native -fomit-frame-pointer -fstrict-aliasing ${MPIP_CFLAG}

#CPROFILE    = -fprofile-generate -ftest-coverage #--coverage -fPIC

# [ my_machine ] ==========================  
ifeq ($(SYSTYPE),my_machine)
CC          =  mpicc

ifeq ($(BUILDTYPE), DEBUG)
COPTS         = $(OMP_FLAG) $(x86_DBG) $(FPIC)
else
COPTS         = $(OMP_FLAG) $(x86_OPT) -fPIC $(FPIC) #-fopt-info-missed
endif

FFTW_LIBR   = -L/usr/lib  -lfftw3_mpi $(FFTW_THREAD_LIB) -lfftw3 $(THREAD_LIB)
FFTW_INCL   = -I/usr/include
PFFT_LIBR   = -L/usr/lib -lpfft
PFFT_INCL   = -I/usr/include
MPI_LIBR    = -lmpi
MPI_INCL    =
GSL_LIBR    = -lgsl -lgslcblas -lm
GSL_INCL    = -I/usr/include
endif


# [ exabed ] ==========================  
ifeq ($(SYSTYPE),exabed)
CC          =  mpicc

ifeq ($(BUILDTYPE), DEBUG)
COPTS         = $(OMP_FLAG) $(ARM_DBG) $(FPIC)
else
COPTS         = $(OMP_FLAG) $(ARM_OPT) -fPIC $(FPIC) #-fopt-info-missed
endif

FFTW_LIBR   = -L$(EBROOTFFTW)/lib -lfftw3_mpi $(FFTW_THREAD_LIB) -lfftw3 $(THREAD_LIB)
FFTW_INCL   = -I$(EBROOTFFTW)/include
PFFT_LIBR   = -L/home/inaf/lib -lnew_pfft
PFFT_INCL   = -I/home/inaf/include
MPI_LIBR    =
MPI_INCL    =
GSL_LIBR    = -L$(EBROOTGSL)/lib -lgsl -lgslcblas -lm
GSL_INCL    = -I$(EBROOTGSL)/include
endif


# [ hotcat ] ==========================  
ifeq ($(SYSTYPE),hotcat)
CC          =  mpicc

ifeq ($(BUILDTYPE), DEBUG)
COPTS         = $(OMP_FLAG) $(x86_DBG) $(FPIC)
else
COPTS         = $(OMP_FLAG) $(x86_OPT) -fPIC $(FPIC) #-fopt-info-missed
endif


FFTW_LIBR   = -L$(FFTW_LIB)  -lfftw3_mpi $(FFTW_THREAD_LIB) -lfftw3 $(THREAD_LIB)
FFTW_INCL   = -I$(FFTW_INC)
PFFT_LIBR   = -L$(PFFT_LIB) -lpfft
PFFT_INCL   = -I$(PFFT_INC)
MPI_LIBR    = -lmpi
MPI_INCL    =
GSL_LIBR    = -L$(GSL_LIB) -lgsl -lgslcblas -lm
GSL_INCL    = -I$(GSL_INC)
endif


# [ CASS ] ==========================  
ifeq ($(SYSTYPE),CASS)
CC          =  mpicc

ifeq ($(BUILDTYPE), DEBUG)
COPTS         = $(OMP_FLAG) $(ARM_DBG) $(FPIC)
else
COPTS         = $(OMP_FLAG) $(ARM_OPT) -fPIC $(FPIC) #-fopt-info-missed
endif

FFTW_LIBR   = -L/u/luca/lib -lfftw3_mpi $(FFTW_THREAD_LIB) -lfftw3
FFTW_INCL   = -I/u/luca/include
PFFT_LIBR   = -L/u/luca -lpfft
PFFT_INCL   = -I/u/luca
MPI_LIBR    = -lmpi
MPI_INCL    =
GSL_LIBR    = -L/usr/lib -lgsl -lgslcblas
GSL_INCL    = -I/usr/include

endif


# included libraries, shared libraries 
INC =  $(FFTW_INCL) $(PFFT_INCL) $(MPI_INCL) $(GSL_INCL)
LIB =  $(EXTRA_LIBS) -lm $(PFFT_LIBR) $(FFTW_LIBR) $(MPI_LIBR) $(GSL_LIBR) 

ifdef PROFILER
ifeq ($(PROFILER), "gperftools")
LIB += -lprofiler
endif
endif



# ==========================
# source files:
OBJECTS = fmax.o variables.o initialization.o collapse_times.o fmax-pfft.o GenIC.o GenIC_large.o \
	ReadParamfile.o ReadWhiteNoise.o write_fields.o allocations.o LPT.o fragment.o \
	build_groups.o distribute.o write_halos.o write_snapshot.o cosmo.o

ifeq (SCALE_DEPENDENT_GROWTH,$(findstring SCALE_DEPENDENT_GROWTH,$(OPTIONS)))
OBJECTS += Pk_from_CAMB.o
endif


# rules
%.o: %.c Makefile pinocchio.h fragment.h
	$(CC) $(INC) $(COPTS) $(OPTIONS) -c $<


pinocchio: $(OBJECTS) pinocchio.o Makefile
	$(CC) -fuse-ld=gold $(INC) $(COPTS) -o $(EXEC) pinocchio.o $(OBJECTS) $(LIB)

memorytest: $(OBJECTS) memorytest.o Makefile
	$(CC) $(INC) $(COPTS) -o memorytest memorytest.o $(OBJECTS) $(LIB)

clean:
	rm -f *.o *~ $(EXEC) memorytest


