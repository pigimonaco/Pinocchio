# Build configuration
#
# Possibly add include directories and library options by calling like:
#   make INC="..." LIB="..."
#
# Override the compiler choice (mpicc) by passing CC="...".


###########
# Options #
###########

# Activate, or not, compile-time options modifying program behavior.
OPTIONS += -DTWO_LPT
OPTIONS += -DTHREE_LPT
OPTIONS += -DELL_CLASSIC
# OPTIONS += -DPLC
# OPTIONS += -DCLASSIC_FRAGMENTATION
# OPTIONS += -DLIGHT_OUTPUT
# OPTIONS += -DSNAPSHOT
OPTIONS += -DNORADIATION
# OPTIONS += -DTABULATED_CT
# OPTIONS += -DELL_SNG
# OPTIONS += -DSCALE_DEPENDENT
# OPTIONS += -DMOD_GRAV_FR
# OPTIONS += -DFR0=1.e-4
# OPTIONS += -DRECOMPUTE_DISPLACEMENTS
# OPTIONS += -DADD_RMAX_TO_SNAPSHOT
# OPTIONS += -DLONGIDS
# OPTIONS += -DLIGHT_OUTPUT
# OPTIONS += -DWHITENOISE
# OPTIONS += -DNO_RANDOM_MODULES
# OPTIONS += -DDOUBLE_PRECISION_PRODUCTS
# OPTIONS += -DUSE_FFT_THREADS
OPTIONS += -DCUSTOM_INTERPOLATION


# Use C compiler with MPI support by default.
CC = mpicc

# Define compiler flags for three different builds: release, debug, test.
CFLAGS_RELEASE = -O2 -Wno-unused-result -Wno-format-overflow
CFLAGS_DEBUG   = -ggdb3 -Wall -fno-omit-frame-pointer
CFLAGS_TEST    = -ggdb3 -Wall -fno-omit-frame-pointer -fsanitize=address

# Specify the libraries we depend on.
# The "override" keyword means we append to whatever value for LIB
# the user may have passed in or defined in the shell environment.
override LIB += -lpfft -lfftw3_mpi -lfftw3 -lgslcblas -lgsl -lm

# Optionally support OpenMP multithreading.
ifeq (USE_FFT_THREADS,$(findstring USE_FFT_THREADS,$(OPTIONS)))
override INC += -fopenmp
override LIB += -lfftw3_omp
endif

# Name all required object files except those of the executables.
OBJECTS = fmax.o variables.o cubic_spline_interpolation.o cosmo.o initialization.o fmax-pfft.o GenIC.o \
	ReadParamfile.o allocations.o LPT.o distribute.o collapse_times.o \
	fragment.o build_groups.o write_halos.o write_snapshot.o
# collapse_times_GPU.o
ifeq (WHITENOISE,$(findstring WHITENOISE,$(OPTIONS)))
OBJECTS += ReadWhiteNoise.o
endif

# We will build three different targets: release, debug, test.
OBJECTS_RELEASE = $(addprefix build/release/objects/, $(OBJECTS))
OBJECTS_DEBUG   = $(addprefix build/debug/objects/,   $(OBJECTS))
OBJECTS_TEST    = $(addprefix build/test/objects/,    $(OBJECTS))

# Define install target using variable names as per the GNU coding standards.
# Installing Pinocchio is not necessary, we can just call the executable
# from the "build" directory. Installation typically only makes sense when
# building a lightweight container that only contains the binaries.
DESTDIR =
PREFIX  = /usr
TARGET  = $(DESTDIR)$(PREFIX)


#########
# Rules #
#########

.PHONY: all release debug test install inspect clean

all: release debug test

release: echo_release \
         build/release \
         build/release/objects \
         $(OBJECTS_RELEASE) \
         build/release/objects/pinocchio.o \
         build/release/objects/run_planner.o \
         build/release/pinocchio \
         build/release/run_planner

debug: echo_debug \
       build/debug \
       build/debug/objects \
       $(OBJECTS_DEBUG) \
       build/debug/objects/pinocchio.o \
       build/debug/objects/run_planner.o \
       build/debug/pinocchio \
       build/debug/run_planner

test: echo_test \
      build/test \
      build/test/objects \
      $(OBJECTS_TEST) \
      build/test/objects/pinocchio.o \
      build/test/objects/run_planner.o \
      build/test/pinocchio \
      build/test/run_planner

install:
	install -d $(TARGET)/bin/
	install -m 777 build/release/pinocchio $(TARGET)/bin/
	install -m 777 build/release/run_planner $(TARGET)/bin/pinocchio_planner

inspect:
	@echo "compile options: $(OPTIONS)"
	@echo "include options: $(INC)"
	@echo "library options: $(LIB)"
	@echo "install target:  $(TARGET)"

clean:
	rm -rf build/


.PHONY: echo_release echo_debug echo_test

echo_release:
	@echo "Building \"release\" version optimized for performance."

echo_debug:
	@echo "Building \"debug\" version for interactive debugging."

echo_test:
	@echo "Building \"test\" version for memory-leak testing."


build:
	mkdir -p $@


build/release build/release/objects: build
	mkdir -p $@

build/release/objects/%.o: src/%.c src/pinocchio.h src/def_splines.h Makefile
	$(CC) -c $(CFLAGS_RELEASE) $(INC) $(OPTIONS) -o $@ $<

build/release/pinocchio: build/release/objects/pinocchio.o $(OBJECTS_RELEASE) Makefile
	$(CC) $(CFLAGS_RELEASE) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_RELEASE) $(LIB)

build/release/run_planner: build/release/objects/run_planner.o $(OBJECTS_RELEASE) Makefile
	$(CC) $(CFLAGS_RELEASE) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_RELEASE) $(LIB)


build/debug build/debug/objects: build
	mkdir -p $@

build/debug/objects/%.o: src/%.c src/pinocchio.h src/def_splines.h Makefile
	$(CC) -c $(CFLAGS_DEBUG) $(INC) $(OPTIONS) -o $@ $<

build/debug/pinocchio: build/debug/objects/pinocchio.o $(OBJECTS_DEBUG) Makefile
	$(CC) $(CFLAGS_DEBUG) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_DEBUG) $(LIB)

build/debug/run_planner: build/debug/objects/run_planner.o $(OBJECTS_DEBUG) Makefile
	$(CC) $(CFLAGS_DEBUG) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_DEBUG) $(LIB)


build/test build/test/objects: build
	mkdir -p $@

build/test/objects/%.o: src/%.c src/pinocchio.h src/def_splines.h Makefile
	$(CC) -c $(CFLAGS_TEST) $(INC) $(OPTIONS) -o $@ $<

build/test/pinocchio: build/test/objects/pinocchio.o $(OBJECTS_TEST) Makefile
	$(CC) $(CFLAGS_TEST) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_TEST) $(LIB)

build/test/run_planner: build/test/objects/run_planner.o $(OBJECTS_TEST) Makefile
	$(CC) $(CFLAGS_TEST) -o $@ $(@D)/objects/$(@F).o $(OBJECTS_TEST) $(LIB)
