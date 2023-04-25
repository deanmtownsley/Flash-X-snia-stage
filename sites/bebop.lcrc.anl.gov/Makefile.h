# Flash-X makefile definitions for using the Intel compiler suite
# tested with intel/17.0.4 and intel-mpi/2017.3
#
#----------------------------------------------------------------------------
# Set the AMReX library path -- manual installation for multiple variants
#----------------------------------------------------------------------------
# AMREX_PATH   =
# MILHOJA_PATH =

#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths
#----------------------------------------------------------------------------
MPI_PATH   = ${I_MPI_ROOT}/intel64
HDF4_PATH  =
HDF5_PATH  = /lcrc/project/Flash-X/soft/intel_17.04/hdf5/1.12.2

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the parallel HDF5 wrappers which use the mpiXX compiler wrappers 
#   -- these will automatically load the proper libraries and include files.
#----------------------------------------------------------------------------
FCOMP   = ${HDF5_PATH}/bin/h5pfc
CCOMP   = ${HDF5_PATH}/bin/h5pcc
CPPCOMP = ${MPI_PATH}/bin/mpiicpc    # hdf5 CXX wrapper isn't installed
LINK    = ${HDF5_PATH}/bin/h5pfc

# pre-processor flag
PP      = -D

#----------------------------------------------------------------------------
# Compilation flags
#
#  Three sets of compilation/linking flags are defined: one for optimized
#  code, one for testing, and one for debugging.  The default is to use the
#  _OPT version.  Specifying -debug to setup will pick the _DEBUG version,
#  these should enable bounds checking.  Specifying _TEST is used for
#  flash_test, and is set for quick code generation, and (sometimes)
#  profiling.  The Makefile generated by setup will assign the generic token
#  (ex. FFLAGS) to the proper set of flags (ex. FFLAGS_OPT).
#----------------------------------------------------------------------------

OPENMP = -qopenmp

FFLAGS_OPT = -c -O2 -r8 -real-size 64 -check uninit
FFLAGS_DEBUG = -ggdb -c -O0 -r8 -real-size 64 \
	-warn all -diag-disable 10120 \
	-check output_conversion -Wunderflow \
	-ffpe-trap=invalid,zero,overflow -check bounds \
	-fimplicit-none -fp-stack-check \
	-traceback -check bounds 
FFLAGS_TEST = -O1 -c -r8 -real-size 64 \
	 -stand f18 -no-wrap-margin \
	 -fminshared -assume buffered_stdout

F90FLAGS = -DHAVE_MPI_MODULE

#The macro _FORTIFY_SOURCE adds some lightweight checks for buffer
#overflows at both compile time and run time (only active at -O1 or higher)
#http://gcc.gnu.org/ml/gcc-patches/2004-09/msg02055.html
CFLAGS_OPT   = -c -O3 -g -D_LARGEFILE64_SOURCE -D_FORTIFY_SOURCE=2 \
-diag-disable 10120
CFLAGS_DEBUG = -ggdb -c -O0 -Wundef \
	-Wconversion -Wstrict-prototypes -Wunreachable-code \
	-Winit-self -Wfloat-equal \
	-fp-stack-check
CFLAGS_TEST = -c

# CFLAGS_NCMPI   = -I$(LIB_NCMPI)/include
# CFLAGS_AMREX   =
# CFLAGS_MILHOJA =
# FFLAGS_HYPRE   = -I${HYPRE_PATH}/include
# CFLAGS_HYPRE   = -I${HYPRE_PATH}/include
# FFLAGS_AMREX   = -I${AMREX_PATH}/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -O0 -o
LFLAGS_TEST  = -o

#----------------------------------------------------------------------------
# Library specific linking
#
#  If a Flash-X source directory has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when Flash-X is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how Flash-X was setup.
#
#  Mostly handled by loading modules with Spack and h5pXX wrappers.
#----------------------------------------------------------------------------

LIB_OPT   =
LIB_DEBUG =
LIB_TEST  =

LIB_HDF5  = -L ${HDF5_PATH}/lib -lhdf5 -DH5_USE_18_API
LIB_PAPI  =
LIB_MATH  =
LIB_MPI   =
LIB_MPE   =
LIB_STDCXX = -lstdc++
LIB_LAPACK = -llapack -lblas

# LIB_HYPRE = -L$(HYPRE_PATH)/lib -lHYPRE
# LIB_AMREX = -L${AMREX_PATH}/lib -lamrex -lpthread

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when Flash-X is built.
#----------------------------------------------------------------------------

MACHOBJ =

#----------------------------------------------------------------------------
# Additional commands
#----------------------------------------------------------------------------

MV = mv -f
AR = ar -r
RM = rm -f
CD = cd
RL = ranlib
ECHO = echo

