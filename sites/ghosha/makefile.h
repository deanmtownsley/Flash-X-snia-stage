# FLASH makefile definitions for x86-64 macOS
#----------------------------------------------------------------------------
# Set the AMReX library path -- manual installation for multiple variants
#----------------------------------------------------------------------------
ifeq      ($(NDIM), 1)
AMREX_PATH=/Users/dubey/amrex_install/1D
else ifeq ($(NDIM), 2)
AMREX_PATH=/Users/dubey/amrex_install/2D
else ifeq ($(NDIM), 3)
AMREX_PATH=/Users/dubey/amrex_install/3D
else
AMREX_PATH=
endif

#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- managed by loading with Spack 
#----------------------------------------------------------------------------
HDF5_PATH = /usr/local
HYPRE_PATH = 
ZLIB_PATH  =
PAPI_PATH  =
PAPI_FLAGS =
LIB_NCMPI = /usr/local

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the parallel HDF5 wrappers which use the mpiXX compiler wrappers 
#   -- these will automatically load the proper libraries and include files.
#----------------------------------------------------------------------------
FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC
LINK    = mpif90  -std=c++11

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

FFLAGS_OPT = -c -O2 -fdefault-real-8 -fdefault-double-8 -Wuninitialized -DHAVE_MPI_MODULE -fallow-argument-mismatch
FFLAGS_DEBUG = -ggdb -c  -fdefault-real-8 -fdefault-double-8 \
	-Wall -Wextra -Waliasing \
	-Wsurprising -Wconversion -Wunderflow \
	-ffpe-trap=invalid,zero,overflow -fbounds-check \
	-fimplicit-none -fstack-protector-all \
	-fbacktrace -fbounds-check \
	-ffree-line-length-none -DHAVE_MPI_MODULE -fallow-argument-mismatch
FFLAGS_TEST = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
	-ffree-line-length-none -DHAVE_MPI_MODULE -fallow-argument-mismatch

FFLAGS_HYPRE =
FFLAGS_AMREX = -I${AMREX_PATH}/include
CFLAGS_AMREX = -I${AMREX_PATH}/include

F90FLAGS =

#The macro _FORTIFY_SOURCE adds some lightweight checks for buffer
#overflows at both compile time and run time (only active at -O1 or higher)
#http://gcc.gnu.org/ml/gcc-patches/2004-09/msg02055.html
CFLAGS_OPT = -c -O2 -Wuninitialized -D_FORTIFY_SOURCE=2 
CFLAGS_DEBUG = -ggdb -c -O0 -Wno-div-by-zero -Wundef \
	-Wconversion -Wstrict-prototypes -Wunreachable-code \
	-pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
	-Wunsafe-loop-optimizations -Wpadded -fstack-protector-all
CFLAGS_TEST = -c

# Platform symbol
CDEFINES += -DDarwin

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I ${HDF5_PATH}/include -DH5_USE_18_API
CFLAGS_NCMPI = -I ${NCMPI_PATH}/include

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
#  If a FLASH module has a 'LIBRARY xxx' line in its Config file, we need to
#  create a macro in this Makefile.h for LIB_xxx, which will be added to the
#  link line when FLASH is built.  This allows us to switch between different
#  (incompatible) libraries.  We also create a _OPT, _DEBUG, and _TEST
#  library macro to add any performance-minded libraries (like fast math),
#  depending on how FLASH was setup.
#
#  Mostly handled by loading modules with Spack and h5pXX wrappers.
#----------------------------------------------------------------------------
 
LIB_OPT   =
LIB_DEBUG =
LIB_TEST  =
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 -lz
LIB_PAPI  =
LIB_MATH  =
LIB_MPI   = 
LIB_MPE   =
LIB_HYPRE =
LIB_AMREX = -L${AMREX_PATH}/lib -lamrex
LIB_STDCXX = -lstdc++

# Uncomment the following line to use electic fence memory debugger.
# Need the following environmental variable (see env.sh):
# export EF_ALLOW_MALLOC_0=1
#CONFIG_LIB = -L/usr/lib64 -lefence

#----------------------------------------------------------------------------
# Additional machine-dependent object files
#
#  Add any machine specific files here -- they will be compiled and linked
#  when FLASH is built.
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

