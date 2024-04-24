# FLASH makefile definitions for x86-64 macOS
#----------------------------------------------------------------------------
# Set the AMReX library path -- manual installation for multiple variants
#----------------------------------------------------------------------------
ifeq      ($(NDIM), 1)
AMREX_PATH=/Users/adubey/FlashStuff/amrex_install/1D
else ifeq ($(NDIM), 2)
AMREX_PATH=/Users/adubey/FlashStuff/amrex_install/2D
else ifeq ($(NDIM), 3)
AMREX_PATH=/Users/adubey/FlashStuff/amrex_install/3D
else
AMREX_PATH=
endif

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------
HDF5_PATH = /opt/homebrew/Cellar/hdf5-mpi/1.14.2
#LIB_HDF5 = ${HDF5_PATH}/lib
MA28_PATH = /Users/adubey/ma28

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

LIB_NCMPI = 
MPE_PATH   =
MPI_PATH = 

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------
FCOMP   = mpif90
CCOMP   = mpicc
CPPCOMP = mpiCC
LINK    = mpif90 -std=c++11 

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

FFLAGS_OPT = -c -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -DHAVE_MPI_MODULE -fallow-argument-mismatch
FFLAGS_DEBUG = -g -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -DHAVE_MPI_MODULE -fallow-argument-mismatch
FFLAGS_TEST = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -DHAVE_MPI_MODULE -fallow-argument-mismatch
FFLAGS_AMREX = -I${AMREX_PATH}/include
FFLAGS_AMREX2D = ${FFLAGS_AMREX} -DN_DIM=2 -DNZB=1

F90FLAGS =

CFLAGS_OPT = -O2  -c
CFLAGS_DEBUG = -g -c
CFLAGS_TEST = -c

# Platform symbol
CDEFINES += -DDarwin 

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I ${HDF5_PATH}/include -DH5_USE_18_API

CFLAGS_NCMPI = -I$(LIB_NCMPI)/include

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -g -o
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
#----------------------------------------------------------------------------

LIB_OPT   = 
LIB_DEBUG =
LIB_TEST  =

#LIB_HDF4  = -lmfhdf -ldf -ljpeg -lz
#LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 /usr/lib64/libz.a
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5

LIB_PAPI  =
#LIB_MATH  = -ldfftw -ldrfftw

LIB_MPI   = 
#LIB_NCMPI = -L $(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   =

LIB_AMREX = -L${AMREX_PATH}/lib -lamrex -lstdc++
LIB_AMREX2D = ${LIB_AMREX}
LIB_STDCXX = -lstdc++ -lc++
LIB_MA28 = -L${MA28_PATH}/lib -lma28
 


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

