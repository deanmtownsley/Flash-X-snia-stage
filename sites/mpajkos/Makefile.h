#Works for Mike's local laptop
# FLASH makefile definitions for ix86-64 Linux (gfortran compiler)

#----------------------------------------------------------------------------
# Set the HDF/HDF5 library paths -- these need to be updated for your system
#----------------------------------------------------------------------------


ifeq      ($(NDIM), 1)
  AMREX_PATH=/usr/local/Cellar/amrex-FlashFluxRegister/install_gnu_1D
else ifeq ($(NDIM), 2)
  AMREX_PATH=/usr/local/Cellar/amrex-FlashFluxRegister/install_gnu_2D
  //^^built w/ just FC=gfortran in mpich/.../mpif90 call
else ifeq ($(NDIM), 3)
  AMREX_PATH=/usr/local/Cellar/amrex-FlashFluxRegister/install_gnu_3D
endif

HDF5_PATH = /usr/local/Cellar/hdf5/1.10.5_1
#HDF5_PATH = /Users/mpajkos/hdf5/hdf5-1.8.18-gnu
HYPRE_PATH = /usr/local

ZLIB_PATH  =

PAPI_PATH  =
PAPI_FLAGS =

LIB_NCMPI =
MPE_PATH  =
MPI_PATH  = /usr/local/Cellar/mpich/3.3.1/bin
#MPI_PATH  = /usr/local/bin

#MESA_DIR  = /Volumes/fs0/smc/software/mesa_5527
#include $(MESA_DIR)/utils/makefile_header
#LOAD_OTHER = -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------
FCOMP   = $(MPI_PATH)/mpif90
CCOMP   = $(MPI_PATH)/mpicc
CPPCOMP = $(MPI_PATH)/mpiCC
LINK    = $(MPI_PATH)/mpif90

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
#-ffpe-trap=invalid,zero,overflow

OPENMP = -fopenmp

FFLAGS_OPT = -c -g -O3 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none \
-Wuninitialized -march=native

#FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
#-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
#-Wsurprising -Wconversion -Wunderflow \
#-fcheck=all \
#-fbacktrace -fdump-core -finit-real=snan \
#-finit-integer=-999999 -fimplicit-none -fstack-protector-all -fopenmp
FFLAGS_DEBUG = -ggdb -c -O0 -fdefault-real-8 -fdefault-double-8 \
-pedantic -Wall -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
  -ffpe-trap=invalid,zero,overflow -fbounds-check \
-fimplicit-none -fstack-protector-all \
-ffree-line-length-none

FFLAGS_HYPRE = -I${HYPRE_PATH}/include

FFLAGS_TEST = -c -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none

F90FLAGS = -I ${HDF5_PATH}/include -DH5_USE_18_API  #-I$(MESA_DIR)/include

CFLAGS_OPT = -O3  -c -DDarwin
CFLAGS_DEBUG = -ggdb -c -Wno-div-by-zero -Wundef  \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Wextra -Winit-self  -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-protector-all
CFLAGS_TEST = -c

# Platform symbol
CDEFINES += -DDarwin

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I ${HDF5_PATH}/include -DH5_USE_18_API

CFLAGS_NCMPI =

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   =  -lgomp -o
LFLAGS_DEBUG =  -lgomp -g -o
LFLAGS_TEST  =  -o


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

LIB_OPT   = #-L${MESA_DIR}/lib $(LOAD_MESA_MICRO)
LIB_DEBUG = #-L${MESA_DIR}/lib $(LOAD_MESA_MICRO)
LIB_TEST  =

#LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 /usr/lib64/libz.a
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5  -lz
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE

LIB_PAPI  =
#LIB_MATH  = -ldfftw -ldrfftw

LIB_MPI   = -L/usr/local/opt/mpich/lib
#LIB_NCMPI = -L $(NCMPI_PATH)/lib -lpnetcdf
LIB_MPE   =

#LIB_MESA = ${LOAD_OTHER}

LIB_STDCXX   = #Standard C++ library
FFLAGS_AMREX = -I${AMREX_PATH}/include -DSKIP_PROBLEMATIC_CODE_FOR_NOW
LIB_AMREX    = ${AMREX_PATH}/lib/libamrex.a

LIB_LAPACK = -I/usr/local/opt/lapack/lib

CONFIG_LIB = ##-L /System/Library/Frameworks/Accelerate.framework -Wl,-framework -Wl,Accelerate

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

ifeq ($(FLASHBINARY),true)

#-openmp option: this is very strange because this file contains no openmp.
FFLAGS_WO_OPENMP = $(patsubst $(OPENMP),-frecursive,$(FFLAGS))
#hy_uhd_unsplit.o : %.o : %.F90
#			     $(FCOMP) $(FFLAGS) $(F90FLAGS) $(FDEFINES)	$<
mpi_amr_1blk_guardcell.o : %.o : %.F90
			     $(FCOMP) $(FFLAGS) -O0 $(F90FLAGS) $(FDEFINES)  $<
endif
