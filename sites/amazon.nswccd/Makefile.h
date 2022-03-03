# FLASH makefile definitions for x86-64 Linux (GNU compilers)
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

#MPI_PATH = /usr
#MPI_PATH = /opt/scyld/mpich2/1.3.2/intel
MPI_PATH = /opt/scyld/mpich2/1.3.2/gnu

HDF4_PATH  =
HDF5_PATH  = /san/home/delaneyk/flash_CO/hdf5/hdf5-1.8.10-linux-x86_64-static
ZLIB_PATH  = 

PAPI_PATH  =
PAPI_FLAGS =

NCMPI_PATH = 
MPE_PATH   =

#HYPRE_PATH = /san/home/delaneyk/flash_CO/HYPRE/hypre-2.9.0b/Hypre2.9
HYPRE_PATH = /san/home/delaneyk/flash_CO/HYPRE_GNU/hypre-2.9.0b/src/hypre

#IFORT_PATH = /opt/intel/Compiler/11.0/083
IFORT_PATH = 

export cur-dir := $(shell pwd)

# Set the location of top directory
export setup_dir = $(cur-dir)


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------

#FCOMP   = ${MPI_PATH}/bin/mpif90
#CCOMP   = ${MPI_PATH}/bin/mpicc
#CPPCOMP = ${MPI_PATH}/bin/mpicxx
#LINK    = ${MPI_PATH}/bin/mpif90

#AMAZON ... f90 is ifort ans rest are gcc-based!!!
#FCOMP   = mpif90 
#CCOMP   = mpicc
#CPPCOMP = mpiCC
#LINK    = mpif90

FCOMP   = ~/mpif90 
CCOMP   = ~/mpicc
CPPCOMP = ~/mpicxx
LINK    = ~/mpif90

#CCOMP   = ${MPI_PATH}/../gnu/bin/mpicc
#CPPCOMP = ${MPI_PATH}/../gnu/bin/mpicxx

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


## IFORT flags:
#FFLAGS_OPT   = -c  -r8 -i4 -O3 -unroll -align -pad -ip \
#-real_size 64 -diag-disable 10120 -vec_report0
#
#FFLAGS_DEBUG = -c -g -r8 -i4 -O0 -check bounds -check format \
#-check output_conversion -real_size 64 -check uninit \
#-traceback -fp-stack-check -diag-disable 10120 -fpe0 -check pointers -extend-source

# GFORT flags
FFLAGS_OPT =  -c -O2 -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -Wuninitialized

FFLAGS_DEBUG = -ggdb -c -fdefault-real-8 -fdefault-double-8 \
-ffree-line-length-none -pedantic -Wall -Wextra -Waliasing \
-Wsurprising -Wconversion -Wunderflow \
-ffpe-trap=invalid,zero,overflow -fbounds-check \
-fbacktrace -fdump-core -finit-real=nan \
-finit-integer=-999999 -fimplicit-none

 

# GCC/G++ FLAGS
CFLAGS_OPT =  -c -O2 -Wuninitialized 


#CFLAGS_DEBUG = -c -O2 -Wuninitialized 

CFLAGS_DEBUG =  -ggdb -c -Wno-div-by-zero -Wundef  \
-Wconversion -Wstrict-prototypes -Wunreachable-code \
-pedantic -Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
-Wunsafe-loop-optimizations -Wpadded -fstack-check -fno-stack-protector
 
CPPLAGS_DEBUG =  -ggdb -c  -pedantic -fstack-check -fstack-protector-all

#CPPLAGS_DEBUG =  -ggdb -c -Wno-div-by-zero -Wundef  \
#-Wconversion -Wstrict-prototypes -Wunreachable-code \
#-Wall -Wextra -Winit-self -ftree-vrp -Wfloat-equal \
#-Wunsafe-loop-optimizations -Wpadded -fstack-check -fstack-protector-all

CFLAGS_TEST = -c

# if we are using HDF5, we need to specify the path to the include files
CFLAGS_HDF5 = -I${HDF5_PATH}/include -DH5_USE_18_API

CFLAGS_HYPRE = -I${HYPRE_PATH}/include 
FFLAGS_HYPRE = -I${HYPRE_PATH}/include 

#----------------------------------------------------------------------------
# Linker flags
#
#  There is a seperate version of the linker flags for each of the _OPT,
#  _DEBUG, and _TEST cases.
#----------------------------------------------------------------------------

LFLAGS_OPT   = -o
LFLAGS_DEBUG = -o
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

LIB_HDF4  = 
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5 -lz -L/san/home/delaneyk/flash_CO/szip/szip-2.1/szip/lib -lsz 
LIB_PAPI  =
LIB_MATH  = 

LIB_MPI   = 
LIB_MPE   =

#LIB_IFORT = -L${IFORT_PATH}/lib/intel64 -limf -lm -lifcoremt
#LIB_IFORT = -L/usr/local/lib -lmpfr 
LIB_IFORT = 

LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE  

LIB_STDCXX = -lstdc++

#Specify TEC_PLOT=YES in order to link the tec plot library.
TEC_PLOT=YES
ifeq ($(TEC_PLOT), YES)
CONFIG_LIB = -I${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -L${setup_dir}/../source/Simulation/SimulationMain/INavierStokes -ltecio -lstdc++
endif

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

CFLAGS_WO_WARNALL = $(CFLAGS_OPT)

bittree_core.o : %.o : %.cxx
	$(CCOMP) $(CFLAGS_WO_WARNALL) $(CDEFINES)          $<	


endif
