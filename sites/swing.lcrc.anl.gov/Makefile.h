# Flash-X Makefile header for swing.lcrc.gov (NVHPC)
# tested with nvhpc/22.3
#   $ module load nvhpc/22.3-zvp5ypz
#
#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH   =
AMREX_PATH =
HDF5_PATH  = /lcrc/project/Flash-X/soft/nvhpc_22.3/hdf5/1.12.2
HYPRE_PATH =
CUDA_PATH  = ${CUDA_HOME}


#----------------------------------------------------------------------------
# Compiler and linker commands
#
#   Use the MPICH wrappers around the compilers -- these will automatically
#   load the proper libraries and include files.  Version of MPICH prior
#   to 1.2.2 (?) do not recognize .F90 as a valid Fortran file extension.
#   You need to edit mpif90 and add .F90 to the test of filename extensions,
#   or upgrade your MPICH.
#----------------------------------------------------------------------------

FCOMP   = mpifort
CCOMP   = mpicc
CPPCOMP = mpicxx -std=c++11
CUCOMP  = nvcc
LINK    = mpifort

# pre-processor flag
MDEFS =
PP    = -D

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

OPENMP = -mp
OMP_CUDA = -Minfo=accel,mp -gpu=cc80

FFLAGS_AMREX = -I${AMREX_PATH}/include
FFLAGS_HDF5  = -I${HDF5_PATH}/include ${MDEFS}${PP}H5_USE_18_API
FFLAGS_HYPRE = -I${HYPRE_PATH}/include ${FFLAGS_LAPACK}
FFLAGS_CUDA  = -I${CUDA_PATH}/include
FFLAGS_OMP_OL= -mp=gpu ${OMP_CUDA}

OPT_FLAGS    = -g -O2 -Mpreprocess
TEST_FLAGS   = -g -O1 -Mpreprocess
DEBUG_FLAGS  = -g -O0 -Mpreprocess -Mbounds -Mnoopenmp

FFLAGS_OPT   = -c ${OPT_FLAGS}
FFLAGS_TEST  = -c ${TEST_FLAGS}
FFLAGS_DEBUG = -c ${DEBUG_FLAGS}

F90FLAGS     = -r8 -i4
f90FLAGS     = ${F90FLAGS}
F77FLAGS     = -r8 -i4 -Mfixed
f77FLAGS     = ${F77FLAGS}

CFLAGS_AMREX = -I${AMREX_PATH}/include
CFLAGS_HDF5  = -I$(HDF5_PATH)/include ${PP}H5_USE_18_API
CFLAGS_HYPRE = -I${HYPRE_PATH}/include ${CFLAGS_LAPACK}
CFLAGS_CUDA  = -I${CUDA_PATH}/include
CFLAGS_OMP_OL= -mp=gpu ${OMP_CUDA}

CFLAGS_OPT   = -c ${OPT_FLAGS}
CFLAGS_TEST  = -c ${TEST_FLAGS}
CFLAGS_DEBUG = -c ${DEBUG_FLAGS}

CU_FLAGS     = -c -g -O2 -m64 -gencode arch=compute_70,code=sm_70

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90 .cu

#----------------------------------------------------------------------------
# Linker flags
#----------------------------------------------------------------------------
LFLAGS_OPT   = ${OPT_FLAGS} -o
LFLAGS_TEST  = ${TEST_FLAGS} -o
LFLAGS_DEBUG = ${DEBUG_FLAGS} -o

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

LIB_MPI =

LIB_STDCXX = -lstdc++
LIB_AMREX = -L${AMREX_PATH}/lib -lamrex -lstdc++
LIB_HDF5  = -L${HDF5_PATH}/lib -rpath=${HDF5_PATH}/lib -lhdf5 ${PP}H5_USE_18_API
LIB_MATH  = ${LIB_ESSL}
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE ${LIB_LAPACK}
LIB_CUDA  = -cudalib=cusparse,cusolver,cublas
LIB_OMP_OL= -mp=gpu ${OMP_CUDA}

LIB_OPT   = -pgc++libs
LIB_DEBUG = -pgc++libs
LIB_TEST  = -pgc++libs

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
AWK = awk
CAT = cat
