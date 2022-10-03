# FLASH makefile definitions for Summit

ifdef OLCF_GCC_ROOT
    PE_ENV=gnu
else ifdef OLCF_PGI_ROOT
    PE_ENV=pgi
else ifdef OLCF_NVHPC_ROOT
    PE_ENV=nvhpc
else ifdef OLCF_XL_ROOT
    PE_ENV=ibm
else
    $(error Cannot determine compiler module. \
	    Load the gcc, pgi, nvhpc, or xl module---e.g. "module load gcc" \
	    ---or specify a specific makefile with the setup argument---e.g. "-makefile=gcc")
endif

#----------------------------------------------------------------------------
# Set the HDF5/MPI library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH   = ${OLCF_OMPI_ROOT}
AMREX_PATH = ${OLCF_AMREX${NDIM}D_ROOT}
BITTREE_PATH = ${OLCF_BITTREE${NDIM}D_ROOT}
HDF5_PATH  = ${OLCF_HDF5_ROOT}
HYPRE_PATH = ${OLCF_HYPRE_ROOT}
CUDA_PATH  = ${OLCF_CUDA_ROOT}
MAGMA_PATH = ${OLCF_MAGMA_ROOT}

ZLIB_PATH  =
PAPI_PATH  = ${OLCF_PAPI_ROOT}
PAPI_FLAGS =

NCMPI_PATH = ${OLCF_PARALLEL_NETCDF_ROOT}
MPE_PATH   =

ifdef OLCF_ESSL_ROOT
    HYPRE_VERSION = essl

    ESSL_PATH     = ${OLCF_ESSL_ROOT}
    CFLAGS_ESSL   = -I${ESSL_PATH}/include
    FFLAGS_ESSL   = -I${ESSL_PATH}/include
    LIB_ESSL      = -L${ESSL_PATH}/lib64 -lessl

    LAPACK_PATH   = ${ESSL_PATH}
    CFLAGS_LAPACK = ${CFLAGS_ESSL}
    FFLAGS_LAPACK = ${FFLAGS_ESSL}
    LIB_LAPACK    = ${LIB_ESSL}
else
    HYPRE_VERSION = default

    LAPACK_PATH   =
    CFLAGS_LAPACK =
    FFLAGS_LAPACK =
    LIB_LAPACK    =
endif

ifdef OLCF_NETLIB_LAPACK_ROOT
    LAPACK_PATH    = ${OLCF_NETLIB_LAPACK_ROOT}
    CFLAGS_LAPACK += -I${LAPACK_PATH}/include
    FFLAGS_LAPACK += -I${LAPACK_PATH}/include
    LIB_LAPACK    += -L${LAPACK_PATH}/lib64 -llapack -lblas
endif

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

## Compiler-specific flags
ifdef OLCF_GCC_ROOT

    # pre-processor flag
    MDEFS        =
    PP           = -D

    # generic flags
    OPENMP       = -fopenmp

    OPT_FLAGS    = -g -Ofast -funroll-loops -fprefetch-loop-arrays
    TEST_FLAGS   = -g -O3
    DEBUG_FLAGS  = -g -Og

    # Fortran-specific flags
    OPT_FFLAGS   =
    TEST_FFLAGS  =
    DEBUG_FFLAGS = -fcheck=bounds,do,mem,pointer -ffpe-trap=invalid,zero,overflow -fbacktrace

    F90FLAGS     = -fdefault-real-8 -fdefault-double-8 -fimplicit-none -ffree-line-length-none -fallow-argument-mismatch -cpp
    f90FLAGS     = ${F90FLAGS}
    F77FLAGS     = -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fallow-argument-mismatch -cpp
    f77FLAGS     = ${F77FLAGS}

    FFLAGS_OACC  = -fopenacc

    # C-specific flags
    OPT_CFLAGS   =
    TEST_CFLAGS  =
    DEBUG_CFLAGS =

    CFLAGS_OACC  = -fopenacc

    # Linker flags
    LIB_OPT      =
    LIB_TEST     =
    LIB_DEBUG    =

    LIB_OACC     =
    LIB_MASS     =

else ifdef OLCF_PGI_ROOT

    # pre-processor flag
    MDEFS        =
    PP           = -D

    # generic flags
    OPENMP       = -mp

    OPT_FLAGS    = -g -O2 -Mpreprocess
    TEST_FLAGS   = -g -O1 -Mpreprocess
    DEBUG_FLAGS  = -g -O0 -Mpreprocess -Mbounds -Mnoopenmp

    # Fortran-specific flags
    OPT_FFLAGS   =
    TEST_FFLAGS  =
    DEBUG_FFLAGS =

    F90FLAGS     = -r8 -i4
    f90FLAGS     = ${F90FLAGS}
    F77FLAGS     = -r8 -i4 -Mfixed
    f77FLAGS     = ${F77FLAGS}

    FFLAGS_OACC  = -acc -ta=tesla:cc70,ptxinfo -Minfo=accel

    # C-specific flags
    OPT_CFLAGS   =
    TEST_CFLAGS  =
    DEBUG_CFLAGS =

    CFLAGS_OACC  = -acc -ta=tesla:cc70,ptxinfo -Minfo=accel

    # Linker flags
    LIB_OPT      = -pgc++libs
    LIB_TEST     = -pgc++libs
    LIB_DEBUG    = -pgc++libs

    LIB_OACC     = -acc -ta=tesla:cc70,ptxinfo -acclibs
    LIB_MASS     =


else ifdef OLCF_NVHPC_ROOT

    # pre-processor flag
    MDEFS        =
    PP           = -D

    # generic flags
    OPENMP       = -mp=multicore

    OPT_FLAGS    = -g -O2 -Mpreprocess
    TEST_FLAGS   = -g -O1 -Mpreprocess
    DEBUG_FLAGS  = -g -O0 -Mpreprocess -Mbounds -Mnoopenmp

    # Fortran-specific flags
    OPT_FFLAGS   =
    TEST_FFLAGS  =
    DEBUG_FFLAGS =

    F90FLAGS     = -r8 -i4
    f90FLAGS     = ${F90FLAGS}
    F77FLAGS     = -r8 -i4 -Mfixed
    f77FLAGS     = ${F77FLAGS}

    FFLAGS_OACC  = -acc -gpu=cc70,ptxinfo -Minfo=accel
    FFLAGS_OMP_OL= ${OPENMP} -mp=gpu -gpu=cc70,ptxinfo -Minfo=accel

    # C-specific flags
    OPT_CFLAGS   =
    TEST_CFLAGS  =
    DEBUG_CFLAGS =

    CFLAGS_OACC  = -acc -gpu=cc70,ptxinfo -Minfo=accel
    CFLAGS_OMP_OL= ${OPENMP} -mp=gpu -gpu=cc70,ptxinfo -Minfo=accel

    # Linker flags
    LIB_OPT      = -pgc++libs
    LIB_TEST     = -pgc++libs
    LIB_DEBUG    = -pgc++libs

    LIB_OACC     = -acc -gpu=cc70,ptxinfo -acclibs
    LIB_OMP_OL   = ${OPENMP} -mp=gpu -gpu=cc70,ptxinfo
    LIB_MASS     =

else ifdef OLCF_XL_ROOT

    # pre-processor flag
    MDEFS        = -WF,
    PP           = -D

    # generic flags
    OPENMP       = -qsmp=omp:noauto

    OPT_FLAGS    = -g -O2 -qarch=pwr9 -qtune=pwr9 -w
    TEST_FLAGS   = -g -O2 -qarch=pwr9 -qtune=pwr9 -w -qstrict=all
    DEBUG_FLAGS  = -g -O0 -qnosmp -qstrict=all \
	           -qfloat=rngchk -qcheck=all:nounset \
                   -qflttrap=enable:invalid:nanq:overflow:zerodivide -qsigtrap=xl__trcedump

    # Fortran-specific flags
    OPENMP_FORTRAN = ${OPENMP} -qnosave -qthreaded

    OPT_FFLAGS   =
    TEST_FFLAGS  =
    DEBUG_FFLAGS = -qflag=i:w

    F90FLAGS     = -qintsize=4 -qrealsize=8 -qzerosize -qport=c_loc -qundef -qsuppress=cmpmsg
    f90FLAGS     = ${F90FLAGS}
    F77FLAGS     = ${F90FLAGS} -qfixed
    f77FLAGS     = ${F77FLAGS}

    FFLAGS_OACC  =
    FFLAGS_OMP_OL= ${OPENMP} -qoffload

    # C-specific flags
    OPENMP_C     = ${OPENMP}

    OPT_CFLAGS   = ${PP}IBM
    TEST_CFLAGS  = ${PP}IBM
    DEBUG_CFLAGS = ${PP}IBM

    CFLAGS_OACC  =
    CFLAGS_OMP_OL= ${OPENMP} -qoffload

    # Linker flags
    OPENMP_LINK  = ${OPENMP} -qnosave -qthreaded

    LIB_OPT      =
    LIB_TEST     = -pg
    LIB_DEBUG    =

    LIB_OACC     =
    LIB_OMP_OL   = ${OPENMP_LINK} -qoffload
    LIB_MASS     = -lmass

endif

## Compiler-independent flags
FFLAGS_OPT   = -c ${OPT_FLAGS} ${OPT_FFLAGS}
FFLAGS_TEST  = -c ${TEST_FLAGS} ${TEST_FFLAGS}
FFLAGS_DEBUG = -c ${DEBUG_FLAGS} ${DEBUG_FFLAGS}

FFLAGS_AMREX = -I${AMREX_PATH}/include
FFLAGS_BITTREE = -I${BITTREE_PATH}/include
FFLAGS_HDF5  = -I${HDF5_PATH}/include 
FFLAGS_NCMPI = -I${NCMPI_PATH}/include
FFLAGS_HYPRE = -I${HYPRE_PATH}/include ${FFLAGS_LAPACK}
FFLAGS_CUDA  = -I${CUDA_PATH}/include
FFLAGS_MAGMA = -I${MAGMA_PATH}/include
FFLAGS_UNIFYFS = -I${UNIFYFS_ROOT}

CFLAGS_OPT   = -c ${OPT_FLAGS} ${OPT_CFLAGS}
CFLAGS_TEST  = -c ${TEST_FLAGS} ${TEST_CFLAGS}
CFLAGS_DEBUG = -c ${DEBUG_FLAGS} ${DEBUG_CFLAGS}

CFLAGS_AMREX = -I${AMREX_PATH}/include
CFLAGS_BITTREE = -I${BITTREE_PATH}/include
CFLAGS_HDF5  = -I$(HDF5_PATH)/include 
CFLAGS_NCMPI = -I$(NCMPI_PATH)/include
CFLAGS_HYPRE = -I${HYPRE_PATH}/include ${CFLAGS_LAPACK}
CFLAGS_CUDA  = -I${CUDA_PATH}/include
CFLAGS_MAGMA = -I${MAGMA_PATH}/include
CFLAGS_UNIFYFS = -I${UNIFYFS_ROOT}

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

LIB_MPI   =

LIB_STDCXX = -lstdc++
LIB_AMREX = -L${AMREX_PATH}/lib -lamrex -lstdc++
LIB_BITTREE = -L${BITTREE_PATH}/lib -lbittree
LIB_HDF5  = -L${HDF5_PATH}/lib -lhdf5_fortran -lhdf5
LIB_NCMPI = -L${NCMPI_PATH}/lib -lpnetcdf
LIB_MATH  = ${LIB_ESSL}
LIB_HYPRE = -L${HYPRE_PATH}/lib -lHYPRE ${LIB_LAPACK}
LIB_CUDA  = -L${CUDA_PATH}/lib64 -lcusparse -lcusolver -lcublas -lcudart -lcuda
LIB_MAGMA = -L$(MAGMA_PATH)/lib -lmagma -Wl,-rpath,$(MAGMA_PATH)/lib
LIB_UNIFYFS = -L${UNIFYFS_ROOT} -lunifyfs_mpi_gotcha

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

ifdef OLCF_PGI_ROOT

## PGI compiler does not like -O2 or above with this file for some reason
mpi_amr_prolong.o : mpi_amr_prolong.F90
	$(ECHO-COMPILING)
	$(FCOMP) $(patsubst -fast,-O1,$(FFLAGS)) $(F90FLAGS) $(FDEFINES) $< -o $(addsuffix .o,$(basename $@))

endif
