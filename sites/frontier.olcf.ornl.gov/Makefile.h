# Flash-X makefile definitions for Frontier

ifeq ($(findstring $(PE_ENV),GNU CRAY),)
$(error Your environment "$(PE_ENV)" is invalid.  It must be "CRAY" or "GNU")
else
$(warning You are using the "$(PE_ENV)" environment)
endif

LC_PE_ENV = $(shell echo ${PE_ENV} | tr A-Z a-z)

#----------------------------------------------------------------------------
# Set the library paths -- these need to be updated for your system
#----------------------------------------------------------------------------

MPI_PATH     =
AMREX_PATH   = ${OLCF_AMREX${NDIM}D_ROOT}
HDF5_PATH    =
HYPRE_PATH   =
HIPFORT_PATH = ${OLCF_HIPFORT_ROOT}

ifdef LAPACK_PATH
  INC_LAPACK = -I${LAPACK_PATH}/include
  LIB_LAPACK = -L${LAPACK_PATH}/lib64 -llapack -lblas
else
  LAPACK_PATH = ${CRAY_LIBSCI_PREFIX_DIR}
  INC_LAPACK =# -I${LAPACK_PATH}/include
  LIB_LAPACK =# -L${LAPACK_PATH}/lib -lsci_cray
endif

ROCM_PATH ?= ${OLCF_ROCM_ROOT}
ifneq (${ROCM_PATH},)
  INC_ROCM = -I$(ROCM_PATH)/include -I$(HIPFORT_PATH)/include/hipfort/amdgcn
  LIB_ROCM = -L$(HIPFORT_PATH)/lib -lhipfort-amdgcn -L$(ROCM_PATH)/lib -lrocsparse -lrocsolver -lrocblas -lhipblas -lhipsparse -lamdhip64
else
  ifeq (${USE_ROCM},TRUE)
    $(error Cannot resolve ROCM path. \
            Load the ROCM module---e.g. "module load rocm"--- \
            or define ROCM_PATH variable for a valid ROCM build.)
  endif
endif

MAGMA_PATH ?= ${OLCF_MAGMA_ROOT}
ifneq (${MAGMA_PATH},)
  INC_MAGMA = -I${MAGMA_PATH}/include
  LIB_MAGMA = -L$(MAGMA_PATH)/lib -lmagma
else
  ifeq (${USE_MAGMA},TRUE)
    $(error Cannot resolve MAGMA path. \
            Load the MAGMA module---e.g. "module load magma"--- \
            or define MAGMA_PATH variable for a valid MAGMA build.)
  endif
endif

#----------------------------------------------------------------------------
# Compiler and linker commands
#----------------------------------------------------------------------------

FCOMP   = ftn
CCOMP   = cc
CPPCOMP = CC -std=c++11
HIPCOMP = hipcc
CUCOMP  = ${HIPCOMP}
LINK    = ftn

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
ifeq ($(PE_ENV),GNU)

    # pre-processor flag
    MDEFS         =
    PP            = -D

    # generic flags
    OPENMP        = -fopenmp

    OPT_FLAGS     = -g -O3
    TEST_FLAGS    = -g -O2
    DEBUG_FLAGS   = -g -Og

    # Fortran-specific flags
    OPT_FFLAGS    =
    TEST_FFLAGS   =
    DEBUG_FFLAGS  = -fcheck=bounds,do,mem,pointer -ffpe-trap=invalid,zero,overflow -fbacktrace

    F90FLAGS      = -fdefault-real-8 -fdefault-double-8 -fimplicit-none -ffree-line-length-none -fallow-argument-mismatch -cpp
    f90FLAGS      = ${F90FLAGS}
    F77FLAGS      = -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fallow-argument-mismatch -cpp
    f77FLAGS      = ${F77FLAGS}

    FFLAGS_OACC   = -fopenacc
    FFLAGS_OMP_OL = -fopenmp

    # C-specific flags
    OPT_CFLAGS    =
    TEST_CFLAGS   =
    DEBUG_CFLAGS  =

    CFLAGS_OACC   = -fopenacc
    CFLAGS_OMP_OL = -fopenmp

    # Linker flags
    LIB_OPT       =
    LIB_TEST      =
    LIB_DEBUG     =

    LIB_OACC      =
    LIB_OMP_OL    =

else ifeq ($(PE_ENV),CRAY)

    # pre-processor flag
    MDEFS         =
    PP            = -D

    # generic flags
    OPENMP        = -fopenmp

    OPT_FLAGS     = -O2
    TEST_FLAGS    = -O1 -em 
    DEBUG_FLAGS   = -g -O0

    # Fortran-specific flags
    OPT_FFLAGS    =
    TEST_FFLAGS   =
    DEBUG_FFLAGS  = -Ktrap=fp

    F90FLAGS      = -s real64 -s integer32 -eI -eZ -ef -f free
    f90FLAGS      = ${F90FLAGS}
    F77FLAGS      = -s real64 -s integer32 -eI -eZ -ef -f fixed
    f77FLAGS      = ${F77FLAGS}

    FFLAGS_OACC   = -hacc
    FFLAGS_OMP_OL = -fopenmp

    # C-specific flags
    OPT_CFLAGS    = -Wno-error=int-conversion
    TEST_CFLAGS   = -Wno-error=int-conversion
    DEBUG_CFLAGS  = -ffpe-trap=fp -Wno-error=int-conversion

    CFLAGS_OACC   = -hacc
    CFLAGS_OMP_OL = -fopenmp

    # Linker flags
    LIB_OPT       =
    LIB_TEST      =
    LIB_DEBUG     =

    LIB_OACC      = -hacc
    LIB_OMP_OL    = -fopenmp

endif

## Compiler-independent flags
FFLAGS_OPT     = -c ${OPT_FLAGS} ${OPT_FFLAGS}
FFLAGS_TEST    = -c ${TEST_FLAGS} ${TEST_FFLAGS}
FFLAGS_DEBUG   = -c ${DEBUG_FLAGS} ${DEBUG_FFLAGS}

FFLAGS_AMREX   = -I${AMREX_PATH}/include
FFLAGS_HDF5    = ${MDEFS}${PP}H5_USE_18_API
FFLAGS_UNIFYFS =# -I${UNIFYFS_ROOT}
FFLAGS_NCMPI   =
FFLAGS_LAPACK  = ${INC_LAPACK}
FFLAGS_HYPRE   =
FFLAGS_ROCM    = ${INC_ROCM}
FFLAGS_CUDA    = ${FFLAGS_ROCM}
FFLAGS_MAGMA   = ${INC_MAGMA}

CFLAGS_OPT     = -c ${OPT_FLAGS} ${OPT_CFLAGS}
CFLAGS_TEST    = -c ${TEST_FLAGS} ${TEST_CFLAGS}
CFLAGS_DEBUG   = -c ${DEBUG_FLAGS} ${DEBUG_CFLAGS}

CFLAGS_AMREX   = -I${AMREX_PATH}/include
CFLAGS_HDF5    = ${PP}H5_USE_18_API
CFLAGS_UNIFYFS =# -I${UNIFYFS_ROOT}
CFLAGS_NCMPI   =
CFLAGS_LAPACK  = ${INC_LAPACK}
CFLAGS_HYPRE   =
CFLAGS_ROCM    = ${INC_ROCM}
CFLAGS_CUDA    = ${CFLAGS_ROCM}
CFLAGS_MAGMA   = ${INC_MAGMA}

HIP_FLAGS    = -c -g -O3 -std=c++11 --amdgpu-target=gfx90a \
               -I${MPICH_DIR}/include -L${MPICH_DIR}/lib -lmpi \
               -L${CRAY_MPICH_ROOTDIR}/gtl/lib -lmpi_gtl_hsa
CU_FLAGS     = ${HIP_FLAGS}

.SUFFIXES: .o .c .f .F .h .fh .F90 .f90 .cu .hip

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

LIB_STDCXX  = -lstdc++
LIB_AMREX   = -L${AMREX_PATH}/lib -lamrex ${LIB_STDCXX}
LIB_HDF5    = -lhdf5_fortran -lhdf5
LIB_UNIFYFS =# -L${UNIFYFS_ROOT} -lunifyfs_mpi_gotcha
LIB_NCMPI   =
LIB_HYPRE   =
LIB_CUDA    = ${LIB_ROCM}

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

ifeq ($(PE_ENV),CRAY)
local_tree_build.o : local_tree_build.F90
	$(ECHO-COMPILING)
	$(FCOMP) $(FFLAGS) $(F90FLAGS) -h ipa1 $(FDEFINES) $< -o $(addsuffix .o,$(basename $@))
endif
