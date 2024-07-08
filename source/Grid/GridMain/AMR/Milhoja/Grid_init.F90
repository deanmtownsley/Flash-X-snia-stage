!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

#include "Milhoja.h"

#include "constants.h"
#include "Simulation.h"

#ifdef GRID_USE_AMREX_BACKEND
#ifndef MILHOJA_AMREX_GRID_BACKEND
#error "Milhoja library is not using the AMReX grid backend"
#endif
#endif

!> @ingroup GridMilhoja
!! @stubref{Grid_init}
!!
!! @brief Concrete implementation of Grid_init
!!
!! @attention
!! While its interface is fixed, this routine is a work in progress.
!! Specifically, it will grow as new initialization work is required to add in
!! new Milhoja Grid unit functionality.  Ideally, this will only be altered when
!! a new test or updated test can confirm that alterations are correct.
!!
!! @todo Confirm that the Milhoja types match the Flash-X types.  Look at the
!!       testing code that used to be in Milhoja's types module.  How does
!!       Milhoja confirm that its types match those of the Grid backend?
!! @todo Figure out the whole max refine deal & implement correctly
!! @todo Pass nrefs to Milhoja
!! @todo What BC values should be passed to Milhoja above NDIM?
!! @todo Are the domain values above NDIM set in gr_globalDomain correct?
!! @todo The setup tool needs to tell us which grid backend has been selected
!!       so that we know which interpolator variable to config Milhoja with.
!! @todo Grid_init needs to be redesigned with constructor inheritance.  For
!! instance, it appears that all Grid implementations should call
!! gr_initGeometry.  However, gr_initGeometry requires that some grid local
!! variables be loaded prior to being called.  This implies a base class with
!! common data members that need to be set by a constructor for the base class.
!! Rather than having some base constructor carry out all this work, we presently
!! expect each grid implementation to do this.  If we were to implement a data
!! member/constructor inheritance scheme, then writing and maintaining Grid_inits at this
!! level would be easier and less error prone.
!! @todo The runtime parameter geometryOverride should probably be made to have
!!       the same effects as in the Amrex Grid implementation (see gr_amrexInit).
!!       For now, code below aborts instead if gr_geometryOverride is TRUE.
subroutine Grid_init()
    use milhoja_types_mod,           ONLY : MILHOJA_INT, &
                                            MILHOJA_REAL
    use milhoja_grid_mod,            ONLY : milhoja_grid_init, &
                                            milhoja_grid_getCoordinateSystem, &
                                            milhoja_grid_getDomainBoundBox, &
                                            milhoja_grid_getMaxFinestLevel, &
                                            milhoja_grid_getBlockSize, &
                                            milhoja_grid_getDomainDecomposition, &
                                            milhoja_grid_getNGuardcells, &
                                            milhoja_grid_getNCcVariables, &
                                            milhoja_grid_getNFluxVariables

    use Grid_data
    use Grid_interface,              ONLY : Grid_getDeltas
    use gr_milhojaInterface,         ONLY : gr_checkMilhojaError, &
                                            gr_fillPhysicalBcCallback, &
                                            gr_markRefineDerefineCallback
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get, &
                                            RuntimeParameters_mapStrToInt
    use Driver_interface,            ONLY : Driver_getMype, &
                                            Driver_getNumProcs, &
                                            Driver_getComm, &
                                            Driver_envGetScalar, &
                                            Driver_abort
    use Logfile_interface,           ONLY : Logfile_stampMessage
    use Simulation_interface,        ONLY : Simulation_getVarnameType

    implicit none

    character(len=2), parameter :: COL_HDR(1:MDIM) = ['dx', 'dy', 'dz']

    character(len=*), parameter :: REF_VAR_NAME      = 'refine_var_'
    character(len=*), parameter :: REF_CUTOFF_NAME   = 'refine_cutoff_'
    character(len=*), parameter :: DEREF_CUTOFF_NAME = 'derefine_cutoff_'
    character(len=*), parameter :: REF_FILTER_NAME   = 'refine_filter_'

    character(len=MAX_STRING_LENGTH) :: paramString
    character(len=MAX_STRING_LENGTH) :: refVarString
    character(len=MAX_STRING_LENGTH) :: refCutoffString
    character(len=MAX_STRING_LENGTH) :: derefCutOffString
    character(len=MAX_STRING_LENGTH) :: refFilterString
    character(len=MAX_STRING_LENGTH) :: xl_bcString, xr_bcString
    character(len=MAX_STRING_LENGTH) :: yl_bcString, yr_bcString
    character(len=MAX_STRING_LENGTH) :: zl_bcString, zr_bcString
    character(len=MAX_STRING_LENGTH) :: eosModeString
    character(len=MAX_STRING_LENGTH) :: ccInterpolatorString

    integer :: nBlocksX, nBlocksY, nBlocksZ
    integer :: lRefineMax
    integer :: nrefs
    real    :: xMin, xMax
    real    :: yMin, yMax
    real    :: zMin, zMax
    integer :: refVar
    integer :: nonrep
    integer :: ccInterpolator

    character(len=20) :: fmtStr
    integer           :: level
    real              :: deltas(1:MDIM)

    real(MILHOJA_REAL)   :: MH_xMin, MH_xMax
    real(MILHOJA_REAL)   :: MH_yMin, MH_yMax
    real(MILHOJA_REAL)   :: MH_zMin, MH_zMax
    integer(MILHOJA_INT) :: MH_logRank
    integer(MILHOJA_INT) :: MH_nxb, MH_nyb, MH_nzb
    integer(MILHOJA_INT) :: MH_nBlocksX, MH_nBlocksY, MH_nBlocksZ
    integer(MILHOJA_INT) :: MH_maxRefinementLevel
    integer(MILHOJA_INT) :: MH_nGuard
    integer(MILHOJA_INT) :: MH_nCcVars
    integer(MILHOJA_INT) :: MH_nFluxVars

    integer(MILHOJA_INT) :: MH_nBlocksX_SC
    integer(MILHOJA_INT) :: MH_nBlocksY_SC
    integer(MILHOJA_INT) :: MH_nBlocksZ_SC

    real(MILHOJA_REAL)   :: MH_domainLo(1:MDIM)
    real(MILHOJA_REAL)   :: MH_domainHi(1:MDIM)

    integer(MILHOJA_INT) :: MH_loBCs(1:MDIM)
    integer(MILHOJA_INT) :: MH_hiBCs(1:MDIM)

    integer(MILHOJA_INT) :: MH_coordSys
    integer(MILHOJA_INT) :: MH_ccInterpolator
    integer(MILHOJA_INT) :: MH_level

    integer(MILHOJA_INT) :: MH_ierr

    integer :: i

    !!!!!----- CONFIRM CORRECT MILHOJA LIBRARY LINKAGE
    if (     (MILHOJA_MDIM /= MDIM) &
        .OR. (MILHOJA_NDIM /= NDIM) &
        .OR. (MILHOJA_K1D  /= K1D)  &
        .OR. (MILHOJA_K2D  /= K2D)  &
        .OR. (MILHOJA_K3D  /= K3D)) then
        CALL Driver_abort("[Grid_init] Simulation linked against incorrect Milhoja library")
    end if

    !!!!!----- MPI & OpenMP Environment
    CALL Driver_getMype(GLOBAL_COMM,     gr_globalMe)
    CALL Driver_getNumProcs(GLOBAL_COMM, gr_globalNumProcs)
    CALL Driver_getComm(GLOBAL_COMM,     gr_globalComm)

    CALL Driver_getMype(MESH_COMM,     gr_meshMe)
    CALL Driver_getNumProcs(MESH_COMM, gr_meshNumProcs)
    CALL Driver_getComm(MESH_COMM,     gr_meshComm)

    CALL Driver_getMype(MESH_ACROSS_COMM,     gr_meshAcrossMe)
    CALL Driver_getNumProcs(MESH_ACROSS_COMM, gr_meshAcrossNumProcs)
    CALL Driver_getComm(MESH_ACROSS_COMM,     gr_meshAcrossComm)

    CALL Driver_envGetScalar("OMP_NUM_THREADS", gr_envOmpNumThreads)

    !!!!!----- Runtime Parameters
    if (gr_globalMe == MASTER_PE) then
        write(*,'(A)') "[Grid] Initializing ..."
    end if

  ! Initialization of gr_geometry etc is done in gr_initGeometry, called below.

  ! DO THIS EARLY - must be before gr_initGeometry is called:
    !get the boundary conditions stored as strings in the par file
    CALL RuntimeParameters_get("xl_boundary_type", xl_bcString)
    CALL RuntimeParameters_get("xr_boundary_type", xr_bcString)
    CALL RuntimeParameters_get("yl_boundary_type", yl_bcString)
    CALL RuntimeParameters_get("yr_boundary_type", yr_bcString)
    CALL RuntimeParameters_get("zl_boundary_type", zl_bcString)
    CALL RuntimeParameters_get("zr_boundary_type", zr_bcString)

    !map the string boundary conditions to integer constants defined in constants.h
    CALL RuntimeParameters_mapStrToInt(xl_bcString, gr_domainBC(LOW, IAXIS))
    CALL RuntimeParameters_mapStrToInt(xr_bcString, gr_domainBC(HIGH,IAXIS))
    CALL RuntimeParameters_mapStrToInt(yl_bcString, gr_domainBC(LOW, JAXIS))
    CALL RuntimeParameters_mapStrToInt(yr_bcString, gr_domainBC(HIGH,JAXIS))
    CALL RuntimeParameters_mapStrToInt(zl_bcString, gr_domainBC(LOW, KAXIS))
    CALL RuntimeParameters_mapStrToInt(zr_bcString, gr_domainBC(HIGH,KAXIS))

    !------------------------------------------------------------------------------
    ! Init geometry first as it can change runtime parameters
    !------------------------------------------------------------------------------
  ! Initialize geometry-related Flash-X runtime parameters,
  ! determine the geometries of  dimensions, and scale
  ! angle value parameters that are expressed in degrees to radians.
    CALL gr_initGeometry()
    if (gr_geometryOverride) then
       CALL Driver_abort("[Grid_init] Runtime parameter geometryOverride not yet supported by Milhoja Grid")
    end if

    !------------------------------------------------------------------------------
    ! Load into local Grid variables all runtime parameters needed by Milhoja
    !------------------------------------------------------------------------------
    ! Don't load any of these into Grid_data variables.  That is done after
    ! initializing Milhoja.
    CALL RuntimeParameters_get('xmin', xMin)
    CALL RuntimeParameters_get('xmax', xMax)
    CALL RuntimeParameters_get('ymin', yMin)
    CALL RuntimeParameters_get('ymax', yMax)
    CALL RuntimeParameters_get('zmin', zMin)
    CALL RuntimeParameters_get('zmax', zMax)

    CALL RuntimeParameters_get("nrefs", nrefs)

    CALL RuntimeParameters_get('lrefine_max', lRefineMax)

    CALL RuntimeParameters_get("nblockx", nBlocksX)
    CALL RuntimeParameters_get("nblocky", nBlocksY)
    CALL RuntimeParameters_get("nblockz", nBlocksZ)

    !----------------------------------------------------------------------------------
    ! Initialize library
    !----------------------------------------------------------------------------------
    CALL RuntimeParameters_get("amrexInterpolator", ccInterpolatorString)
    CALL RuntimeParameters_mapStrToInt(ccInterpolatorString, ccInterpolator)
    if (ccInterpolator == NONEXISTENT) then
      CALL Driver_abort("[Grid_init] Unknown amrexInterpolator runtime parameter value")
    end if

    MH_xMin               = REAL(xMin,           kind=MILHOJA_REAL)
    MH_xMax               = REAL(xMax,           kind=MILHOJA_REAL)
    MH_yMin               = REAL(yMin,           kind=MILHOJA_REAL)
    MH_yMax               = REAL(yMax,           kind=MILHOJA_REAL)
    MH_zMin               = REAL(zMin,           kind=MILHOJA_REAL)
    MH_zMax               = REAL(zMax,           kind=MILHOJA_REAL)

    ! Assume that these have already been confirmed to be
    ! non-negative
    MH_logRank            =  INT(MASTER_PE,      kind=MILHOJA_INT)
    MH_nxb                =  INT(NXB,            kind=MILHOJA_INT)
    MH_nyb                =  INT(NYB,            kind=MILHOJA_INT)
    MH_nzb                =  INT(NZB,            kind=MILHOJA_INT)
    MH_nBlocksX           =  INT(nBlocksX,       kind=MILHOJA_INT)
    MH_nBlocksY           =  INT(nBlocksY,       kind=MILHOJA_INT)
    MH_nBlocksZ           =  INT(nBlocksZ,       kind=MILHOJA_INT)
    MH_maxRefinementLevel =  INT(lRefineMax,     kind=MILHOJA_INT)
    MH_nGuard             =  INT(NGUARD,         kind=MILHOJA_INT)
    MH_nCcVars            =  INT(NUNK_VARS,      kind=MILHOJA_INT)
    MH_nFluxVars          =  INT(NFLUXES,        kind=MILHOJA_INT)
    MH_ccInterpolator     =  INT(ccInterpolator, kind=MILHOJA_INT)

    select case (gr_geometry)
    case (CARTESIAN)
        MH_coordSys = MILHOJA_CARTESIAN
    case (CYLINDRICAL)
        MH_coordSys = MILHOJA_CYLINDRICAL
    case (SPHERICAL)
        MH_coordSys = MILHOJA_SPHERICAL
    case default
        CALL Driver_abort("[Grid_init] Unsupported coord sys")
    end select

    MH_loBCs(:) = INT(MILHOJA_PERIODIC, kind=MILHOJA_INT)
    MH_hiBCs(:) = INT(MILHOJA_PERIODIC, kind=MILHOJA_INT)
    do i = 1, NDIM
        if(gr_domainBC(LOW,  i) /= PERIODIC) then
            MH_loBCs(i) = INT(MILHOJA_EXTERNAL_BC, kind=MILHOJA_INT)
        end if
        if(gr_domainBC(HIGH, i) /= PERIODIC) then
            MH_hiBCs(i) = INT(MILHOJA_EXTERNAL_BC, kind=MILHOJA_INT)
        end if
    end do

    CALL milhoja_grid_init(gr_globalComm, MH_logRank,             &
                           MH_coordSys,                           &
                           MH_xMin, MH_xMax,                      &
                           MH_yMin, MH_yMax,                      &
                           MH_zMin, MH_zMax,                      &
                           MH_loBCs, MH_hiBCs,                    &
                           gr_fillPhysicalBcCallback,             &
                           MH_nxb, MH_nyb, MH_nzb,                &
                           MH_nBlocksX, MH_nBlocksY, MH_nBlocksZ, &
                           MH_maxRefinementLevel,                 &
                           MH_ccInterpolator,                     &
                           MH_nGuard, MH_nCcVars, MH_nFluxVars,   &
                           gr_markRefineDerefineCallback,         &
                           MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)

    ! Sanity check configuration
    ! NOTE: Some sanity checks don't act in the direct interest of Flash-X,
    ! but rather as a unit test for those developing and maintaining the
    ! Milhoja Fortran interface.
    CALL milhoja_grid_getCoordinateSystem(MH_coordSys, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    select case (MH_coordSys)
    case (MILHOJA_CARTESIAN)
        if (gr_geometry /= CARTESIAN) then
            CALL Driver_abort("[Grid_init] Milhoja CoordSys is Cartesian")
        end if
    case (MILHOJA_CYLINDRICAL)
        if (gr_geometry /= CYLINDRICAL) then
            CALL Driver_abort("[Grid_init] Milhoja CoordSys is Cylindrical")
        end if
    case (MILHOJA_SPHERICAL)
        if (gr_geometry /= SPHERICAL) then
            CALL Driver_abort("[Grid_init] Milhoja CoordSys is Spherical")
        end if
    case default
        CALL Driver_abort("[Grid_init] Unknown Milhoja CoordSys")
    end select

    ! Set garbage values to confirm that they are not overwritten above NDIM
    MH_domainLo(:) = -1.23456
    MH_domainHi(:) =  6.54321
    MH_nBlocksX_SC = -1
    MH_nBlocksY_SC = -1
    MH_nBlocksZ_SC = -1
    CALL milhoja_grid_getDomainBoundBox(MH_domainLo, MH_domainHi, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    CALL milhoja_grid_getDomainDecomposition(MH_nBlocksX_SC, &
                                             MH_nBlocksY_SC, &
                                             MH_nBlocksZ_SC, &
                                             MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (     (MH_domainLo(IAXIS) /= MH_xMin) &
        .OR. (MH_domainHi(IAXIS) /= MH_xMax)) then
        CALL Driver_abort("[Grid_init] Bad xMin/xMax config in Milhoja")
    else if (MH_nBlocksX_SC /= MH_nBlocksX) then
        CALL Driver_abort("[Grid_init] Bad nBlocksX config in Milhoja")
    end if
#if NDIM >= 2
    if (     (MH_domainLo(JAXIS) /= MH_yMin) &
        .OR. (MH_domainHi(JAXIS) /= MH_yMax)) then
        CALL Driver_abort("[Grid_init] Bad yMin/yMax config in Milhoja")
    else if (MH_nBlocksY_SC /= MH_nBlocksY) then
        CALL Driver_abort("[Grid_init] Bad nBlocksY config in Milhoja")
    end if
#else
    if (     (MH_domainLo(JAXIS) /= -1.23456) &
        .OR. (MH_domainHi(JAXIS) /=  6.54321)) then
        CALL Driver_abort("[Grid_init] Bad yMin/yMax config in Milhoja")
    else if (MH_nBlocksY_SC /= -1) then
        CALL Driver_abort("[Grid_init] Bad nBlocksY config in Milhoja")
    end if
#endif
#if NDIM == 3
    if (     (MH_domainLo(KAXIS) /= MH_zMin) &
        .OR. (MH_domainHi(KAXIS) /= MH_zMax)) then
        CALL Driver_abort("[Grid_init] Bad zMin/zMax config in Milhoja")
    else if (MH_nBlocksZ_SC /= MH_nBlocksZ) then
        CALL Driver_abort("[Grid_init] Bad nBlocksZ config in Milhoja")
    end if
#else
    if (     (MH_domainLo(KAXIS) /= -1.23456) &
        .OR. (MH_domainHi(KAXIS) /=  6.54321)) then
        CALL Driver_abort("[Grid_init] Bad zMin/zMax config in Milhoja")
    else if (MH_nBlocksZ_SC /= -1) then
        CALL Driver_abort("[Grid_init] Bad nBlocksZ config in Milhoja")
    end if
#endif

    CALL milhoja_grid_getMaxFinestLevel(MH_level, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (MH_level /= MH_maxRefinementLevel) then
        CALL Driver_abort("[Grid_init] Bad max refine level in Milhoja")
    end if

    ! Set garbage values to confirm that they are not overwritten above NDIM
    MH_nxb = -1
    MH_nyb = -1
    MH_nzb = -1
    CALL milhoja_grid_getBlockSize(MH_nxb, MH_nyb, MH_nzb, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (MH_nxb /= NXB) then
        CALL Driver_abort("[Grid_init] Bad nxb config in Milhoja")
    end if
#if NDIM >= 2
    if (MH_nyb /= NYB) then
        CALL Driver_abort("[Grid_init] Bad nyb config in Milhoja")
    end if
#else
    if (MH_nyb /= -1) then
        CALL Driver_abort("[Grid_init] Bad nyb config in Milhoja")
    end if
#endif
#if NDIM == 3
    if (MH_nzb /= NZB) then
        CALL Driver_abort("[Grid_init] Bad nzb config in Milhoja")
    end if
#else
    if (MH_nzb /= -1) then
        CALL Driver_abort("[Grid_init] Bad nzb config in Milhoja")
    end if
#endif

    CALL milhoja_grid_getNGuardcells(MH_nGuard, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (MH_nGuard /= NGUARD) then
        CALL Driver_abort("[Grid_init] Bad nGuard config in Milhoja")
    end if

    CALL milhoja_grid_getNCcVariables(MH_nCcVars, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (MH_nCcVars /= NUNK_VARS) then
        CALL Driver_abort("[Grid_init] Bad nCcVars config in Milhoja")
    end if

    CALL milhoja_grid_getNFluxVariables(MH_nFluxVars, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    if (MH_nFluxVars /= NFLUXES) then
        CALL Driver_abort("[Grid_init] Bad nFluxVars config in Milhoja")
    end if

    !----------------------------------------------------------------------------------
    ! Cache Milhoja-"owned" data as local Grid data variables for optimization
    !----------------------------------------------------------------------------------
    ! The goal is *not* to cache all variable as pre-optimization, but rather to
    ! cache based on a profiling-proven need to cache.  Indeed, it is better to
    ! only have one copy of a configuration value in existence rather than two.
    !
    ! Note that we may be forced to cache some data as more general Grid code
    ! may be written under the assumption that some variables are initialized
    ! by all Grid unit implementations.
    !
    ! Load and cache the values from Milhoja rather than from the runtime
    ! parameters.  This is inline with the notion that Milhoja owns this data.

    ! The following routine does not set values above NDIM.  Therefore, set
    ! the correct value beforehand
    gr_globalDomain(:, :) = 0.0
    CALL milhoja_grid_getDomainBoundBox(MH_domainLo, MH_domainHi, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    !The following assignments should be unnecessary -- basically no-ops --
    !since gr_globalDomain is already being set in gr_initGeometry,
    !and those better be the same values!
    do i = 1, NDIM
        gr_globalDomain(LOW,  i) = REAL(MH_domainLo(i))
        gr_globalDomain(HIGH, i) = REAL(MH_domainHi(i))
    end do

    CALL milhoja_grid_getMaxFinestLevel(MH_level, MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)
    gr_lRefineMax = INT(MH_level)

    gr_maxRefine = gr_lRefineMax

    CALL Grid_getDeltas(gr_lRefineMax, gr_minCellSizes)
    gr_minCellSize = gr_minCellSizes(IAXIS)
#if NDIM >= 2
    if (.NOT. gr_dirIsAngular(JAXIS)) then
        gr_minCellSize = MIN(gr_minCellSize, gr_minCellSizes(JAXIS))
    end if
#endif
#if NDIM == 3
    if (.NOT. gr_dirIsAngular(KAXIS)) then
        gr_minCellSize = MIN(gr_minCellSize, gr_minCellSizes(KAXIS))
    end if
#endif

    gr_allPeriodic = .TRUE.
    do i = 1, NDIM
        if (gr_domainBC(LOW,  i) /= PERIODIC) then
            gr_allPeriodic = .FALSE.
        end if
        if (gr_domainBC(HIGH, i) /= PERIODIC) then
            gr_allPeriodic = .FALSE.
        end if
    end do

    !----------------------------------------------------------------------------------
    ! Setup all remaining local Grid data variables
    !----------------------------------------------------------------------------------
    do i = UNK_VARS_BEGIN, UNK_VARS_END
        CALL Simulation_getVarnameType(i, gr_vartypes(i))
    end do

    CALL RuntimeParameters_get("gr_enableTiling", gr_enableTiling)
    if (gr_enableTiling) then
        CALL Driver_abort("[Grid_init] Tiling not implemented yet")
    end if
    CALL RuntimeParameters_get("gr_useTiling",    gr_useTiling)
    if (gr_useTiling .AND. (.NOT. gr_enableTiling)) then
       if (gr_meshMe == MASTER_PE) then
          write(*,*) "[Grid_init] WARNING: Tiling is disabled by gr_enableTiling=.FALSE."
          write(*,*) "                     Therefore gr_useTiling will not take effect."
       end if
    end if
    CALL RuntimeParameters_get("gr_tileSizeX", gr_tileSize(IAXIS))
    CALL RuntimeParameters_get("gr_tileSizeY", gr_tileSize(JAXIS))
    CALL RuntimeParameters_get("gr_tileSizeZ", gr_tileSize(KAXIS))

    ! DEV: Unable to set an upper limit on these runtime parameters using
    ! N[XYZ]B.  Therefore, we are checking here.
    if (gr_tileSize(IAXIS) > NXB) then
        if (gr_meshMe == MASTER_PE) then
           print*,'WARNING: Tile size along x-axis cannot exceed block size along x'
           print*,'         Setting tile size to block size in x'
        end if
        gr_tileSize(IAXIS) = NXB
    end if
    if (gr_tileSize(JAXIS) > NYB) then
        if (gr_meshMe == MASTER_PE) then
           print*,'WARNING: Tile size along y-axis cannot exceed block size along y'
           print*,'         Setting tile size to block size in y'
        end if
        gr_tileSize(JAXIS) = NYB
    end if
    if (gr_tileSize(KAXIS) > NZB) then
        if (gr_meshMe == MASTER_PE) then
           print*,'WARNING: Tile size along z-axis cannot exceed block size along z'
           print*,'         Setting tile size to block size in z'
        end if
        gr_tileSize(KAXIS) = NZB
    end if

    CALL RuntimeParameters_get("smallx", gr_smallx)
    CALL RuntimeParameters_get("smalle", gr_smalle)
    CALL RuntimeParameters_get("smlrho", gr_smallrho)

    CALL RuntimeParameters_get("eosMode", eosModeString)
    CALL RuntimeParameters_mapStrToInt(eosModeString, gr_eosMode)

    CALL RuntimeParameters_get("eosModeInit", eosModeString)
    CALL RuntimeParameters_mapStrToInt(eosModeString, gr_eosModeInit)

    call RuntimeParameters_get('useOrchestration',gr_useOrchestration)

#ifdef GRID_WITH_MONOTONIC
    CALL Driver_abort("[Grid_init] Monotonic not implemented with Milhoja yet")
    ! Could possibly be less if gr_intpol < 2  - KW
    !gr_intpolStencilWidth = 2     
#endif

    ! This section of the code identifies the variables to use in
    ! the refinement criterion. If a variable is a refinement variable
    ! then the corresponding refinement/derefinement cutoff and filter
    ! values also have to be fetched. The config file defines
    ! refinement variables as strings, names as "refine_var_1",
    ! "refine_var_2" etc, with the current maximum being 4. The
    ! general utility routine takes the base "refine_var_" and appends
    ! the index at the end of the string to generate the parameter
    ! name and the routine Simulation_mapStrToInt finds its index into UNK.
    CALL RuntimeParameters_get("refine_var_count", gr_numRefineVarsMax)
    gr_refine_var = NONEXISTENT
    gr_numRefineVars = 0

    do i = 1, gr_numRefineVarsMax
        CALL concatStringWithInt(REF_VAR_NAME, i, refVarString)
        CALL RuntimeParameters_get(refVarString, paramString)
        if (paramString /= "none") then
            CALL Simulation_mapStrToInt(paramString, refVar, MAPBLOCK_UNK)
            if (refVar > 0) then
                nonrep = 0
                ! DEV: FIXME This is segfaulting, but isn't needed right now
!                CALL Grid_getVarNonRep(MAPBLOCK_UNK, refVar, nonrep)
                if (nonrep > 0) then
                    refVar = 0
                else
                    gr_numRefineVars = gr_numRefineVars + 1
                    gr_refine_var(gr_numRefineVars) = refVar
                    CALL concatStringWithInt(REF_CUTOFF_NAME, gr_numRefineVars, &
                                             refCutoffString)
                    CALL concatStringWithInt(DEREF_CUTOFF_NAME, gr_numRefineVars, &
                                             derefCutOffString)
                    CALL concatStringWithInt(REF_FILTER_NAME, gr_numRefineVars, &
                                             refFilterString)
                    CALL RuntimeParameters_get(refCutoffString, &
                                               gr_refine_cutoff(gr_numRefineVars))
                    CALL RuntimeParameters_get(derefCutoffString, &
                                               gr_derefine_cutoff(gr_numRefineVars))
                    CALL RuntimeParameters_get(refFilterString, &
                                               gr_refine_filter(gr_numRefineVars))
                end if
            end if

            if (refVar <= 0) then
                if (gr_globalMe == MASTER_PE) then
                   write(*,*) 'WARNING: Unrecognized or non-replicated variable name in ', &
                              REF_VAR_NAME, i,' treating it as "none"'
                end if
                CALL Logfile_stampMessage( &
                   'WARNING: Unrecognized or non-replicated variable name in refine_var, treating it as "none"')
            end if
        end if
    end do

    if (gr_numRefineVars == 0) then
        if (gr_meshMe == MASTER_PE) then
            write(*,*) 'WARNING : Adaptive Grid did not find any refinement variables'
        end if
        CALL Logfile_stampMessage("WARNING : Adaptive Grid did not find any variable to refine")
    end if

    gr_enforceMaxRefinement = .FALSE.

    gr_justExchangedGC = .FALSE.
    gr_gcellsUpToDate  = .FALSE.

    !----------------------------------------------------------------------------------
    ! Log Grid setup information to stdout
    !----------------------------------------------------------------------------------
    if (gr_globalMe == MASTER_PE) then
        write(*,'(A)')                 "[Grid] Global Domain"
        write(*,'(A,F15.7,A,F15.7,A)') "       X in (", gr_globalDomain(LOW,  IAXIS), ", ", &
                                                        gr_globalDomain(HIGH, IAXIS), ")"
        write(*,'(A,F15.7,A,F15.7,A)') "       Y in (", gr_globalDomain(LOW,  JAXIS), ", ", &
                                                        gr_globalDomain(HIGH, JAXIS), ")"
        write(*,'(A,F15.7,A,F15.7,A)') "       Z in (", gr_globalDomain(LOW,  KAXIS), ", ", &
                                                        gr_globalDomain(HIGH, KAXIS), ")"

        write(*,'(A,I0)') "[Grid] lRefineMax = ", gr_lRefineMax

        write(*,'(A)') '[Grid] Resolution based on runtime parameters'
        write(*,'(A9,3(A12:4x))')  'lrefine', (COL_HDR(i), i=1,MDIM)
        do level = 1, gr_lRefineMax
            CALL Grid_getDeltas(level, deltas)
            if (maxval(deltas(IAXIS:NDIM)) .GT. 999999999999.999) then
                fmtStr = '(I7,2x,1P,3G16.3)'
            else if (minval(deltas(IAXIS:NDIM)) .LE. 0.0009) then
                fmtStr = '(I7,2x,1P,3G16.3)'
            else
                fmtStr = '(I7,2x,3F16.3)'
            end if

            write(*,fmtStr) level, (deltas(i), i=1,MDIM)
        end do
    end if

    !----------------------------------------------------------------------------------
    ! Triumphantly notify the world that all is done!
    !----------------------------------------------------------------------------------
    if (gr_globalMe == MASTER_PE) then
        write(*,'(A)') "[Grid] Initialized and ready for use"
    end if
end subroutine Grid_init

