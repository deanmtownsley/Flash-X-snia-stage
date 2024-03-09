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

#include "Simulation.h"

#include "Milhoja.h"

#ifdef ORCHESTRATION_USE_DATAPACKETS
#ifndef RUNTIME_SUPPORT_DATAPACKETS
#error "Milhoja library does not have support for datapackets"
#endif
#endif

#ifdef ORCHESTRATION_USE_GPUS
#ifndef MILHOJA_GPUS_SUPPORTED
#error "Milhoja library does not have GPU-capable runtime backend"
#endif
#endif

!> @ingroup OrchestrationMilhoja
!! @stubref{Orchestration_init}
!!
!! @brief Concrete implementation of Orchestration_init
!!
!! @todo Should this confirm matching types between Milhoja and Flash-X?
!! @todo Error check that nBytesInMemoryPools cast doesn't overflow
subroutine Orchestration_init()
    use milhoja_types_mod,           ONLY : MILHOJA_INT, &
                                            MILHOJA_SIZE_T
    use milhoja_runtime_mod,         ONLY : milhoja_runtime_init

    use Orchestration_interface,     ONLY : Orchestration_checkInternalError
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get

    implicit none

    integer :: nThreadTeams
    integer :: nThreadsPerTeam
    integer :: nStreams
    real    :: nBytesInCpuMemoryPool
    real    :: nBytesInGpuMemoryPools

    integer(MILHOJA_INT)    :: MH_nThreadTeams
    integer(MILHOJA_INT)    :: MH_nThreadsPerTeam
    integer(MILHOJA_INT)    :: MH_nStreams
    integer(MILHOJA_SIZE_T) :: MH_nBytesInCpuMemoryPool
    integer(MILHOJA_SIZE_T) :: MH_nBytesInGpuMemoryPools
    integer(MILHOJA_INT)    :: MH_ierr

#ifndef FULL_MILHOJAGRID
    call fake_milhoja_grid_init
#endif
    !!!!!----- Runtime Parameters
    CALL RuntimeParameters_get("or_nThreadTeams",           nThreadTeams)
    CALL RuntimeParameters_get("or_nThreadsPerTeam",        nThreadsPerTeam)
    CALL RuntimeParameters_get("or_nBytesInCpuMemoryPool",  nBytesInCpuMemoryPool)
#ifdef ORCHESTRATION_USE_GPUS
    CALL RuntimeParameters_get("or_nStreams",               nStreams)
    CALL RuntimeParameters_get("or_nBytesInGpuMemoryPools", nBytesInGpuMemoryPools)
#else
    !!DEV: Maybe use some of those RPs when ORCHESTRATION_USE_PACKETS is defined? - KW
    nStreams = 0
    nBytesInGpuMemoryPools = 0.0
#endif

    !!!!!----- CAST TO MILHOJA TYPES
    MH_nThreadTeams    = INT(nThreadTeams,    kind=MILHOJA_INT)
    MH_nThreadsPerTeam = INT(nThreadsPerTeam, kind=MILHOJA_INT)
    MH_nStreams        = INT(nStreams,        kind=MILHOJA_INT)

    ! The runtime parameters unit does not allow for declaring
    ! or_nBytesInMemoryPools as the equivalent of size_t or long long.  However,
    ! this parameter will likely need values that are too large for 32-bit.
    ! Therefore, as a work around, this parameter is specified as a real.
    !
    ! Cast it to the type required by the runtime.
    MH_nBytesInCpuMemoryPool  = NINT(nBytesInCpuMemoryPool,  kind=MILHOJA_SIZE_T)
    MH_nBytesInGpuMemoryPools = NINT(nBytesInGpuMemoryPools, kind=MILHOJA_SIZE_T)

    !!!!!----- Initialize library
    CALL milhoja_runtime_init(MH_nThreadTeams,           &
                              MH_nThreadsPerTeam,        &
                              MH_nStreams,               &
                              MH_nBytesInCpuMemoryPool,  &
                              MH_nBytesInGpuMemoryPools, &
                              MH_ierr)
    CALL Orchestration_checkInternalError("Orchestration_init", MH_ierr)

  contains
    subroutine fake_milhoja_grid_init()
#if(0)
      ! initially copied from GridMain/AMR/Milhoja/Grid_init.F90
    use milhoja_types_mod,           ONLY : MILHOJA_REAL
#ifdef FULL_MILHOJAGRID
    use milhoja_grid_mod,            ONLY : milhoja_grid_init, &
                                            milhoja_grid_getCoordinateSystem, &
                                            milhoja_grid_getDomainBoundBox, &
                                            milhoja_grid_getMaxFinestLevel, &
                                            milhoja_grid_getBlockSize, &
                                            milhoja_grid_getDomainDecomposition, &
                                            milhoja_grid_getNGuardcells, &
                                            milhoja_grid_getNCcVariables, &
                                            milhoja_grid_getNFluxVariables
#else
    use milhoja_grid_mod,            ONLY : milhoja_grid_init!, &
!!$                                            milhoja_grid_getCoordinateSystem, &
!!$                                            milhoja_grid_getDomainBoundBox, &
!!$                                            milhoja_grid_getMaxFinestLevel, &
!!$                                            milhoja_grid_getBlockSize, &
!!$                                            milhoja_grid_getDomainDecomposition, &
!!$                                            milhoja_grid_getNGuardcells, &
!!$                                            milhoja_grid_getNCcVariables, &
!!$                                            milhoja_grid_getNFluxVariables
#endif
    
    use Grid_data,                ONLY: gr_globalComm, &
                                        gr_geometry
    use Orchestration_data,          ONLY : or_globalComm, &
                                            or_globalMe, &
                                            or_meshComm, &
                                            or_meshMe, &
                                            or_meshNumProcs
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
#include "constants.h"

    integer :: nBlocksX, nBlocksY, nBlocksZ
    integer :: lRefineMax

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
    
    !!!!!----- MPI & OpenMP Environment
    CALL Driver_getMype(GLOBAL_COMM,     or_globalMe)
!    CALL Driver_getNumProcs(GLOBAL_COMM, or_globalNumProcs)
    CALL Driver_getComm(GLOBAL_COMM,     or_globalComm)

    CALL Driver_getMype(MESH_COMM,     or_meshMe)
    CALL Driver_getNumProcs(MESH_COMM, or_meshNumProcs)
    CALL Driver_getComm(MESH_COMM,     or_meshComm)

    !!!!!----- Runtime Parameters
    if (or_globalMe == MASTER_PE) then
        write(*,'(A)') "[FAKE MILHOJA GRID] Initializing ..."
    end if

    CALL RuntimeParameters_get('lrefine_max', lRefineMax)

    CALL RuntimeParameters_get("nblockx", nBlocksX)
    CALL RuntimeParameters_get("nblocky", nBlocksY)
    CALL RuntimeParameters_get("nblockz", nBlocksZ)

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
    MH_ccInterpolator     =  INT(UNKNOWN,        kind=MILHOJA_INT)

    CALL milhoja_grid_init(gr_globalComm, MH_logRank,             &
                           MH_coordSys,                           &
                           MH_xMin, MH_xMax,                      &
                           MH_yMin, MH_yMax,                      &
                           MH_zMin, MH_zMax,                      &
                           MH_loBCs, MH_hiBCs,                    &
#ifdef FULL_MILHOJAGRID
                           gr_fillPhysicalBcCallback,             &
#endif
                           MH_nxb, MH_nyb, MH_nzb,                &
                           MH_nBlocksX, MH_nBlocksY, MH_nBlocksZ, &
                           MH_maxRefinementLevel,                 &
                           MH_ccInterpolator,                     &
                           MH_nGuard, MH_nCcVars, MH_nFluxVars,   &
#ifdef FULL_MILHOJAGRID
                           gr_markRefineDerefineCallback,         &
#endif
                           MH_ierr)
    CALL gr_checkMilhojaError("Grid_init", MH_ierr)

    ! Sanity check configuration
    ! NOTE: Some sanity checks don't act in the direct interest of Flash-X,
    ! but rather as a unit test for those developing and maintaining the
    ! Milhoja Fortran interface.
#ifdef FULL_MILHOJAGRID
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
#endif
    if (MH_nFluxVars /= NFLUXES) then
        CALL Driver_abort("[Grid_init] Bad nFluxVars config in Milhoja")
    end if
#endif
  end subroutine fake_milhoja_grid_init
end subroutine Orchestration_init

