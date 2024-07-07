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
    use Orchestration_data
    use milhoja_types_mod,           ONLY : MILHOJA_INT, &
                                            MILHOJA_SIZE_T
    use milhoja_runtime_mod,         ONLY : milhoja_runtime_init

    use Orchestration_interface,     ONLY : Orchestration_checkInternalError
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get

#include "constants.h"
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

    !!!!!----- Runtime Parameters
    CALL RuntimeParameters_get("or_nThreadTeams",           nThreadTeams)
    CALL RuntimeParameters_get("or_nThreadsPerTeam",        nThreadsPerTeam)
    CALL RuntimeParameters_get("or_nBytesInCpuMemoryPool",  nBytesInCpuMemoryPool)
#if defined(ORCHESTRATION_USE_GPUS) || defined(ORCHESTRATION_USE_DATAPACKETS)
    CALL RuntimeParameters_get("or_nStreams",               nStreams)
    CALL RuntimeParameters_get("or_nBytesInGpuMemoryPools", nBytesInGpuMemoryPools)
#else
    !!DEV: OK to use those RPs whenever ORCHESTRATION_USE_DATAPACKETS is defined? - KW
    nStreams = 0
    nBytesInGpuMemoryPools = 0.0
#endif
    CALL RuntimeParameters_get("or_nThreads_1",             or_nThreads_1)
    CALL RuntimeParameters_get("or_nThreads_2",             or_nThreads_2)
    CALL RuntimeParameters_get("or_nThreads_3",             or_nThreads_3)
    CALL RuntimeParameters_get("or_nThreads_4",             or_nThreads_4)
    CALL RuntimeParameters_get("or_nThreads_5",             or_nThreads_5)

    CALL RuntimeParameters_get("or_nTilesPerPacket_1",      or_nTilesPerPacket_1)
    CALL RuntimeParameters_get("or_nTilesPerPacket_2",      or_nTilesPerPacket_2)
    CALL RuntimeParameters_get("or_nTilesPerPacket_3",      or_nTilesPerPacket_3)
    CALL RuntimeParameters_get("or_nTilesPerPacket_4",      or_nTilesPerPacket_4)
    CALL RuntimeParameters_get("or_nTilesPerPacket_5",      or_nTilesPerPacket_5)

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

end subroutine Orchestration_init

