!!****if* source/numericalTools/MoL/MoLMemory/ml_memAlloc
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!      ml_memAlloc
!!
!!  SYNOPSIS
!!
!!      call ml_memAlloc()
!!
!!  DESCRIPTION
!!
!!      Allocate memory at startup (after variable registration)
!!
!! ARGUMENTS
!!
!!
!!***
subroutine ml_memAlloc()

    use ml_memInterface, only: ml_memFree
    use ml_memData,      only: scratch_data
    use MoL_data,        only: MoL_nscratch_total
    use MoL_variables,   only: MoL_nvars

    use Grid_interface, only: Grid_getBlkIndexLimits

#include "Simulation.h"
#include "constants.h"

    implicit none

    integer, dimension(LOW:HIGH,MDIM) :: lim, limGC

    call Grid_getBlkIndexLimits(1, lim, limGC)

    if (allocated(scratch_data)) call ml_memFree()

    ! Don't need guard-cells here
    allocate(scratch_data(MoL_nvars,                      &
                          lim(LOW,IAXIS):lim(HIGH,IAXIS), &
                          lim(LOW,JAXIS):lim(HIGH,JAXIS), &
                          lim(LOW,KAXIS):lim(HIGH,KAXIS), &
                          MAXBLOCKS,                      &
                          MoL_nscratch_total))

    scratch_data = 0d0
end subroutine ml_memAlloc
