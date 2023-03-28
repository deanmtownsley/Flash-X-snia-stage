!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief ml_memAlloc implementation

!> @ingroup MoLMemory
!!
!! @brief Implements ml_memAlloc
!!
!! @note This implementation is specific to UG and Paramesh and allocates
!!       a blocklist of the same length and blocksize as found in Grid
!!
!! @stubref{ml_memAlloc}
subroutine ml_memAlloc()

   use ml_memInterface, only: ml_memFree
   use ml_memData, only: ml_scratch_data
   use MoL_data, only: ml_nscratch_total
   use ml_variables, only: ml_nvars

   use Grid_interface, only: Grid_getBlkIndexLimits

#include "Simulation.h"
#include "constants.h"

   implicit none

   integer, dimension(LOW:HIGH, MDIM) :: lim, limGC

   call Grid_getBlkIndexLimits(1, lim, limGC)

   if (allocated(ml_scratch_data)) call ml_memFree()

   ! Don't need guard-cells here
#ifdef MOL_REORDER
   allocate (ml_scratch_data(lim(LOW, IAXIS):lim(HIGH, IAXIS), &
                          lim(LOW, JAXIS):lim(HIGH, JAXIS), &
                          lim(LOW, KAXIS):lim(HIGH, KAXIS), &
                          ml_nvars, &
                          MAXBLOCKS, &
                          ml_nscratch_total))
#else
   allocate (ml_scratch_data(ml_nvars, &
                             lim(LOW, IAXIS):lim(HIGH, IAXIS), &
                             lim(LOW, JAXIS):lim(HIGH, JAXIS), &
                             lim(LOW, KAXIS):lim(HIGH, KAXIS), &
                             MAXBLOCKS, &
                             ml_nscratch_total))
#endif
   ml_scratch_data = 0.0
end subroutine ml_memAlloc
