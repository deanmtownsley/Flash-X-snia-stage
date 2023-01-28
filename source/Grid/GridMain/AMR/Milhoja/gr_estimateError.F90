#include "Simulation.h"

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
!! For each block, estimate the error associated with the given variable to
!! help determine if the block needs refinement or derefinement.  Update the
!! corresponding value in the error array to be the maximum of the incoming
!! value and the value calculated here.
!!
!! This subroutine has not been implemented and aborts if called.
!!
!! @todo Paramesh and AMReX also have versions of this, but there is no stub
!! nor explicit interface.  This implies that all three versions can be
!! distinct.  Is this as intended?  If not, then do we need at least a stub that
!! defines the generic interface of this subroutine and has more detailed
!! official documentation?  Note that it is called from
!!             GridMain/AMR/Grid_markRefineDerefine.F90.
!! If this usage is correct, then this implies the need for a generic interface.
!! Should that Grid_markRefineDerefine be moved to Paramesh and at the same time
!! make gr_estimateError local to Paramesh?
!!
!! @param error          indexed by block IDs
!! @param iref           index of the refinement variable in data structure
!!                       "unk"
!! @param refine_filter  makes sure that error calculations to determine
!!                       refinement don't diverge numerically 
subroutine gr_estimateError(error, iref, refine_filter)
  use Driver_interface, ONLY : Driver_abort

  implicit none

  real,    intent(INOUT) :: error(MAXBLOCKS)
  integer, intent(IN)    :: iref
  real,    intent(IN)    :: refine_filter
 
  call Driver_abort("[gr_estimateError] Not implemented with Milhoja")
end subroutine gr_estimateError

