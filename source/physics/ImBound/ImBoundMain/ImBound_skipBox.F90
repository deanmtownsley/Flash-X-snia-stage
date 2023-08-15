!!****if* source/physics/ImBound/ImBoundMain/ImBound_skipBox
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
!!
!!
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine ImBound_skipBox(tileDesc, bodyInfo, skipBox)

   use Grid_tile, ONLY: Grid_tile_t
   use ImBound_type, ONLY: ImBound_type_t
   use ImBound_data, ONLY: ib_enableSelectiveMapping

   implicit none

   type(Grid_tile_t), INTENT(IN) :: tileDesc
   type(ImBound_type_t), intent(in) :: bodyInfo
   logical, intent(out) :: skipBox

   logical :: boxInside, boxOutside
   real :: boundBox(LOW:HIGH, MDIM)
   real :: del(MDIM), threshold(MDIM)

   call tileDesc%boundBox(boundBox)
   call tileDesc%deltas(del)

   threshold = (/NXB*del(IAXIS), NYB*del(JAXIS), 0./)

   boxOutside = &
      boundBox(HIGH, IAXIS) + threshold(IAXIS) .le. bodyInfo%boundBox(LOW, IAXIS) .or. &
      boundBox(HIGH, JAXIS) + threshold(JAXIS) .le. bodyInfo%boundBox(LOW, JAXIS) .or. &
      boundBox(LOW, IAXIS) - threshold(IAXIS) .ge. bodyInfo%boundBox(HIGH, IAXIS) .or. &
      boundBox(LOW, JAXIS) - threshold(JAXIS) .ge. bodyInfo%boundBox(HIGH, JAXIS)

   boxInside = &
      boundBox(HIGH, IAXIS) + threshold(IAXIS) .gt. bodyInfo%boundBox(LOW, IAXIS) .and. &
      boundBox(LOW, IAXIS) - threshold(IAXIS) .lt. bodyInfo%boundBox(HIGH, IAXIS) .and. &
      boundBox(LOW, IAXIS) - threshold(IAXIS) .gt. bodyInfo%boundBox(LOW, IAXIS) .and. &
      boundBox(HIGH, IAXIS) + threshold(IAXIS) .lt. bodyInfo%boundBox(HIGH, IAXIS) .and. &
      boundBox(HIGH, JAXIS) + threshold(JAXIS) .gt. bodyInfo%boundBox(LOW, JAXIS) .and. &
      boundBox(LOW, JAXIS) - threshold(JAXIS) .lt. bodyInfo%boundBox(HIGH, JAXIS) .and. &
      boundBox(LOW, JAXIS) - threshold(JAXIS) .gt. bodyInfo%boundBox(LOW, JAXIS) .and. &
      boundBox(HIGH, JAXIS) + threshold(JAXIS) .lt. bodyInfo%boundBox(HIGH, JAXIS)

   skipBox = (boxInside .or. boxOutside) .and. ib_enableSelectiveMapping

end subroutine ImBound_skipBox
