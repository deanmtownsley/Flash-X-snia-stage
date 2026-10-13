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
   real :: threshold

   call tileDesc%boundBox(boundBox)

   threshold = 0.5*sqrt((bodyInfo%boundBox(HIGH, IAXIS) - bodyInfo%boundBox(LOW, IAXIS))**2 + &
                        (bodyInfo%boundBox(HIGH, JAXIS) - bodyInfo%boundBox(LOW, JAXIS))**2)

   boxOutside = &
      boundBox(HIGH, IAXIS) + threshold .le. bodyInfo%boundBox(LOW, IAXIS) .or. &
      boundBox(HIGH, JAXIS) + threshold .le. bodyInfo%boundBox(LOW, JAXIS) .or. &
      boundBox(LOW, IAXIS) - threshold .ge. bodyInfo%boundBox(HIGH, IAXIS) .or. &
      boundBox(LOW, JAXIS) - threshold .ge. bodyInfo%boundBox(HIGH, JAXIS)

   boxInside = &
      boundBox(HIGH, IAXIS) + threshold .gt. bodyInfo%boundBox(LOW, IAXIS) .and. &
      boundBox(LOW, IAXIS) - threshold .lt. bodyInfo%boundBox(HIGH, IAXIS) .and. &
      boundBox(LOW, IAXIS) - threshold .gt. bodyInfo%boundBox(LOW, IAXIS) .and. &
      boundBox(HIGH, IAXIS) + threshold .lt. bodyInfo%boundBox(HIGH, IAXIS) .and. &
      boundBox(HIGH, JAXIS) + threshold .gt. bodyInfo%boundBox(LOW, JAXIS) .and. &
      boundBox(LOW, JAXIS) - threshold .lt. bodyInfo%boundBox(HIGH, JAXIS) .and. &
      boundBox(LOW, JAXIS) - threshold .gt. bodyInfo%boundBox(LOW, JAXIS) .and. &
      boundBox(HIGH, JAXIS) + threshold .lt. bodyInfo%boundBox(HIGH, JAXIS)

   skipBox = (boxInside .or. boxOutside) .and. ib_enableSelectiveMapping

end subroutine ImBound_skipBox
