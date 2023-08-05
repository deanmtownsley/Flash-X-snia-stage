!!****if* source/physics/ImBound/ImBoundMain/ImBound_advance
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
!! NAME
!!
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the ImBound module.
!!
!!***

#include "constants.h"

subroutine ImBound_advance(bodyInfo, time, dt)

   use ImBound_data, ONLY: ib_bruteForceMapping
   use ImBound_type, ONLY: ImBound_type_t
   use ib_interface, ONLY: ib_annBuildTree

   implicit none
   type(ImBound_type_t), intent(inout) :: bodyInfo
   real, intent(in) :: time, dt
   bodyInfo%elems(:)%xA = bodyInfo%elems(:)%xA + dt*bodyInfo%velx
   bodyInfo%elems(:)%yA = bodyInfo%elems(:)%yA + dt*bodyInfo%vely

   bodyInfo%elems(:)%xB = bodyInfo%elems(:)%xB + dt*bodyInfo%velx
   bodyInfo%elems(:)%yB = bodyInfo%elems(:)%yB + dt*bodyInfo%vely

   bodyInfo%elems(:)%xCenter = bodyInfo%elems(:)%xCenter + dt*bodyInfo%velx
   bodyInfo%elems(:)%yCenter = bodyInfo%elems(:)%yCenter + dt*bodyInfo%vely

   !bodyInfo%elems(:)%xCenter = (bodyInfo%elems(:)%xA + bodyInfo%elems(:)%xB)/2
   !bodyInfo%elems(:)%yCenter = (bodyInfo%elems(:)%yA + bodyInfo%elems(:)%yB)/2

   bodyInfo%boundBox(:, IAXIS) = (/minval(bodyInfo%elems(:)%xCenter), &
                                   maxval(bodyInfo%elems(:)%xCenter)/)

   bodyInfo%boundBox(:, JAXIS) = (/minval(bodyInfo%elems(:)%yCenter), &
                                   maxval(bodyInfo%elems(:)%yCenter)/)

   if (.not. ib_bruteForceMapping) call ib_annBuildTree(bodyInfo)

end subroutine ImBound_advance
