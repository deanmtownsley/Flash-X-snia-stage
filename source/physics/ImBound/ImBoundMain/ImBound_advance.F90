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
#include "Simulation.h"

subroutine ImBound_advance(bodyInfo, time, dt)

   use ImBound_data, ONLY: ib_bruteForceMapping
   use ImBound_type, ONLY: ImBound_type_t
   use ib_interface, ONLY: ib_annBuildTree

   implicit none
   type(ImBound_type_t), intent(inout) :: bodyInfo

   real, intent(in) :: time, dt
   real, dimension(2, 2) :: rotate
   real, dimension(2) :: offset
   real, dimension(2) :: vector
   integer :: panelIndex

   rotate(:, 1) = (/cos(dt*bodyInfo%thetaz), -sin(dt*bodyInfo%thetaz)/)
   rotate(:, 2) = (/sin(dt*bodyInfo%thetaz), cos(dt*bodyInfo%thetaz)/)

   offset(:) = (/dt*bodyInfo%velx, dt*bodyInfo%vely/)

   do panelIndex = 1, bodyInfo%numElems
      vector(1) = bodyInfo%elems(panelIndex)%xA
      vector(2) = bodyInfo%elems(panelIndex)%yA
      vector = matmul(vector, rotate) + offset
      bodyInfo%elems(panelIndex)%xA = vector(1)
      bodyInfo%elems(panelIndex)%yA = vector(2)

      vector(1) = bodyInfo%elems(panelIndex)%xB
      vector(2) = bodyInfo%elems(panelIndex)%yB
      vector = matmul(vector, rotate) + offset
      bodyInfo%elems(panelIndex)%xB = vector(1)
      bodyInfo%elems(panelIndex)%yB = vector(2)

      vector(1) = bodyInfo%elems(panelIndex)%xCenter
      vector(2) = bodyInfo%elems(panelIndex)%yCenter
      vector = matmul(vector, rotate) + offset
      bodyInfo%elems(panelIndex)%xCenter = vector(1)
      bodyInfo%elems(panelIndex)%yCenter = vector(2)

      vector(1) = bodyInfo%elems(panelIndex)%xNorm
      vector(2) = bodyInfo%elems(panelIndex)%yNorm
      vector = matmul(vector, rotate)
      bodyInfo%elems(panelIndex)%xNorm = vector(1)
      bodyInfo%elems(panelIndex)%yNorm = vector(2)
   end do

   bodyInfo%boundBox(:, IAXIS) = (/minval(bodyInfo%elems(:)%xCenter), &
                                   maxval(bodyInfo%elems(:)%xCenter)/)

   bodyInfo%boundBox(:, JAXIS) = (/minval(bodyInfo%elems(:)%yCenter), &
                                   maxval(bodyInfo%elems(:)%yCenter)/)

   call ib_annBuildTree(bodyInfo)

end subroutine ImBound_advance
