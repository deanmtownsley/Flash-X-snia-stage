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

#if NDIM==MDIM

#else
   rotate(:, 1) = (/cos(dt*bodyInfo%theta(3)), -sin(dt*bodyInfo%theta(3))/)
   rotate(:, 2) = (/sin(dt*bodyInfo%theta(3)), cos(dt*bodyInfo%theta(3))/)

   offset(:) = (/dt*bodyInfo%velc(1), dt*bodyInfo%velc(2)/)

   do panelIndex = 1, bodyInfo%numElems

      bodyInfo%elems(panelIndex)%pA(1:2) = &
         matmul(bodyInfo%elems(panelIndex)%pA(1:2), rotate) + offset

      bodyInfo%elems(panelIndex)%pB(1:2) = &
         matmul(bodyInfo%elems(panelIndex)%pB(1:2), rotate) + offset

      bodyInfo%elems(panelIndex)%center(1:2) = &
         matmul(bodyInfo%elems(panelIndex)%center(1:2), rotate) + offset

      bodyInfo%elems(panelIndex)%normal(1:2) = &
         matmul(bodyInfo%elems(panelIndex)%normal(1:2), rotate)

   end do
#endif

   bodyInfo%boundBox(:, IAXIS) = (/minval(bodyInfo%elems(:)%center(1)), &
                                   maxval(bodyInfo%elems(:)%center(1))/)

   bodyInfo%boundBox(:, JAXIS) = (/minval(bodyInfo%elems(:)%center(2)), &
                                   maxval(bodyInfo%elems(:)%center(2))/)
   call ib_annBuildTree(bodyInfo)

end subroutine ImBound_advance
