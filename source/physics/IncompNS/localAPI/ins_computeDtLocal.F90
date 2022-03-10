!!****if* source/physics/localAPI/ins_computeDtLocal
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
!!  ins_computeDtLocal
!!
!!
!! SYNOPSIS
!!
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***
subroutine ins_computeDtLocal(blockID, &
                              isize, jsize, ksize, &
                              dx, dy, dz, &
                              blkLimits, blkLimitsGC, &
                              facexData, faceyData, &
                              facezData, &
                              dtLocal, lminloc)
   implicit none
   integer, intent(IN) :: blockID
   integer, dimension(2, MDIM), intent(IN) :: blkLimits, blkLimitsGC
   integer, intent(IN) :: isize, jsize, ksize
   real, intent(IN) :: dx, dy, dz
   real, pointer, dimension(:, :, :, :)  :: facexData, faceyData, facezData
   real, intent(INOUT) :: dtLocal
   integer, intent(INOUT) :: lminloc(5)
end subroutine ins_computeDtLocal
