!!***if* source/physics/sourceTerms/Heater/localAPI/htr_checkSitesBlk
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!***

#include "constants.h"
#include "Simulation.h"

subroutine htr_checkSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2, lblock)
   implicit none
   real, dimension(:, :, :), intent(in)      :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)        :: boundBox
   integer, intent(in)                    :: ix1, ix2, jy1, jy2, lblock
end subroutine htr_checkSitesBlk2d

subroutine htr_checkSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2, lblock)
   implicit none
   real, dimension(:, :, :), intent(in)  :: phi
   real, dimension(:), intent(in)      :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)    :: boundBox
   integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2, lblock
end subroutine htr_checkSitesBlk3d
