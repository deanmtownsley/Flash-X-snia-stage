!!***if* source/Simulation/localAPI/sim_heaterCheckSitesBlk
!!
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
!!***

subroutine sim_heaterCheckSitesBlk2d(phi, xcell, ycell, boundBox, ix1, ix2, jy1, jy2)
   implicit none
   real, dimension(:, :, :), intent(in)      :: phi
   real, dimension(:), intent(in)          :: xcell, ycell
   real, dimension(:, :), intent(in)        :: boundBox
   integer, intent(in)                    :: ix1, ix2, jy1, jy2
end subroutine sim_heaterCheckSitesBlk2d

subroutine sim_heaterCheckSitesBlk3d(phi, xcell, ycell, zcell, boundBox, ix1, ix2, jy1, jy2, kz1, kz2)
   implicit none
   real, dimension(:, :, :), intent(in)  :: phi
   real, dimension(:), intent(in)      :: xcell, ycell, zcell
   real, dimension(:, :), intent(in)    :: boundBox
   integer, intent(in)                :: ix1, ix2, jy1, jy2, kz1, kz2
end subroutine sim_heaterCheckSitesBlk3d
