!!****if* source/physics/localAPI/ins_indicators
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
!!*****
subroutine ins_indicators(u, v, w, pres, divv, ix1, ix2, jy1, jy2, kz1, kz2, vecminaux, vecmaxaux)
   implicit none
   real, dimension(:, :, :), intent(in) :: u, v, w
   real, dimension(:, :, :), intent(in) :: pres, divv
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, dimension(5), intent(inout) :: vecminaux, vecmaxaux
end subroutine ins_indicators
