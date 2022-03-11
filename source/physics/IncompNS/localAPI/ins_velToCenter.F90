subroutine ins_velToCenter(uf, vf, wf, uc, vc, wc, ix1, ix2, jy1, jy2, kz1, kz2)
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
   implicit none
   real, dimension(:, :, :), intent(in)  :: uf, vf, wf
   real, dimension(:, :, :), intent(out) :: uc, vc, wc
   integer, intent(in) :: ix1, jy1, kz1
   integer, intent(in) :: ix2, jy2, kz2
end subroutine ins_velToCenter
