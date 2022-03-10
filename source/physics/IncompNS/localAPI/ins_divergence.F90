!!***
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
SUBROUTINE ins_divergence(uni, vni, wni, ix1, ix2, jy1, jy2, kz1, kz2, &
                          dx, dy, dz, divv)
   implicit none
   INTEGER, INTENT(IN) :: ix1, ix2, jy1, jy2, kz1, kz2
   REAL, INTENT(IN) :: dx, dy, dz
   REAL, DIMENSION(:, :, :), INTENT(IN) :: uni, vni, wni
   REAL, DIMENSION(:, :, :), INTENT(OUT) :: divv
END SUBROUTINE ins_divergence
