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
SUBROUTINE ins_diffusion2d_vardens(uni, vni, ru1, ix1, ix2, jy1, jy2, dx, dy, ru, rv, &
                                   visc, rhox, rhoy)
   implicit none
   INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2
   REAL, INTENT(IN):: ru1, dx, dy
   REAL, DIMENSION(:, :, :), INTENT(IN):: uni, vni, visc, rhox, rhoy
   REAL, DIMENSION(:, :, :), INTENT(OUT):: ru, rv
END SUBROUTINE ins_diffusion2d_vardens

SUBROUTINE ins_diffusion3d_vardens(uni, vni, wni, tv, ru1, &
                                   ix1, ix2, jy1, jy2, kz1, kz2, &
                                   dx, dy, dz, ru, rv, rw, visc, &
                                   rhox, rhoy, rhoz)
   implicit none
   INTEGER, INTENT(IN):: ix1, ix2, jy1, jy2, kz1, kz2
   REAL, INTENT(IN):: ru1, dx, dy, dz
   REAL, DIMENSION(:, :, :), INTENT(IN):: uni, vni, wni, tv, visc, rhox, rhoy
   REAL, DIMENSION(:, :, :), INTENT(IN):: rhoz
   REAL, DIMENSION(:, :, :), INTENT(OUT):: ru, rv, rw
END SUBROUTINE ins_diffusion3d_vardens
