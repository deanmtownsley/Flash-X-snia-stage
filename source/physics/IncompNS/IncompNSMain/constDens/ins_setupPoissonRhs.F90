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
!!**
#include "Simulation.h"
#include "constants.h"

subroutine ins_setupPoissonRhs_constdens(divu, dt)
   implicit none
   real, dimension(:, :, :), intent(inout) :: divu
   real, intent(in) :: dt
   divu = divu/dt
end subroutine ins_setupPoissonRhs_constdens

subroutine ins_setupPoissonRhs_vardens(divu, &
                                       sigx, sigy, sigz, &
                                       pxn1, pyn1, pzn1, &
                                       pxn2, pyn2, pzn2, &
                                       rhox, rhoy, rhoz, &
                                       rhoGas, dt, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)

   implicit none
   real, dimension(:, :, :), intent(inout) :: divu
   real, dimension(:, :, :), intent(in) :: sigx, sigy, sigz
   real, dimension(:, :, :), intent(in) :: rhox, rhoy, rhoz, pxn1, pyn1, pzn1, &
                                           pxn2, pyn2, pzn2
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: dt, dx, dy, dz, rhoGas
end subroutine ins_setupPoissonRhs_vardens

