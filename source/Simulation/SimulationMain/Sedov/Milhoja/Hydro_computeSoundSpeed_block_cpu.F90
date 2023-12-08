!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

#include "constants.h"
#include "Simulation.h"

!> This code was ported from the Milhoja/Sedov C++ code rather than the Flash-X
!! Fortran code so that performance results from both software can be compared
!! cleanly.  Therefore, if any changes are made here, we should determine if
!! similar changes need to be made in the C++ code.
!!
!! @todo Compare this code against C++ before launching a study.
subroutine Hydro_computeSoundSpeed_block_cpu(lo, hi, U, loU, auxC, loAux)
    implicit none

    integer, intent(IN)  :: lo(1:MDIM)
    integer, intent(IN)  :: hi(1:MDIM)
    integer, intent(IN)  :: loU(1:MDIM+1)
    real,    intent(IN)  :: U(loU(1):, loU(2):, loU(3):, loU(4):)
    integer, intent(IN)  :: loAux(1:MDIM)
    real,    intent(OUT) :: auxC(loAux(1):, loAux(2):, loAux(3):)

    integer :: i, j, k

    do         k = lo(KAXIS)-K3D, hi(KAXIS)+K3D
        do     j = lo(JAXIS)-K2D, hi(JAXIS)+K2D
            do i = lo(IAXIS)-K1D, hi(IAXIS)+K1D
                auxC(i, j, k) = SQRT(  U(i, j, k, GAMC_VAR)   &
                                     * U(i, j, k, PRES_VAR)   &
                                     / U(i, j, k, DENS_VAR) )
            end do
        end do
    end do

end subroutine Hydro_computeSoundSpeed_block_cpu

