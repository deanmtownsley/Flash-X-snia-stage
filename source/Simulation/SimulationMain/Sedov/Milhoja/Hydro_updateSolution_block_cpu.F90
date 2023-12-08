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

#include "UHD.h"

!> This code was ported from the Milhoja/Sedov C++ code rather than the Flash-X
!! Fortran code so that performance results from both software can be compared
!! cleanly.  Therefore, if any changes are made here, we should determine if
!! similar changes need to be made in the C++ code.
!!
!! @todo Compare this code against C++ before launching a study.
subroutine Hydro_updateSolution_block_cpu(lo, hi, flX, flY, flZ, loFl, U, loU)
    implicit none

    integer, intent(IN)    :: lo(1:MDIM)
    integer, intent(IN)    :: hi(1:MDIM)
    integer, intent(IN)    :: loFl(1:MDIM+1)
    real,    intent(IN)    :: flX(loFl(1):, loFl(2):, loFl(3):, loFl(4):)
    real,    intent(IN)    :: flY(loFl(1):, loFl(2):, loFl(3):, loFl(4):)
    real,    intent(IN)    :: flZ(loFl(1):, loFl(2):, loFl(3):, loFl(4):)
    integer, intent(IN)    :: loU(1:MDIM+1)
    real,    intent(INOUT) :: U(loU(1):, loU(2):, loU(3):, loU(4):)

#ifdef EINT_VAR
    real :: norm2_sqr
#endif
    real :: densOld
    real :: densNew
    real :: densNew_inv

    integer :: i, j, k

    do         k = lo(KAXIS), hi(KAXIS)
        do     j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)
                ! Update density first
                densOld = U(i, j, k, DENS_VAR)
#if NDIM == 1
                densNew =   densOld                          &
                          + flX(i,   j, k, HY_DENS_FLUX)     &
                          - flX(i+1, j, k, HY_DENS_FLUX)
#elif NDIM == 2
                densNew =   densOld                          &
                          + flX(i,   j,   k, HY_DENS_FLUX)   &
                          - flX(i+1, j,   k, HY_DENS_FLUX)   &
                          + flY(i,   j,   k, HY_DENS_FLUX)   &
                          - flY(i,   j+1, k, HY_DENS_FLUX)
#elif NDIM == 3
                densNew =   densOld                          &
                          + flX(i,   j,   k,   HY_DENS_FLUX) &
                          - flX(i+1, j,   k,   HY_DENS_FLUX) &
                          + flY(i,   j,   k,   HY_DENS_FLUX) &
                          - flY(i,   j+1, k,   HY_DENS_FLUX) &
                          + flZ(i,   j,   k,   HY_DENS_FLUX) &
                          - flZ(i,   j,   k+1, HY_DENS_FLUX)
#endif
                U(i, j, k, DENS_VAR) = densNew
                densNew_inv = 1.0 / densNew

                ! velocities and total energy can be updated independently
                ! using density result
#if NDIM == 1
                U(i, j, k, VELX_VAR) = (    U(i,   j, k, VELX_VAR) * densOld     &
                                        + flX(i,   j, k, HY_XMOM_FLUX)           &
                                        - flX(i+1, j, k, HY_XMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELY_VAR) = (    U(i,   j, k, VELY_VAR) * densOld     &
                                        + flX(i,   j, k, HY_YMOM_FLUX)           &
                                        - flX(i+1, j, k, HY_YMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELZ_VAR) = (    U(i,   j, k, VELZ_VAR) * densOld     &
                                        + flX(i,   j, k, HY_ZMOM_FLUX)           &
                                        - flX(i+1, j, k, HY_ZMOM_FLUX) ) * densNew_inv

                U(i, j, k, ENER_VAR) = (    U(i,   j, k, ENER_VAR) * densOld     &
                                        + flX(i,   j, k, HY_ENER_FLUX)           &
                                        - flX(i+1, j, k, HY_ENER_FLUX) ) * densNew_inv
#elif NDIM == 2
                U(i, j, k, VELX_VAR) = (    U(i,   j,   k, VELX_VAR) * densOld   &
                                        + flX(i,   j,   k, HY_XMOM_FLUX)         &
                                        - flX(i+1, j,   k, HY_XMOM_FLUX)         &
                                        + flY(i,   j,   k, HY_XMOM_FLUX)         &
                                        - flY(i,   j+1, k, HY_XMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELY_VAR) = (    U(i,   j,   k, VELY_VAR) * densOld   &
                                        + flX(i,   j,   k, HY_YMOM_FLUX)         &
                                        - flX(i+1, j,   k, HY_YMOM_FLUX)         &
                                        + flY(i,   j,   k, HY_YMOM_FLUX)         &
                                        - flY(i,   j+1, k, HY_YMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELZ_VAR) = (    U(i,   j,   k, VELZ_VAR) * densOld   &
                                        + flX(i,   j,   k, HY_ZMOM_FLUX)         &
                                        - flX(i+1, j,   k, HY_ZMOM_FLUX)         &
                                        + flY(i,   j,   k, HY_ZMOM_FLUX)         &
                                        - flY(i,   j+1, k, HY_ZMOM_FLUX) ) * densNew_inv

                U(i, j, k, ENER_VAR) = (    U(i,   j,   k, ENER_VAR) * densOld   &
                                        + flX(i,   j,   k, HY_ENER_FLUX)         &
                                        - flX(i+1, j,   k, HY_ENER_FLUX)         &
                                        + flY(i,   j,   k, HY_ENER_FLUX)         &
                                        - flY(i,   j+1, k, HY_ENER_FLUX) ) * densNew_inv
#elif NDIM == 3
                U(i, j, k, VELX_VAR) = (    U(i,   j,   k,   VELX_VAR) * densOld &
                                        + flX(i,   j,   k,   HY_XMOM_FLUX)       &
                                        - flX(i+1, j,   k,   HY_XMOM_FLUX)       &
                                        + flY(i,   j,   k,   HY_XMOM_FLUX)       &
                                        - flY(i,   j+1, k,   HY_XMOM_FLUX)       &
                                        + flZ(i,   j,   k,   HY_XMOM_FLUX)       &
                                        - flZ(i,   j,   k+1, HY_XMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELY_VAR) = (    U(i,   j,   k,   VELY_VAR) * densOld &
                                        + flX(i,   j,   k,   HY_YMOM_FLUX)       &
                                        - flX(i+1, j,   k,   HY_YMOM_FLUX)       &
                                        + flY(i,   j,   k,   HY_YMOM_FLUX)       &
                                        - flY(i,   j+1, k,   HY_YMOM_FLUX)       &
                                        + flZ(i,   j,   k,   HY_YMOM_FLUX)       &
                                        - flZ(i,   j,   k+1, HY_YMOM_FLUX) ) * densNew_inv

                U(i, j, k, VELZ_VAR) = (    U(i,   j,   k,   VELZ_VAR) * densOld &
                                        + flX(i,   j,   k,   HY_ZMOM_FLUX)       &
                                        - flX(i+1, j,   k,   HY_ZMOM_FLUX)       &
                                        + flY(i,   j,   k,   HY_ZMOM_FLUX)       &
                                        - flY(i,   j+1, k,   HY_ZMOM_FLUX)       &
                                        + flZ(i,   j,   k,   HY_ZMOM_FLUX)       &
                                        - flZ(i,   j,   k+1, HY_ZMOM_FLUX) ) * densNew_inv

                U(i, j, k, ENER_VAR) = (    U(i,   j,   k,   ENER_VAR) * densOld &
                                        + flX(i,   j,   k,   HY_ENER_FLUX)       &
                                        - flX(i+1, j,   k,   HY_ENER_FLUX)       &
                                        + flY(i,   j,   k,   HY_ENER_FLUX)       &
                                        - flY(i,   j+1, k,   HY_ENER_FLUX)       &
                                        + flZ(i,   j,   k,   HY_ENER_FLUX)       &
                                        - flZ(i,   j,   k+1, HY_ENER_FLUX) ) * densNew_inv
#endif

#ifdef EINT_VAR
                ! Compute energy correction from new velocities and energy
                norm2_sqr =   U(i, j, k, VELX_VAR) * U(i, j, k, VELX_VAR) &
                            + U(i, j, k, VELY_VAR) * U(i, j, k, VELY_VAR) &
                            + U(i, j, k, VELZ_VAR) * U(i, j, k, VELZ_VAR)
                U(i, j, k, EINT_VAR) = U(i, j, k, ENER_VAR) - (0.5 * norm2_sqr);
#endif
            end do
        end do
    end do

end subroutine Hydro_updateSolution_block_cpu

