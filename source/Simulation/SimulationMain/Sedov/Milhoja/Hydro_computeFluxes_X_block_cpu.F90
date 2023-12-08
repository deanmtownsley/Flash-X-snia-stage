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
!! @todo Klaus has previously pointed out that in certain cases we define
!! variables that are never used.  Get rid of these and see if this can explain
!! the fact that walltimes grow as the Sedov simulation evolves when we run on
!! only the CPU.  Same changes are likely needed for both Y and Z.
!! @todo Should these replace the current, official simpleUnsplit implementations?
!! If not, do the same corrections need to be made to the official code?
subroutine Hydro_computeFluxes_X_block_cpu(dt, lo, hi, deltas, U, loU, auxC, loAux, flX, loFl)
    implicit none

    real,    intent(IN)  :: dt
    integer, intent(IN)  :: lo(1:MDIM)
    integer, intent(IN)  :: hi(1:MDIM)
    real,    intent(IN)  :: deltas(1:MDIM)
    integer, intent(IN)  :: loU(1:MDIM+1)
    real,    intent(IN)  :: U(loU(1):, loU(2):, loU(3):, loU(4):)
    integer, intent(IN)  :: loAux(1:MDIM)
    real,    intent(IN)  :: auxC(loAux(1):, loAux(2):, loAux(3):)
    integer, intent(IN)  :: loFl(1:MDIM+1)
    real,    intent(OUT) :: flX(loFl(1):, loFl(2):, loFl(3):, loFl(4):)

    real    :: sL
    real    :: sR
    real    :: sRsL
    real    :: vn
    real    :: vL
    real    :: vR
    integer :: is
    integer :: iL
    integer :: iR
    real    :: dtdx

    integer :: i, j, k

    dtdx = dt / deltas(IAXIS)

    do         k = lo(KAXIS), hi(KAXIS)
        do     j = lo(JAXIS), hi(JAXIS)
            do i = lo(IAXIS), hi(IAXIS)+K1D
                sL = MIN(U(i-1, j, k, VELX_VAR) - auxC(i-1, j, k), &
                         U(i,   j, k, VELX_VAR) - auxC(i,   j, k))
                sR = MAX(U(i-1, j, k, VELX_VAR) + auxC(i-1, j, k), &
                         U(i,   j, k, VELX_VAR) + auxC(i,   j, k))
                sRsL = sR - sL
                if (sL > 0.0) then
                    vn = U(i-1, j, k, VELX_VAR)
                    is = i - 1
                    iL = i - 1
                    iR = i - 1
                else if (sR < 0.0) then
                    vn = U(i, j, k, VELX_VAR)
                    is = i
                    iL = i
                    iR = i
                else
                    vn = 0.5 * (  U(i-1, j, k, VELX_VAR)  &
                                + U(i,   j, k, VELX_VAR))
                    is = i
                    iL = i-1
                    iR = i
                    if (vn > 0.0) then
                        is = is - 1
                    end if 
                end if

                vL = U(iL, j, k, VELX_VAR)
                vR = U(iR, j, k, VELX_VAR)
                if (iL == iR) then
                    flX(i, j, k, HY_DENS_FLUX) =   vn * U(is, j, k, DENS_VAR)
                    flX(i, j, k, HY_XMOM_FLUX) =   vn * U(is, j, k, DENS_VAR) &
                                                      * U(is, j, k, VELX_VAR) &
                                                 +      U(is, j, k, PRES_VAR)
                    flX(i, j, k, HY_YMOM_FLUX) =   vn * U(is, j, k, DENS_VAR) &
                                                      * U(is, j, k, VELY_VAR)
                    flX(i, j, k, HY_ZMOM_FLUX) =   vn * U(is, j, k, DENS_VAR) &
                                                      * U(is, j, k, VELZ_VAR)
                    flX(i, j, k, HY_ENER_FLUX) =   vn * U(is, j, k, DENS_VAR) &
                                                      * U(is, j, k, ENER_VAR) &
                                                 + vn * U(is, j, k, PRES_VAR)
                else
                    flX(i, j, k, HY_DENS_FLUX) = (  sR * vL * U(iL, j, k, DENS_VAR) &
                                                  - sL * vR * U(iR, j, k, DENS_VAR) &
                                                  + sR*sL*(   U(iR, j, k, DENS_VAR) &
                                                            - U(iL, j, k, DENS_VAR)) ) / sRsL
                    flX(i, j, k, HY_XMOM_FLUX) = (  sR * vL * U(iL, j, k, DENS_VAR) * U(iL, j, k, VELX_VAR)  &
                                                  - sL * vR * U(iR, j, k, DENS_VAR) * U(iR, j, k, VELX_VAR)  &
                                                  + sR*sL*(   U(iR, j, k, DENS_VAR) * U(iR, j, k, VELX_VAR)  &
                                                            - U(iL, j, k, DENS_VAR) * U(iL, j, k, VELX_VAR)) &
                                                  + sR *      U(iL, j, k, PRES_VAR)                          &
                                                  - sL *      U(iR, j, k, PRES_VAR) ) / sRsL
                    flX(i, j, k, HY_YMOM_FLUX) = (  sR * vL * U(iL, j, k, DENS_VAR) * U(iL, j, k, VELY_VAR) &
                                                  - sL * vR * U(iR, j, k, DENS_VAR) * U(iR, j, k, VELY_VAR) &
                                                  + sR*sL*(   U(iR, j, k, DENS_VAR) * U(iR, j, k, VELY_VAR) &
                                                            - U(iL, j, k, DENS_VAR) * U(iL, j, k, VELY_VAR)) ) / sRsL
                    flX(i, j, k, HY_ZMOM_FLUX) = (  sR * vL * U(iL, j, k, DENS_VAR) * U(iL, j, k, VELZ_VAR) &
                                                  - sL * vR * U(iR, j, k, DENS_VAR) * U(iR, j, k, VELZ_VAR) &
                                                  + sR*sL*(   U(iR, j, k, DENS_VAR) * U(iR, j, k, VELZ_VAR) &
                                                            - U(iL, j, k, DENS_VAR) * U(iL, j, k, VELZ_VAR)) ) / sRsL
                    flX(i, j, k, HY_ENER_FLUX) = (  sR * vL * U(iL, j, k, DENS_VAR) * U(iL, j, k, ENER_VAR)  &
                                                  - sL * vR * U(iR, j, k, DENS_VAR) * U(iR, j, k, ENER_VAR)  &
                                                  + sR*sL*(   U(iR, j, k, DENS_VAR) * U(iR, j, k, ENER_VAR)  &
                                                            - U(iL, j, k, DENS_VAR) * U(iL, j, k, ENER_VAR)) &
                                                  + sR * vL * U(iL, j, k, PRES_VAR)                          &
                                                  - sL * vR * U(iR, j, k, PRES_VAR) ) / sRsL
                end if

                flX(i, j, k, HY_DENS_FLUX) = flX(i, j, k, HY_DENS_FLUX) * dtdx
                flX(i, j, k, HY_XMOM_FLUX) = flX(i, j, k, HY_XMOM_FLUX) * dtdx
                flX(i, j, k, HY_YMOM_FLUX) = flX(i, j, k, HY_YMOM_FLUX) * dtdx
                flX(i, j, k, HY_ZMOM_FLUX) = flX(i, j, k, HY_ZMOM_FLUX) * dtdx
                flX(i, j, k, HY_ENER_FLUX) = flX(i, j, k, HY_ENER_FLUX) * dtdx
            end do
        end do
    end do

end subroutine Hydro_computeFluxes_X_block_cpu

