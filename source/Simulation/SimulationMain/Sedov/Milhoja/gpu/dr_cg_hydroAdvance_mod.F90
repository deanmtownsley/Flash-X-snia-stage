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

#ifndef MILHOJA_ENABLE_OPENACC_OFFLOAD
#error "This file should only be compiled if using OpenACC offloading"
#endif

#include "constants.h"
#include "Simulation.h"
#include "UHD.h"

!> A module that encapsulates all compilable variants of the PUD-developed
!! static Fortran routines that are finalized by the CODE GENERATOR for use in
!! hydro advance solution task functions.  They have been "finalized" by the
!! CODE GENERATOR in the sense that macros in the original code have been
!! replaced with the appropriate target-/problem-/platform-specific keys.
!! In addition, the CODE GENERATOR will have effectively removed Flash-X
!! offloading hints and directives in CPU variants; replaced with concrete,
!! compilable offloading directives for GPU variants.
module dr_cg_hydroAdvance_mod
    implicit none
    private

    public :: hy_computeSoundSpeedHll_gpu_oacc
    public :: hy_computeFluxesHll_X_gpu_oacc
    public :: hy_computeFluxesHll_Y_gpu_oacc
    public :: hy_computeFluxesHll_Z_gpu_oacc
    public :: hy_updateSolutionHll_gpu_oacc
    public :: eos_idealGammaDensIe_gpu_oacc 

contains

    !> @todo lo/hi must be local indices.  Make this index-space-agnostic.
    subroutine hy_computeSoundSpeedHll_gpu_oacc(lo, hi, U, auxC)
        !$acc routine vector

        integer, intent(IN)  :: lo(1:MDIM)
        integer, intent(IN)  :: hi(1:MDIM)
        real,    intent(IN)  :: U(:, :, :, :)
        real,    intent(OUT) :: auxC(:, :, :)

        integer :: i, j, k

        !$acc loop vector collapse(3)
        do         k = lo(KAXIS)-K3D, hi(KAXIS)+K3D
            do     j = lo(JAXIS)-K2D, hi(JAXIS)+K2D
                do i = lo(IAXIS)-K1D, hi(IAXIS)+K1D
                    auxC(i, j, k) = SQRT(  U(i, j, k, GAMC_VAR)   &
                                         * U(i, j, k, PRES_VAR)   &
                                         / U(i, j, k, DENS_VAR) )
                end do
            end do
        end do
        !$acc end loop
    end subroutine hy_computeSoundSpeedHll_gpu_oacc

    !> @todo lo/hi must be local indices.  Make this index-space-agnostic.
    subroutine hy_computeFluxesHll_X_gpu_oacc(dt, lo, hi, deltas, U, auxC, flX)
        !$acc routine vector

        real,    intent(IN)  :: dt
        integer, intent(IN)  :: lo(1:MDIM)
        integer, intent(IN)  :: hi(1:MDIM)
        real,    intent(IN)  :: deltas(1:MDIM)
        real,    intent(IN)  :: U(:, :, :, :)
        real,    intent(IN)  :: auxC(:, :, :)
        real,    intent(OUT) :: flX(:, :, :, :)

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
    
        !$acc loop vector collapse(3)
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
                        ! TODO: Klaus mentioned that this can be removed.  Confirm
                        ! and remove.  This extra, unnecessary branch could
                        ! potentially explain a progressive slowing of execution on
                        ! CPU as the simulation evolves.  For example, branch
                        ! prediction on this could be good at the start and
                        ! progressively worse as the shock moves out.
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
        !$acc end loop
    
    end subroutine hy_computeFluxesHll_X_gpu_oacc

    !> @todo lo/hi must be local indices.  Make this index-space-agnostic.
    subroutine hy_computeFluxesHll_Y_gpu_oacc(dt, lo, hi, deltas, U, auxC, flY)
        !$acc routine vector
    
        real,    intent(IN)  :: dt
        integer, intent(IN)  :: lo(1:MDIM)
        integer, intent(IN)  :: hi(1:MDIM)
        real,    intent(IN)  :: deltas(1:MDIM)
        real,    intent(IN)  :: U(:, :, :, :)
        real,    intent(IN)  :: auxC(:, :, :)
        real,    intent(OUT) :: flY(:, :, :, :)
    
        real    :: sL
        real    :: sR
        real    :: sRsL
        real    :: vn
        real    :: vL
        real    :: vR
        integer :: js
        integer :: jL
        integer :: jR
        real    :: dtdy
    
        integer :: i, j, k
    
        dtdy = dt / deltas(JAXIS)
    
        !$acc loop vector collapse(3)
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)+K2D
                do i = lo(IAXIS), hi(IAXIS)
                    sL = MIN(U(i, j-1, k, VELY_VAR) - auxC(i, j-1, k), &
                             U(i, j,   k, VELY_VAR) - auxC(i, j,   k))
                    sR = MAX(U(i, j-1, k, VELY_VAR) + auxC(i, j-1, k), &
                             U(i, j,   k, VELY_VAR) + auxC(i, j,   k))
                    sRsL = sR - sL
                    if (sL > 0.0) then
                        vn = U(i, j-1, k, VELY_VAR)
                        js = j - 1
                        jL = j - 1
                        jR = j - 1
                    else if (sR < 0.0) then
                        vn = U(i, j, k, VELY_VAR)
                        js = j
                        jL = j
                        jR = j
                    else
                        vn = 0.5 * (  U(i, j-1, k, VELY_VAR)  &
                                    + U(i, j,   k, VELY_VAR))
                        js = j
                        jL = j - 1
                        jR = j
                        ! TODO: Klaus mentioned that this can be removed.  Confirm
                        ! and remove.  This extra, unnecessary branch could
                        ! potentially explain a progressive slowing of execution on
                        ! CPU as the simulation evolves.  For example, branch
                        ! prediction on this could be good at the start and
                        ! progressively worse as the shock moves out.
                        if (vn > 0.0) then
                            js = js - 1
                        end if
                    end if
    
                    vL = U(i, jL, k, VELY_VAR)
                    vR = U(i, jR, k, VELY_VAR)
                    if (jL == jR) then
                        flY(i, j, k, HY_DENS_FLUX) =   vn * U(i, js, k, DENS_VAR)
                        flY(i, j, k, HY_XMOM_FLUX) =   vn * U(i, js, k, DENS_VAR) &
                                                          * U(i, js, k, VELX_VAR)
                        flY(i, j, k, HY_YMOM_FLUX) =   vn * U(i, js, k, DENS_VAR) &
                                                          * U(i, js, k, VELY_VAR) &
                                                     +      U(i, js, k, PRES_VAR)
                        flY(i, j, k, HY_ZMOM_FLUX) =   vn * U(i, js, k, DENS_VAR) &
                                                          * U(i, js, k, VELZ_VAR)
                        flY(i, j, k, HY_ENER_FLUX) =   vn * U(i, js, k, DENS_VAR) &
                                                          * U(i, js, k, ENER_VAR) &
                                                     + vn * U(i, js, k, PRES_VAR)
                    else
                        flY(i, j, k, HY_DENS_FLUX) = (  sR * vL * U(i, jL, k, DENS_VAR) &
                                                      - sL * vR * U(i, jR, k, DENS_VAR) &
                                                      + sR*sL*(   U(i, jR, k, DENS_VAR) &
                                                               -  U(i, jL, k, DENS_VAR))) / sRsL
                        flY(i, j, k, HY_XMOM_FLUX) = (  sR * vL * U(i, jL, k, DENS_VAR) * U(i, jL, k, VELX_VAR) &
                                                      - sL * vR * U(i, jR, k, DENS_VAR) * U(i, jR, k, VELX_VAR) &
                                                      + sR*sL*(   U(i, jR, k, DENS_VAR) * U(i, jR, k, VELX_VAR) &
                                                               -  U(i, jL, k, DENS_VAR) * U(i, jL, k, VELX_VAR)) ) / sRsL
                        flY(i, j, k, HY_YMOM_FLUX) = (  sR * vL * U(i, jL, k, DENS_VAR) * U(i, jL, k, VELY_VAR)  &
                                                      - sL * vR * U(i, jR, k, DENS_VAR) * U(i, jR, k, VELY_VAR)  &
                                                      + sR*sL*(   U(i, jR, k, DENS_VAR) * U(i, jR, k, VELY_VAR)  &
                                                               -  U(i, jL, k, DENS_VAR) * U(i, jL, k, VELY_VAR)) &
                                                      + sR *      U(i, jL, k, PRES_VAR)                          &
                                                      - sL *      U(i, jR, k, PRES_VAR) ) / sRsL;
                        flY(i, j, k, HY_ZMOM_FLUX) = (  sR * vL * U(i, jL, k, DENS_VAR) * U(i, jL, k, VELZ_VAR) &
                                                      - sL * vR * U(i, jR, k, DENS_VAR) * U(i, jR, k, VELZ_VAR) &
                                                      + sR*sL*(   U(i, jR, k, DENS_VAR) * U(i, jR, k, VELZ_VAR) &
                                                               -  U(i, jL, k, DENS_VAR) * U(i, jL, k, VELZ_VAR)) ) / sRsL
                        flY(i, j, k, HY_ENER_FLUX) = (  sR * vL * U(i, jL, k, DENS_VAR) * U(i, jL, k, ENER_VAR)  &
                                                      - sL * vR * U(i, jR, k, DENS_VAR) * U(i, jR, k, ENER_VAR)  &
                                                      + sR*sL*(   U(i, jR, k, DENS_VAR) * U(i, jR, k, ENER_VAR)  &
                                                               -  U(i, jL, k, DENS_VAR) * U(i, jL, k, ENER_VAR)) &
                                                      + sR * vL * U(i, jL, k, PRES_VAR)                          &
                                                      - sL * vR * U(i, jR, k, PRES_VAR) ) /sRsL
                    end if
    
                    flY(i, j, k, HY_DENS_FLUX) = flY(i, j, k, HY_DENS_FLUX) * dtdy
                    flY(i, j, k, HY_XMOM_FLUX) = flY(i, j, k, HY_XMOM_FLUX) * dtdy
                    flY(i, j, k, HY_YMOM_FLUX) = flY(i, j, k, HY_YMOM_FLUX) * dtdy
                    flY(i, j, k, HY_ZMOM_FLUX) = flY(i, j, k, HY_ZMOM_FLUX) * dtdy
                    flY(i, j, k, HY_ENER_FLUX) = flY(i, j, k, HY_ENER_FLUX) * dtdy
                end do
            end do
        end do
        !$acc end loop
    
    end subroutine hy_computeFluxesHll_Y_gpu_oacc

    subroutine hy_computeFluxesHll_Z_gpu_oacc(dt, lo, hi, deltas, U, auxC, flZ)
        !$acc routine vector
    
        real,    intent(IN)  :: dt
        integer, intent(IN)  :: lo(1:MDIM)
        integer, intent(IN)  :: hi(1:MDIM)
        real,    intent(IN)  :: deltas(1:MDIM)
        real,    intent(IN)  :: U(:, :, :, :)
        real,    intent(IN)  :: auxC(:, :, :)
        real,    intent(OUT) :: flZ(:, :, :, :)
    
        real    :: sL
        real    :: sR
        real    :: sRsL
        real    :: vn
        real    :: vL
        real    :: vR
        integer :: ks
        integer :: kL
        integer :: kR
        real    :: dtdz
    
        integer :: i, j, k
    
        dtdz = dt / deltas(KAXIS)
    
        !$acc loop vector collapse(3)
        do         k = lo(KAXIS), hi(KAXIS)+K3D
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    sL = MIN(U(i, j, k-1, VELZ_VAR) - auxC(i, j, k-1), &
                             U(i, j, k,   VELZ_VAR) - auxC(i, j, k))
                    sR = MAX(U(i, j, k-1, VELZ_VAR) + auxC(i, j, k-1), &
                             U(i, j, k,   VELZ_VAR) + auxC(i, j, k))
                    sRsL = sR - sL
                    if (sL > 0.0) then
                        vn = U(i, j, k-1, VELZ_VAR)
                        ks = k - 1
                        kL = k - 1
                        kR = k - 1
                    else if (sR < 0.0) then
                        vn = U(i, j, k, VELZ_VAR)
                        ks = k
                        kL = k
                        kR = k
                    else
                        vn = 0.5 * (  U(i, j, k-1, VELZ_VAR)  &
                                    + U(i, j, k,   VELZ_VAR))
                        ks = k
                        kL = k-1
                        kR = k
                        ! TODO: Klaus mentioned that this can be removed.  Confirm
                        ! and remove.  This extra, unnecessary branch could
                        ! potentially explain a progressive slowing of execution on
                        ! CPU as the simulation evolves.  For example, branch
                        ! prediction on this could be good at the start and
                        ! progressively worse as the shock moves out.
                        if (vn > 0.0) then
                          ks = ks - 1
                        end if
                    end if
    
                    vL = U(i, j, kL, VELZ_VAR)
                    vR = U(i, j, kR, VELZ_VAR)
                    if (kL == kR) then
                        flZ(i, j, k, HY_DENS_FLUX) =   vn * U(i, j, ks, DENS_VAR)
                        flZ(i, j, k, HY_XMOM_FLUX) =   vn * U(i, j, ks, DENS_VAR) &
                                                          * U(i, j, ks, VELX_VAR)
                        flZ(i, j, k, HY_YMOM_FLUX) =   vn * U(i, j, ks, DENS_VAR) &
                                                          * U(i, j, ks, VELY_VAR)
                        flZ(i, j, k, HY_ZMOM_FLUX) =   vn * U(i, j, ks, DENS_VAR) &
                                                          * U(i, j, ks, VELZ_VAR) &
                                                     +      U(i, j, ks, PRES_VAR)
                        flZ(i, j, k, HY_ENER_FLUX) =   vn * U(i, j, ks, DENS_VAR) &
                                                          * U(i, j, ks, ENER_VAR) &
                                                     + vn * U(i, j, ks, PRES_VAR)
                    else
                        flZ(i, j, k, HY_DENS_FLUX) = (  sR * vL * U(i, j, kL, DENS_VAR) &
                                                      - sL * vR * U(i, j, kR, DENS_VAR) &
                                                      + sR*sL*(   U(i, j, kR, DENS_VAR) &
                                                               -  U(i, j, kL, DENS_VAR))) / sRsL
                        flZ(i, j, k, HY_XMOM_FLUX) = (  sR * vL * U(i, j, kL, DENS_VAR) * U(i, j, kL, VELX_VAR) &
                                                      - sL * vR * U(i, j, kR, DENS_VAR) * U(i, j, kR, VELX_VAR) &
                                                      + sR*sL*(   U(i, j, kR, DENS_VAR) * U(i, j, kR, VELX_VAR) &
                                                               -  U(i, j, kL, DENS_VAR) * U(i, j, kL, VELX_VAR)) ) / sRsL
                        flZ(i, j, k, HY_YMOM_FLUX) = (  sR * vL * U(i, j, kL, DENS_VAR) * U(i, j, kL, VELY_VAR) &
                                                      - sL * vR * U(i, j, kR, DENS_VAR) * U(i, j, kR, VELY_VAR) &
                                                      + sR*sL*(   U(i, j, kR, DENS_VAR) * U(i, j, kR, VELY_VAR) &
                                                               -  U(i, j, kL, DENS_VAR) * U(i, j, kL, VELY_VAR)) ) / sRsL
                        flZ(i, j, k, HY_ZMOM_FLUX) = (  sR * vL * U(i, j, kL, DENS_VAR) * U(i, j, kL, VELZ_VAR)  &
                                                      - sL * vR * U(i, j, kR, DENS_VAR) * U(i, j, kR, VELZ_VAR)  &
                                                      + sR*sL*(   U(i, j, kR, DENS_VAR) * U(i, j, kR, VELZ_VAR)  &
                                                               -  U(i, j, kL, DENS_VAR) * U(i, j, kL, VELZ_VAR)) &
                                                      + sR *      U(i, j, kL, PRES_VAR)                          &
                                                      - sL *      U(i, j, kR, PRES_VAR) ) / sRsL
                        flZ(i, j, k, HY_ENER_FLUX) = (  sR * vL * U(i, j, kL, DENS_VAR) * U(i, j, kL, ENER_VAR)  &
                                                      - sL * vR * U(i, j, kR, DENS_VAR) * U(i, j, kR, ENER_VAR)  &
                                                      + sR*sL*(   U(i, j, kR, DENS_VAR) * U(i, j, kR, ENER_VAR)  &
                                                               -  U(i, j, kL, DENS_VAR) * U(i, j, kL, ENER_VAR)) &
                                                      + sR * vL * U(i, j, kL, PRES_VAR)                          &
                                                      - sL * vR * U(i, j, kR, PRES_VAR) ) / sRsL
                    end if
    
                    flZ(i, j, k, HY_DENS_FLUX) = flZ(i, j, k, HY_DENS_FLUX) * dtdz
                    flZ(i, j, k, HY_XMOM_FLUX) = flZ(i, j, k, HY_XMOM_FLUX) * dtdz
                    flZ(i, j, k, HY_YMOM_FLUX) = flZ(i, j, k, HY_YMOM_FLUX) * dtdz
                    flZ(i, j, k, HY_ZMOM_FLUX) = flZ(i, j, k, HY_ZMOM_FLUX) * dtdz
                    flZ(i, j, k, HY_ENER_FLUX) = flZ(i, j, k, HY_ENER_FLUX) * dtdz
                end do
            end do
        end do
        !$acc end loop
    
    end subroutine hy_computeFluxesHll_Z_gpu_oacc

    subroutine hy_updateSolutionHll_gpu_oacc(lo, hi, flX, flY, flZ, U)
        !$acc routine vector
    
        integer, intent(IN)    :: lo(1:MDIM)
        integer, intent(IN)    :: hi(1:MDIM)
        real,    intent(IN)    :: flX(:, :, :, :)
        real,    intent(IN)    :: flY(:, :, :, :)
        real,    intent(IN)    :: flZ(:, :, :, :)
        real,    intent(INOUT) :: U(:, :, :, :)

    #ifdef EINT_VAR
        real :: norm2_sqr
    #endif
        real :: densOld
        real :: densNew
        real :: densNew_inv
    
        integer :: i, j, k
    
        !$acc loop vector collapse(3)
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
        !$acc end loop
    
    !> @todo lo/hi must be local indices.  Make this index-space-agnostic.
    end subroutine hy_updateSolutionHll_gpu_oacc

    !> @todo lo/hi must be local indices.  Make this index-space-agnostic.
    subroutine eos_idealGammaDensIe_gpu_oacc(lo, hi, U)
        !$acc routine vector
    
        ! Taken from Flash-X Sedov/setup_params file 
        ! TODO: Get this from Eos_data as eos_gamma?
        real, parameter :: EOS_GAMMA = 1.4
    
        ! Taken from Flash-X Sedov/setup_params file 
        ! TODO: Get this from Eos_data as eos_singleSpeciesA?
        real, parameter :: SINGLE_SPECIES_A = 1.0
    
        ! Taken from Flash-X Physical Constants
        ! TODO: Get this from Eos_data as eos_gasConstant?
        real, parameter :: GAS_CONSTANT = 8.3144598e7  ! J/mol/K
    
        integer, intent(IN)    :: lo(1:MDIM)
        integer, intent(IN)    :: hi(1:MDIM)
        real,    intent(INOUT) :: U(:, :, :, :)
    
        real    :: ggprod_inv
        real    :: gamma_m_1_inv
        integer :: i, j, k
    
        ggprod_inv    = (EOS_GAMMA - 1.0) / GAS_CONSTANT
        gamma_m_1_inv =  EOS_GAMMA - 1.0
    
        !$acc loop vector collapse(3)
        do         k = lo(KAXIS), hi(KAXIS)
            do     j = lo(JAXIS), hi(JAXIS)
                do i = lo(IAXIS), hi(IAXIS)
                    U(i, j, k, PRES_VAR) =   U(i, j, k, DENS_VAR) &
                                           * U(i, j, k, EINT_VAR) &
                                           * gamma_m_1_inv
                    U(i, j, k, TEMP_VAR) =   U(i, j, k, EINT_VAR) &
                                           * ggprod_inv           &
                                           * SINGLE_SPECIES_A
                end do
            end do
        end do
        !$acc end loop
    
    end subroutine eos_idealGammaDensIe_gpu_oacc

end module dr_cg_hydroAdvance_mod

