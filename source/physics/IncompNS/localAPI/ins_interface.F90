!!****if* source/physics/IncompNS/localAPI/ins_interface
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
!! NAME
!!  ins_inteface
!!
!! SYNOPSIS
!!  use ins_interface
!!
!! DESCRIPTION
!!  This is an interface module that defines explicit interfaces for private procedures
!!  within the Incompressible Navier Stokes unit.
!!
!!***
Module ins_interface

#include "constants.h"
#include "Simulation.h"
#include "IncompNS.h"

   interface
      subroutine ins_computeDtLocal(blockID, &
                                    isize, jsize, ksize, &
                                    dx, dy, dz, &
                                    blkLimits, blkLimitsGC, &
                                    facexData, faceyData, &
                                    facezData, &
                                    dtLocal, lminloc)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(2, MDIM), intent(IN) :: blkLimits, blkLimitsGC
         integer :: isize, jsize, ksize
         real :: dx, dy, dz
         real, intent(INOUT) :: dtLocal
         integer, intent(INOUT) :: lminloc(5)
         real, pointer, dimension(:, :, :, :)  :: facexData, faceyData, facezData
      end subroutine ins_computeDtLocal
   end interface

   interface
      SUBROUTINE ins_divergence(uni, vni, wni, ix1, ix2, jy1, jy2, kz1, kz2, &
                                dx, dy, dz, divv)
         implicit none
         INTEGER, INTENT(IN) :: ix1, ix2, jy1, jy2, kz1, kz2
         REAL, INTENT(IN) :: dx, dy, dz
         REAL, DIMENSION(:, :, :), INTENT(IN) :: uni, vni, wni
         REAL, DIMENSION(:, :, :), INTENT(OUT) :: divv
      END SUBROUTINE ins_divergence
   end interface

   interface ins_setupPoissonRhs
      subroutine ins_setupPoissonRhs_constdens(divu, dt)
         implicit none
         real, dimension(:, :, :), intent(inout) :: divu
         real, intent(in) :: dt
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
   end interface

   interface ins_corrector
      subroutine ins_corrector_constdens(uni, vni, wni, pxn1, pyn1, pzn1, p, ix1, ix2, jy1, jy2, kz1, kz2, &
                                         dt, dx, dy, dz)
         implicit none
         INTEGER, INTENT(IN) :: ix1, ix2, jy1, jy2, kz1, kz2
         REAL, INTENT(IN) :: dt, dx, dy, dz
         REAL, DIMENSION(:, :, :), INTENT(IN) :: p
         REAL, DIMENSION(:, :, :), INTENT(INOUT) :: uni, vni, wni, pxn1, pyn1, pzn1
      end subroutine ins_corrector_constdens

      subroutine ins_corrector_vardens(uni, vni, wni, sigx, sigy, sigz, pxn1, pyn1, pzn1, &
                                       pxn2, pyn2, pzn2, &
                                       rhox, rhoy, rhoz, p, rhoGas, dt, dx, dy, dz, &
                                       ix1, ix2, jy1, jy2, kz1, kz2)

         implicit none
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: dt, dx, dy, dz, rhoGas
         real, dimension(:, :, :), intent(in) :: p
         real, dimension(:, :, :), intent(in) :: rhox, rhoy, rhoz
         real, dimension(:, :, :), intent(in) :: sigx, sigy, sigz
         real, dimension(:, :, :), intent(inout) :: pxn1, pxn2, pyn1, pyn2, pzn1, pzn2
         real, dimension(:, :, :), intent(inout) :: uni, vni, wni
      end subroutine ins_corrector_vardens
   end interface

   interface ins_diffusion_vardens
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
   end interface

   interface
      subroutine ins_indicators(u, v, w, pres, divv, ix1, ix2, jy1, jy2, kz1, kz2, vecminaux, vecmaxaux)
         implicit none
         real, dimension(:, :, :), intent(in) :: u, v, w
         real, dimension(:, :, :), intent(in) :: pres, divv
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, dimension(5), intent(inout) :: vecminaux, vecmaxaux
      end subroutine ins_indicators
   end interface

   interface
      subroutine ins_init()
         implicit none
      end subroutine ins_init
   end interface

end Module ins_interface
