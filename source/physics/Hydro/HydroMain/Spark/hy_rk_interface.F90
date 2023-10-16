!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup HydroSpark
!!
!! @brief Internal interface for Spark Hydro subroutines
!<
module hy_rk_interface

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"

  interface
     subroutine hy_rk_initSolnScratch(Uin, starState, tmpState, blkLimitsGC, stage)
       real, dimension(:,:,:,:), pointer :: Uin, starState, tmpState
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: blkLimitsGC
       integer, intent(IN) :: stage
     end subroutine hy_rk_initSolnScratch
  end interface

  interface
     subroutine hy_rk_permutateLimits(blkLimits, blkLimitsGC, permLimits, permLimitsGC)
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: blkLimits, blkLimitsGC
       integer, dimension(LOW:HIGH, MDIM, NDIM, MAXSTAGE), intent(OUT) :: permLimits, permLimitsGC
     end subroutine hy_rk_permutateLimits
  end interface

  interface
     subroutine hy_rk_calcLimits(stage, blkLimits, limits)
       integer, intent(IN) :: stage
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: blkLimits
       integer, dimension(LOW:HIGH, MDIM), intent(OUT) :: limits
     end subroutine hy_rk_calcLimits
  end interface

  interface
     subroutine hy_rk_getFaceFlux (starState, flat3d, flx, fly, flz, &
                                   lim, limgc, stage, &
                                   hybridRiemann, cvisc, C_hyp, tinyZero, smalld, smallp, smallx, &
                                   a_flux, a_flat, a_shck, a_rope, a_uPlus, a_uMinus)
       real, dimension(:,:,:,:), pointer :: starState, flx, fly, flz
       real, dimension(:,:,:), pointer :: flat3d
       integer, dimension(LOW:HIGH, MDIM, NDIM, MAXSTAGE), intent(IN) :: lim, limgc
       integer, intent(IN) :: stage
       logical, intent(IN) :: hybridRiemann
       real, intent(IN) :: cvisc, C_hyp, smalld, smallp, smallx, tinyZero
       real, dimension(:), target, intent(IN) :: a_flux, a_flat, a_shck, a_rope, a_uPlus, a_uMinus
     end subroutine  hy_rk_getFaceFlux
  end interface

  interface
     subroutine hy_rk_updateSoln (starState, tmpState, rk_coeffs, &
                                  grav, flx, fly, flz, &
                                  deltas, fareaX, fareaY, fareaZ, cvol, xCenter, &
                                  xLeft, xRight, yLeft, yRight, &
                                  geometry, &
                                  smalle, smalld, alphaGLM, C_hyp, &
                                  dt, dtOld, limits)
       real, dimension(:,:,:,:), pointer :: starState, tmpState, flx, fly, flz
       real, dimension(:,:,:,:), pointer :: grav
       real, dimension(:,:,:), pointer :: fareaX, fareaY, fareaZ, cvol
       real, dimension(:), pointer :: xCenter, xLeft, xRight, yLeft, yRight
       real, dimension(3), intent(IN) :: rk_coeffs
       real, dimension(MDIM), intent(IN)  :: deltas
       integer, intent(IN) :: geometry
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
       real, intent(IN) :: smalle, smalld, alphaGLM, C_hyp, dt, dtOld
     end subroutine hy_rk_updateSoln
  end interface

  interface
     subroutine hy_rk_getGraveAccel(starState, grav, radCenter, thtCenter, deltas, geometry, blkLimitsGC)
       real, dimension(:,:,:,:), pointer :: starState, grav
       real, dimension(:), pointer :: radCenter, thtCenter
       real, dimension(MDIM), intent(IN)  :: deltas
       integer, intent(IN) :: geometry
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: blkLimitsGC
     end subroutine hy_rk_getGraveAccel
  end interface

  interface
     subroutine hy_rk_correctFluxes(Uin, blkLimits, &
                                    fluxBufX, fluxBufY, fluxBufZ, &
                                    deltas, fareaX, fareaY, fareaZ, cvol, xCenter, &
                                    xLeft, xRight, yLeft, yRight, &
                                    geometry, &
                                    smalle, smalld, &
                                    dt, isFlux)
       real, pointer :: Uin(:,:,:,:)
       real, pointer :: fluxBufX(:,:,:,:), fluxBufY(:,:,:,:), fluxBufZ(:,:,:,:)
       integer, dimension(LOW:HIGH, MDIM), intent(IN) :: blkLimits
       real, dimension(:,:,:), pointer :: fareaX, fareaY, fareaZ, cvol
       real, dimension(:), pointer :: xCenter, xLeft, xRight, yLeft, yRight
       real, dimension(MDIM), intent(IN)  :: deltas
       integer, intent(IN) :: geometry
       real, intent(IN) :: smalle, smalld, dt
       logical, intent(IN) :: isFlux
     end subroutine hy_rk_correctFluxes
  end interface


  interface
     subroutine hy_rk_shockDetect(Uin, blkLimitsGC, tinyZero)
       real, pointer, dimension(:, :, :, :) :: Uin
       integer, intent(IN) :: blkLimitsGC(LOW:HIGH, MDIM)
       real, intent(IN) :: tinyZero
     end subroutine hy_rk_shockDetect
  end interface

  interface
     subroutine hy_rk_getFlatteningLimiter(is_flattening, starState, flat3d, limits)
        logical, intent(IN) :: is_flattening
        real, dimension(:,:,:,:), pointer :: starState
        real, dimension(:,:,:), pointer :: flat3d
        integer, intent(IN), dimension(LOW:HIGH, MDIM) :: limits
     end subroutine hy_rk_getFlatteningLimiter
   end interface

end module hy_rk_interface
