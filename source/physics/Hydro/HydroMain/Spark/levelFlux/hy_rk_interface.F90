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
!!REORDER(4):fluxBuf[XYZ]

module hy_rk_interface

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"

  interface
     subroutine  hy_rk_getFaceFlux (hy_starState, blklimits,blkLimitsGC, limits)
       real, pointer,dimension(:,:,:,:) :: hy_starState       
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits,blkLimits,blklimitsGC
       real,dimension(MDIM) :: hy_del
     end subroutine  hy_rk_getFaceFlux
  end interface

  interface
     subroutine hy_rk_updateSoln (hy_starState, hy_tmpState, blkLimits,blklimitsGC,hy_del, dt, dtOld, limits, hy_coeffs)
       real, pointer,dimension(:,:,:,:) :: hy_starState, hy_tmpState
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits,blkLimits,blkLimitsGC
       real,dimension(MDIM) :: hy_del
       real, intent(IN) :: dt, dtOld
       real, dimension(3), intent(IN) :: hy_coeffs
     end subroutine hy_rk_updateSoln
  end interface

  interface
     subroutine hy_rk_getGraveAccel (hy_starState, hy_del,limits,blkLimitsGC)
       real, pointer,dimension(:,:,:,:) :: hy_starState
       real,dimension(MDIM),intent(IN)  :: hy_del
       integer,dimension(LOW:HIGH,MDIM), intent(IN) :: limits, blkLimitsGC
     end subroutine hy_rk_getGraveAccel
  end interface
  
  
  interface
     subroutine hy_rk_correctFluxes(Uin,blkLimits,&
                                    fluxBufX, fluxBufY, fluxBufZ, &
                                    blklimitsGC,level,hy_del, dt, isFlux)
       implicit none
       real, pointer :: Uin(:,:,:,:)
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       real, CONTIGUOUS,TARGET, &
            dimension(1:, blkLimits(LOW,1):, blkLimits(LOW,2):, blkLimits(LOW,3):) :: fluxBufX, fluxBufY, fluxBufZ
       integer,intent(IN) :: level
       real,dimension(MDIM) :: hy_del
       real, intent(IN) :: dt
       logical, intent(IN) :: isFlux ! handle as fluxes rather than flux densities?
     end subroutine hy_rk_correctFluxes
  end interface

end module hy_rk_interface
