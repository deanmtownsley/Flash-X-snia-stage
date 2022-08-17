!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_interface
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
!!   hy_rk_interface
!!
!! SYNOPSIS
!!   use hy_rk_interface : ONLY
!!
!!  DESCRIPTION
!!    Interface for internal Spark Hydro subroutines
!!
!!***
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
     subroutine hy_rk_updateSoln (hy_starState, hy_tmpState, blkLimits,blklimitsGC,level,hy_del, dt, dtOld, limits, coeffs)
       real, pointer,dimension(:,:,:,:) :: hy_starState, hy_tmpState
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits,blkLimits,blkLimitsGC
       integer,intent(IN) :: level
       real,dimension(MDIM) :: hy_del
       real, intent(IN) :: dt, dtOld
       real, dimension(3), intent(IN) :: coeffs
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
     subroutine hy_rk_correctFluxes(Uin,blkLimits,blklimitsGC,level,hy_del, dt)
       real, pointer :: Uin(:,:,:,:)
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       integer,intent(IN) :: level
       real,dimension(MDIM) :: hy_del
       real, intent(IN) :: dt
     end subroutine hy_rk_correctFluxes
  end interface

end module hy_rk_interface
