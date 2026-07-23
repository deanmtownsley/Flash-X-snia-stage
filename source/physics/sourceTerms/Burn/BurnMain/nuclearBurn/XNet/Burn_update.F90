!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Burn_update
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
!!
!!  Burn_update
!!
!!
!! SYNOPSIS
!!
!!   call Burn_update(
!!      real, intent(IN OUT)  ::  Uin(:,:,:,:),
!!      integer, intent(IN)   :: loGC(:),
!!      integer, intent(IN)   :: blkLimits(:,:),
!!      real, intent(IN)      :: dt
!!   )
!!
!! DESCRIPTION
!!
!!   Update the solution data using the outputs from the burner.
!!
!! ARGUMENTS
!!
!!   Uin       --
!!   loGC      --
!!   blkLimits --
!!   dt        --
!!
!! NOTES
!!
!!  This subroutine assumes that the nuclear energy generation rate is available
!!             at the Uin(ENUC_VAR,:,:,:) calculated from the burner.
!!
!!***

!!REORDER(4): Uin

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

subroutine Burn_update (Uin, loGC, blkLimits, dt)

  implicit none

  !args
  integer, dimension(MDIM), intent(IN) ::  loGC
  real, dimension(1:,loGC(IAXIS):, loGC(JAXIS):, loGC(KAXIS):), intent(IN OUT) :: Uin
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt

   ! TODO: use this in Burn.F90 as like Approx13

end subroutine Burn_Update

