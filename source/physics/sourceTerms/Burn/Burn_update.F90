!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn
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
!!   call Burn_update (Uin, lo, hi, loGC, hiGC,
!!              real, intent(IN) ::  dt  )
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   dt  --       passed to the internal bn_burner module
!!
!!
!!***

!!REORDER(4): Uin

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

subroutine Burn_update (Uin, loGC, blkLimits, dt)

  implicit none


  !args
  integer, dimension(MDIM),intent(in) ::  loGC
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits
  real, intent(in) :: dt
  real,dimension(1:,loGC(IAXIS):, loGC(JAXIS):, loGC(KAXIS):),intent(inout) :: Uin

  return
end subroutine Burn_Update
