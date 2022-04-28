!!****if* source/physics/RadTrans/RadTrans_computeDt
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
!!  NAME 
!!
!!   RadTrans
!!
!!  SYNOPSIS
!!
!!   RadTrans_computeDt(Grid_tile_t(IN) :: tileDesc,
!!                      real,pointer    :: solnData(:,:,:,:),   
!!                      real,(INOUT)    :: dt_radtrans, 
!!                      integer(INOUT)  :: dtMinLoc(:))
!!
!!  DESCRIPTION 
!!    Compute radiative transfer time step
!!
!!  ARGUMENTS
!!    tileDesc      --  meta-information about the tile/block
!!    solnData      --  the physical, solution data from grid
!!    dt_radtrans   --  variable to hold timestep constraint
!!    dt_minloc(5)  --  array to hold limiting zone info:  zone indices
!!
!!***

!!REORDER(4): solnData

subroutine RadTrans_computeDt(tileDesc, solnData, dt_radtrans, dtMinLoc)
  
  use Grid_tile, ONLY : Grid_tile_t

#include "constants.h"

  implicit none

  type(Grid_tile_t), intent(IN) :: tileDesc
  real, pointer :: solnData(:,:,:,:) 
  real, intent(INOUT) :: dt_radtrans
  integer, intent(INOUT)  :: dtMinLoc(5)

  ! Stub implementation
  return

end subroutine RadTrans_computeDt
