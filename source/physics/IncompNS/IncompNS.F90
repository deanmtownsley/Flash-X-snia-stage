!!****f* source/physics/IncompNS/IncompNS
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!! NAME
!!
!!  IncompNS
!!
!!
!! SYNOPSIS
!!
!!  IncompNS(integer(IN) :: blockCount, 
!!      integer(IN) :: blockList(blockCount)
!!      real(IN)    :: timeEndAdv,
!!      real(IN)    :: dt,
!!      real(IN)    :: dtOld,
!!      integer(IN) :: sweepOrder)
!!
!!
!! DESCRIPTION
!! 
!!  Performs INS timestep advancement.
!!
!!  The blockList and blockCount arguments tell this routine on 
!!  which blocks and on how many to operate.  blockList is an 
!!  integer array of size blockCount that contains the local 
!!  block numbers of blocks on which to advance.
!!
!!  dt gives the timestep through which this update should advance,
!!  and timeEndAdv tells the time that this update will reach when
!!  it finishes.  dtOld gives the previously taken timestep.
!!
!! ARGUMENTS
!!
!!  timeEndAdv - dummy consistent with toplayer stub function
!!  dt         - timestep
!!  dtOld      - dummy consistent with toplayer stub function
!!  sweepOrder - dummy argument for the unsplit scheme, just a dummy
!!               variable to be consistent with a top-layer stub function
!!
!!***

subroutine IncompNS(timeEndAdv,dt,dtOld,sweepOrder)

  implicit none

#include "Simulation.h"

  integer, INTENT(IN) :: sweepOrder
  real,    INTENT(IN) :: timeEndAdv, dt, dtOld

end subroutine IncompNS
