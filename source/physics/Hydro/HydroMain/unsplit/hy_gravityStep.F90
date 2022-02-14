#include "Simulation.h"
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
#include "constants.h"

subroutine hy_gravityStep(simTime, dt, dtOld)

  use Grid_interface,      ONLY : Grid_getTileIterator, &
                                  Grid_releaseTileIterator
  use Timers_interface,    ONLY : Timers_start, Timers_stop

  use hy_interface,        ONLY : hy_gravityStepBlk
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Grid_tile,           ONLY : Grid_tile_t

  implicit none

  real, intent(IN) ::  simTime, dt, dtOld

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC
  real,pointer,dimension(:,:,:,:) :: Uout, Uin
  real,dimension(MDIM) :: del

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  logical, save :: fakeTimer = .TRUE.

  nullify(Uin)
  nullify(Uout)

#ifdef DEBUG_DRIVER
  print*,' ***************   HYDRO LEVEL  **********************'
#endif

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  call Timers_start("loop5")
  if (.not. itor%isValid() .AND. fakeTimer) then
     call Timers_start("loop5 body")
     call Timers_stop ("loop5 body")
  end if
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     blkLimits(:,:)   = tileDesc%limits
     blkLimitsGC(:,:) = tileDesc%blkLimitsGC

     call tileDesc%getDataPtr(Uout, CENTER)
     call tileDesc%deltas(del)
     Uin => Uout
     call hy_gravityStepBlk(tileDesc,blkLimitsGC,Uin, blkLimits, Uout, del,simTime, dt, dtOld)
     call tileDesc%releaseDataPtr(Uout, CENTER)
     nullify(Uout)
     fakeTimer = .FALSE.

     call itor%next()
  end do
  call Timers_stop("loop5")
  call Grid_releaseTileIterator(itor)
end subroutine hy_gravityStep
