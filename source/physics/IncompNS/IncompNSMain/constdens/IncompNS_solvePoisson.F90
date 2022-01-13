!!****if* source/physics/IncompNS/IncompNSMain/constdens/IncompNS_solvePoisson
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
!!
!!***
!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"
#include "IncompNS.h"

subroutine IncompNS_solvePoisson(dt)

  use Grid_interface,   ONLY : Grid_fillGuardCells,Grid_getTileIterator,&
                               Grid_releaseTileIterator,Grid_solvePoisson
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use ins_interface,    ONLY : ins_setupPoissonRhs_constdens
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Driver_interface, ONLY : Driver_getNStep
  use IncompNS_data

  implicit none
  include "Flash_mpi.h"
  !---Argument List-------
  real,    INTENT(IN) :: dt

!------------------------------------------------------------------------------------------
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  real, pointer, dimension(:,:,:,:) :: solnData
  integer TA(2),count_rate
  real*8  ET
  real del(MDIM)
  integer :: NStep
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor

!------------------------------------------------------------------------------------------
  nullify(solnData)

  CALL SYSTEM_CLOCK(TA(1),count_rate)

  !---POISSON RHS:-------------------------------------------------------------------------------------
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC
     call tileDesc%deltas(del)
     call tileDesc%getDataPtr(solnData,  CENTER)

     ! Poisson RHS source vector
     call ins_setupPoissonRhs_constdens(solnData(DUST_VAR,:,:,:),dt)

     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData,  CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  !------------------------------------------------------------------------------------------------------
  ! SOLUTION OF POISSON EQUATION FOR PRESSURE:
  ! -------- -- ------- -------- --- --------
  !------------------------------------------------------------------------------------------------------
  call Timers_start("Grid_solvePoisson")
  call Grid_solvePoisson (iSoln=PRES_VAR, iSrc=DUST_VAR, &
                          bcTypes=ins_pressureBC_types, &
                          bcValues=ins_pressureBC_values, &
                          poisfact=ins_poisfact)
  call Timers_stop("Grid_solvePoisson")

  gcMask = .FALSE.
  gcMask(PRES_VAR) = .TRUE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,  &
      maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask, &
      selectBlockType=ACTIVE_BLKS)

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ins_meshMe .eq. MASTER_PE)  write(*,*) 'Total IncompNS Poisson Solve Time =',ET

  return
end subroutine IncompNS_solvePoisson
