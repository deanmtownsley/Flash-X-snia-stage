!!****if* source/Grid/GridMain/Grid_addToVar
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
!!  Grid_addToVar
!!
!! SYNOPSIS
!!
!!  call Grid_addToVar(integer(in) :: srcVar,
!!                     integer(in) :: destVar,
!!                     real(in) :: multFactor,
!!                     logical(in) :: reset)
!!
!! DESCRIPTION
!!   Compute solnData(srcVar,:,:,:)*multFactor and save in solnData(destVar,:,:,:).
!!
!!   If reset is true, the destination variable is first zeroed;
!!   otherwise the product is added to the existing values of destVar.
!!
!!   The operation is applied to interior cells of all LEAF blocks.
!!
!! ARGUMENTS
!!
!!
!!   srcVar : the state variables to be used in the RHS of the expression
!!
!!   destVar : the state variables to be used in the LHS of the expression
!!
!!   multFactor : multiplication factor
!!
!!   reset : indicates whether the destination variable should be zeroed first
!!
!! NOTES
!!
!!   srcVar == destVar is allowed and behaves as expected iff reset is .FALSE.
!!
!!   For a copy call Grid_addToVar(srcVar, destVar, 1.0, .true.)
!!
!!***

!!REORDER(4): solnData

#include "Simulation.h"
#include "constants.h"  
#include "Particles.h"

subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)
  use Grid_interface,   ONLY : Grid_getTileIterator, &
                               Grid_releaseTileIterator
  use Grid_tile,        ONLY : Grid_tile_t
  use Grid_iterator,    ONLY : Grid_iterator_t
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer, intent(in) :: srcVar, destVar
  real,    intent(in) :: multFactor
  logical, intent(in) :: reset

  real, dimension(:,:,:,:), pointer :: solnData

  type(Grid_tile_t)     :: tileDesc
  type(Grid_iterator_t) :: itor

  integer :: i, j, k
 
  nullify(solnData)

  call Driver_abort("[Grid_addToVar] This update has not been tested")

  call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%getDataPtr(solnData, CENTER)

     if (reset) then
        ! DEV: TODO Why no just set the destVar by overwriting rather than
        !       adding?  We can save an entire loop over interiors with that.
        do       k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
           do    j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
              do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
                 solnData(destVar,i,j,k) = 0.0
              end do
           end do
        end do
     end if
     do       k = tileDesc%limits(LOW, KAXIS), tileDesc%limits(HIGH, KAXIS)
        do    j = tileDesc%limits(LOW, JAXIS), tileDesc%limits(HIGH, JAXIS)
           do i = tileDesc%limits(LOW, IAXIS), tileDesc%limits(HIGH, IAXIS)
              solnData(destVar,i,j,k) =              solnData(destVar, i, j, k) & 
                                        + multFactor*solnData(srcVar,  i, j, k)
           end do
        end do
     end do

     call tileDesc%releaseDataPtr(solnData, CENTER)
  end do
  call Grid_releaseTileIterator(itor)
end subroutine Grid_addToVar

