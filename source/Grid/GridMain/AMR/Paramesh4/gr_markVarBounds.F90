!!****if* source/Grid/GridMain/Paramesh4/gr_markVarBounds
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
!!  gr_markVarBounds
!!
!! SYNOPSIS
!!  gr_markVarBounds(  integer(in) :: Var
!!                     real(in)    :: var_bnd_min,
!!                     real(in)    :: var_bnd_max,
!!                     integer(in) :: lref )
!!
!! PURPOSE
!!  Refine all blocks for which a given variable (Var) exceeds or falls
!!  between bounds defined by var_bnd_min and var_bnd_max.
!!
!! ARGUMENTS
!!  Var -    the variable of interest
!!
!!  var_bnd_min  -  lower limit for bound
!!
!!  var_bnd_max  -  upper limit for bound
!!
!!   lref -       If > 0, bring all qualifying blocks to this level of refinement.
!!
!!               If <= 0, refine qualifying blocks once.
!!
!! NOTES
!!
!!  This routine has not yet been tested and should be used only as a guideline for
!!  a user's implementation.
!!
!!***

!!REORDER(5):unk

subroutine gr_markVarBounds(Var, var_bnd_min, var_bnd_max, lref)

   use tree, ONLY: refine, derefine, lrefine, lnblocks, nodetype, lrefine_min, stay
   use physicaldata, ONLY: unk
   implicit none
#include "constants.h"
#include "Simulation.h"
! Arguments

   integer, intent(IN) :: Var
   real, intent(IN) :: var_bnd_min, var_bnd_max
   integer, intent(IN) :: lref

! Local data

   integer :: b
   logical :: Grid_mark

!-------------------------------------------------------------------------------

! Loop over all leaf-node blocks.

   do b = 1, lnblocks
      if (nodetype(b) == LEAF .or. nodetype(b) == PARENT_BLK) then

! Compare the variable against bounds

         Grid_mark = ( maxval(unk(var, GRID_ILO:GRID_IHI, &
                                       GRID_JLO:GRID_JHI, &
                                       GRID_KLO:GRID_KHI, b)) > &
                       var_bnd_min ) &
                     .AND. &
                     ( minval(unk(var, GRID_ILO:GRID_IHI, &
                                       GRID_JLO:GRID_JHI, &
                                       GRID_KLO:GRID_KHI, b)) < &
                       var_bnd_max )

! If the test passed, Grid_mark the block for refinement.
         if (Grid_mark) then

            if (lrefine(b) < lref) then
               refine(b) = .true.
               derefine(b) = .false.
            else if (lrefine(b) == lref) then
               derefine(b) = .false.
               stay(b) = .true.
            else if (lref <= 0) then
               refine(b) = .true.
            end if

         else

            if (lrefine(b) > lrefine_min) then
               if (.not. stay(b) .and. nodetype(b) /= PARENT_BLK) then
                  derefine(b) = .true.
               end if
            end if

         end if

      end if
   end do
!-------------------------------------------------------------------------------

   return
end subroutine gr_markVarBounds
