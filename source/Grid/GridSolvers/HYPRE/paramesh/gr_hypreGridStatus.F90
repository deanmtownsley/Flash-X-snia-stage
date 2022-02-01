!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreGridStatus
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
!!  gr_hypreGridStatus
!!
!!  SYNOPSIS
!!
!!  call gr_hypreGridStatus (integer(IN)  :: blockCount,
!!                           integer(IN), dimension(blockCount) :: blockList,
!!                  OPTIONAL,integer(IN)  :: nvars)
!!
!!  DESCRIPTION 
!!      With AMR mesh verifies if the grid has been modified, if
!!      modified it resets the HYPRE grid object. If called first time 
!!      it sets up the HYPRE grid. 
!!
!! ARGUMENTS
!!
!!   blockCount     : The number of blocks in the list.   
!!   blockList      : The list of blocks on which the solution must be updated.
!!   nvars          : Number of variables, also number of equations, for a
!!                    system. Default is 1.
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!   HYPRE grid is setup only once in UG. 
!!
!!***

#include "Simulation.h"

subroutine gr_hypreGridStatus (blockCount, blockType, nvars)
 
  use gr_hypreData,     ONLY : gr_hypreGridIsSetUp, &
                               gr_hypreNVars,       &
                               gr_hypreNStep
  use gr_hypreLocalInterface,     ONLY : gr_hypreSetupGrid,gr_hypreDestroyGrid
  use Driver_interface, ONLY : Driver_getNStep

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use tree, ONLY: grid_changed
  
  implicit none
  
#include "constants.h" 

  integer, intent(IN):: blockCount
  integer, intent(IN):: blockType
  integer, intent(IN),OPTIONAL :: nvars

  integer :: NStep
  logical :: forceSetup

  call Timers_start ("gr_hypreGridStatus")

  
  if (.not. gr_hypreGridIsSetUp) then 
     call gr_hypreSetupGrid (blockCount, blockType, nvars)
     call Driver_getNStep(gr_hypreNStep)
  else
     if (present(nvars)) then
        forceSetup = ( gr_hypreNVars .NE. nvars )
     else
        forceSetup = ( gr_hypreNVars .NE. 1 )
     end if
     if (grid_changed == 1 .AND. .NOT. forceSetup) then
        call Driver_getNStep(NStep)
        forceSetup = (Nstep /= gr_hypreNStep)
     end if
     if (forceSetup) then
        call gr_hypreDestroyGrid ()
        call gr_hypreSetupGrid (blockCount, blockType, nvars)
        gr_hypreNStep = NStep
     end if
  end if
 
  call Timers_stop ("gr_hypreGridStatus")
 
  return
  
end subroutine gr_hypreGridStatus
