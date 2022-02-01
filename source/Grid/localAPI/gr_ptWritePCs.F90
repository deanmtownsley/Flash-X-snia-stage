!!****if* source/Grid/localAPI/gr_ptWritePCs
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
!!  gr_ptWritePCs
!!
!! SYNOPSIS
!!
!!  gr_ptWritePCs(char(IN)    :: fname(:,:),
!!            logical(IN) :: isCheckpoint)
!!
!!  
!! DESCRIPTION 
!!  
!!  This routine writes data in the AMReX's particleContainer object for 
!!  checkpoint or plot
!!
!!  With time integration only a small fraction of particles move out
!!  of a block at any timestep. However, all particles must be examined
!!  to determine if they moved out of their curret block. With refinement
!!  all particles of one block move together to a new block. The logistics
!!  of moving the data between processors is the same in both situations.
!!  Therefore this routine can be used in both modes. 
!! 
!! ARGUMENTS 
!!
!!  fname : Name of the folder in which data will be written
!!
!! isCheckpoint : Is data for checkpoint? See AMReX user guide for details
!!
!!
!! NOTES
!! called from gr_writeData
!!
!! SEE ALSO
!! gr_writeData
!!
!!
!!***
#include "Simulation.h"

subroutine gr_ptWritePCs(filename, is_checkpoint)

    implicit none

    character(len=*), intent(IN) :: filename
    logical,  intent(IN) :: is_checkpoint

end subroutine gr_ptWritePCs
