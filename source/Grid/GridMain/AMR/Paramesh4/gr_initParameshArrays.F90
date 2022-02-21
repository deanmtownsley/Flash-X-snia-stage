!!****if* source/Grid/GridMain/AMR/Paramesh4/gr_initParameshArrays
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
!!  gr_initParameshArrays
!!
!! SYNOPSIS
!!
!!  call gr_initParameshArrays(logical(IN) :: restart,
!!                             integer(IN) :: xlboundary,
!!                             integer(IN) :: xrboundary,
!!                             integer(IN) :: ylboundary,
!!                             integer(IN) :: yrboundary,
!!                             integer(IN) :: zlboundary,
!!                             integer(IN) :: zrboundary
!!                             )
!!
!! DESCRIPTION
!!
!!  Perform early initialization of some Grid data structures.
!!
!!  This routine prepares the Grid for being filled with
!!  meaningful data.
!!
!! ARGUMENTS
!!
!!   restart -   Is the grid being prepared for initialization with
!!               data from a checkpoint file?
!!   xlboundary - boundary condition type of outer domain boundary in lower X direction.
!!   xrboundary - boundary condition type of outer domain boundary in upper X direction.
!!   ylboundary - boundary condition type of outer domain boundary in lower Y direction.
!!   yrboundary - boundary condition type of outer domain boundary in upper Y direction.
!!   zlboundary - boundary condition type of outer domain boundary in lower Z direction.
!!   zrboundary - boundary condition type of outer domain boundary in upper Z direction.
!!
!! NOTES
!!
!!  This files contains a stub implementation that aborts execution.
!!
!! HISTORY
!!
!!  2022 Created as stub for Flash-X    Klaus Weide
!!***
subroutine gr_initParameshArrays(restart,&
                                     &  xlboundary, xrboundary, &
                                     &  ylboundary, yrboundary, &
                                     &  zlboundary, zrboundary)

   use Driver_interface, only: Driver_abort

   implicit none

#include "Simulation.h"
   logical,intent(IN) :: restart
   integer,intent(IN) :: xlboundary, xrboundary
   integer,intent(IN) :: ylboundary, yrboundary
   integer,intent(IN) :: zlboundary, zrboundary

#if defined(BITTREE)
   call Driver_abort(&
        "Flash-X was configured for a Paramesh4 Grid implementation, &
        &but a proper implementation of gr_initParameshArrays was not&
        & found. Maybe the directory Paramesh4/PM4_package/bittree is&
        &  empty?")
#else
   call Driver_abort(&
        "Flash-X was configured for a Paramesh4 Grid implementation, &
        &but a proper implementation of gr_initParameshArrays was not&
        & found. Maybe a properly checked-out git submodule is missing&
        & in the source directory Paramesh4/PM4_package ?")
#endif   

end subroutine gr_initParameshArrays

