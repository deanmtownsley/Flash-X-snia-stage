!!****if* source/Grid/GridSolvers/gr_solversInit
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
!!  gr_solversInit
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversInit()
!!
!!
!! DESCRIPTION
!!
!!  This routine initializes all grid solvers; solvers in use will
!!  have something more exciting than a stub.
!!
!!***

subroutine gr_solversInit()    

  use Grid_interface, ONLY : Grid_setSolverDbgContextInfo
!  use Grid_data, ONLY : gr_isolatedBoundaries
  use gr_hypreInterface, only: gr_hypreInit
  use gr_amrexMultigridInterface, only: gr_amrexMultigridInit

  implicit none 

  ! reset debug context info
  call Grid_setSolverDbgContextInfo()


  ! Multipole
  call gr_mpoleInit()

  call gr_hypreInit()
  call gr_amrexMultigridInit()

  ! solver unit testing
  call gr_solversTestInit()

  return
end subroutine gr_solversInit     
