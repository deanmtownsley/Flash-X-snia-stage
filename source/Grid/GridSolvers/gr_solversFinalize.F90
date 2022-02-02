!!****if* source/Grid/GridSolvers/gr_solversFinalize
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
!!  gr_solversFinalize
!!
!! 
!! SYNOPSIS
!!
!!  call gr_solversFinalize()
!!
!!
!! DESCRIPTION
!!
!!  This routine finalizes all solvers; ones in 
!!  use will have something more interesting than a stub
!!  compiled in.
!!
!!***

subroutine gr_solversFinalize()

  use gr_hypreInterface, ONLY: gr_hypreFinalize
  use gr_amrexMultigridInterface, only: gr_amrexMultigridFinalize

  implicit none 

  ! Multipole
  call gr_mpoleFinalize()

  call gr_hypreFinalize()

  call gr_amrexMultigridFinalize()
return
end subroutine gr_solversFinalize
