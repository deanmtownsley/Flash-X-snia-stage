!!****if* source/Grid/GridMain/AMR/Amrex/gr_auxFluxData
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  gr_auxFluxData
!!
!! SYNOPSIS
!!
!!  use gr_auxFluxData
!!
!!  call gr_auxFluxInit()
!!  call gr_auxFluxFree()
!!
!! DESCRIPTION
!!
!!  This module holds some auxiliary arrays that store fluxes.
!!
!!  The arrays here, together with FlashFluxRegisters provided by
!!  AMReX, implement the SPFS (semi-permanent flux storage) in the
!!  case of the GridMain/AMR/Amrex Grid implementation with no
!!  level-wide fluxes.
!!
!! NOTES
!!
!!  The only non-stub implementation provided is for the Amrex Grid
!!  implementation. This may become useful and used for other Grid
!!  implementations, too.
!!  The Paramesh4 Grid implementation has the same functionality with
!!  somwhat different names and source locations. In fact, the existing
!!  Paramesh4 implementation is the origin of this code. It would makes
!!  sense to unify the implementations.
!!
!!
!! HISTORY
!!  2023-10-28 K. Weide  Created stub
!!
!!
!!***

Module gr_auxFluxData
   implicit none

contains
   subroutine gr_auxFluxInit()
      implicit none
   end subroutine gr_auxFluxInit

   subroutine gr_auxFluxFree()
      implicit none
   end subroutine gr_auxFluxFree
end Module gr_auxFluxData
! Local Variables:
! f90-program-indent: 3
! indent-tabs-mode: nil
! End:
