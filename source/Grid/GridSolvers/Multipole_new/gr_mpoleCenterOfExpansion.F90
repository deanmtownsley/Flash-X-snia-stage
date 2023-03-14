!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCenterOfExpansion
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
!!
!!  gr_mpoleCenterOfExpansion
!!
!! SYNOPSIS
!!
!!  gr_mpoleCenterOfExpansion (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of expansion for the multipoles.
!!  This routine is just the wrapper to call the appropriate routine according
!!  to the geometry present.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

subroutine gr_mpoleCenterOfExpansion (idensvar)

  use gr_mpoleData,      ONLY : gr_mpoleGeometry

  use gr_mpoleInterface, ONLY : gr_mpoleCen3Dcartesian,   &
                                gr_mpoleCen2Dcylindrical, &
                                gr_mpoleCen2Dspherical, &
                                gr_mpoleCen3Dspherical, &
                                gr_mpoleCen1Dspherical
  
  use Driver_interface, ONLY : Driver_abort
  implicit none

#include "gr_mpole.h"
  
  integer, intent (in) :: idensvar
!  
!
!     ...Select the appropriate subroutine.
!
!
  select case (gr_mpoleGeometry)

    case (GRID_3DCARTESIAN)

          call gr_mpoleCen3Dcartesian   (idensvar)

    case (GRID_3DCYLINDRICAL)

          call Driver_abort("this geometry is not supported")

    case (GRID_2DCYLINDRICAL)

          call gr_mpoleCen2Dcylindrical (idensvar)

    case (GRID_3DSPHERICAL)

          call gr_mpoleCen3Dspherical   (idensvar)

    case (GRID_2DSPHERICAL)

          call gr_mpoleCen2Dspherical   (idensvar)

    case (GRID_1DSPHERICAL)

          call gr_mpoleCen1Dspherical   (idensvar)

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCenterOfExpansion
