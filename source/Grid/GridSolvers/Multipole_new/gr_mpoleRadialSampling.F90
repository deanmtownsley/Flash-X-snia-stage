!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleRadialSampling
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
!!  gr_mpoleRadialSampling
!!
!! SYNOPSIS
!!
!!  gr_mpoleRadialSampling ()
!!
!! DESCRIPTION
!!
!!  Determines the radial sampling for accumulating the moments.
!!  It defines the radial bins into which each moment falls. If the
!!  radial sampling is chosen unwisely, lots of memory for holding
!!  empty moment bins might be wasted. This routine calls the appropriate
!!  subroutines according to the geometry specified.
!!
!!***

subroutine gr_mpoleRadialSampling ()

  use gr_mpoleInterface, ONLY : gr_mpoleRad3Dcartesian,   &
                                gr_mpoleRad2Dcylindrical, &
                                gr_mpoleRad2Dspherical, &
                                gr_mpoleRad3Dspherical, &
                                gr_mpoleRad1Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleGeometry

  implicit none

#include "gr_mpole.h"
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (gr_mpoleGeometry)

    case (GRID_3DCARTESIAN)

          call gr_mpoleRad3Dcartesian   ()

    case (GRID_3DCYLINDRICAL)

          call Driver_abort("this geometry is not supported")

    case (GRID_2DCYLINDRICAL)

          call gr_mpoleRad2Dcylindrical ()

    case (GRID_3DSPHERICAL)

          call gr_mpoleRad3Dspherical   ()

    case (GRID_2DSPHERICAL)

          call gr_mpoleRad2Dspherical   ()

    case (GRID_1DSPHERICAL)

          call gr_mpoleRad1Dspherical   ()

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleRadialSampling
