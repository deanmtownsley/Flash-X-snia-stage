!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMoments
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
!!  gr_mpoleMoments
!!
!! SYNOPSIS
!!
!!  gr_mpoleMoments (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution, assuming
!!  the center of mass and the total mass have first been computed. On
!!  output, the gr_mpoleMomentR and gr_mpoleMomentI arrays contain the mass
!!  moments over the regular and irregular solid harmonics. This routine
!!  calls the appropriate subroutines according to the domain geometry specified.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

subroutine gr_mpoleMoments (idensvar)

  use gr_mpoleInterface, ONLY : gr_mpoleMom3Dcartesian,   &
                                gr_mpoleMom2Dcylindrical, &
                                gr_mpoleMom2Dspherical, &
                                gr_mpoleMom3Dspherical, &
                                gr_mpoleMom1Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleGeometry

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

          call gr_mpoleMom3Dcartesian (idensvar)

    case (GRID_3DCYLINDRICAL)

          call  Driver_abort("this geometry is not supported")

    case (GRID_2DCYLINDRICAL)

          call gr_mpoleMom2Dcylindrical (idensvar)

    case (GRID_3DSPHERICAL)

          call gr_mpoleMom3Dspherical(idensvar)

    case (GRID_2DSPHERICAL)

          call gr_mpoleMom2Dspherical(idensvar)

    case (GRID_1DSPHERICAL)

          call gr_mpoleMom1Dspherical (idensvar)

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleMoments
