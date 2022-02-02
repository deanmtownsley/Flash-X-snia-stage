!!****if* source/Grid/localAPI/gr_mpoleMoments
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
!!  gr_mpoleMoments
!!
!! SYNOPSIS
!!
!!  gr_mpoleMoments (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution,
!!  assuming the center of mass and the total mass have first been
!!  computed.  On output, the Moment_R () and Moment_I () arrays
!!  contain the mass moments over the regular and irregular
!!  solid harmonics. This routine calls the appropriate subroutines
!!  according to the geometry specified.
!!
!! ARGUMENTS
!!
!!  idensvar -- the index of the density variable
!!
!!***

subroutine gr_mpoleMoments (idensvar)

  implicit none
    
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMoments
