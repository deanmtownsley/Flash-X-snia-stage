!!****if* source/Grid/localAPI/gr_mpoleCen1Dspherical
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
!!  gr_mpoleCM1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen1Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 1D spherical
!!  geometry. For a 1D spherical problem the center of expansion for the multipoles
!!  is always at the radial domain origin and does not need to be computed. The
!!  following is computed here:
!!
!!                  1) total mass (aborts, if <= 0)
!!                  2) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen1Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCen1Dspherical
