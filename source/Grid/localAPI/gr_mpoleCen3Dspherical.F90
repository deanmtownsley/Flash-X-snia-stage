!!****if* source/Grid/localAPI/gr_mpoleCen3Dspherical
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
!!  gr_mpoleCen3Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleCen3Dspherical (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of multipole expansion for 3D spherical
!!  geometry. It computes the center of expansion for the multipoles for 3D spherical
!!  geometries. The center is calculated using the position weighted by the square
!!  density:
!!
!!
!!                          integral (r * rho * rho  dr)
!!                Cen (z) = ----------------------------
!!                            integral (rho * rho  dr)
!!
!!
!!  which, due to uniform density in each cell, becomes:
!!
!!
!!                   sum cells (cell center r * cell mass * cell rho)
!!         Cen (z) = ------------------------------------------------
!!                           sum cells (cell mass * cell rho)
!!
!!
!!  After the initial Cen (z) has been determined, it is placed on the nearest
!!  nearest cell corner. The following is computed here:
!!
!!                  1) multipole expansion center (placed on nearest cell corner)
!!                  2) total mass (aborts, if <= 0)
!!                  3) the 'atomic' inner zone length (and its inverse)
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!! NOTES
!!
!!  gr_mpoleZcenter denotes the location of the center of multipole expansion
!!  in the 3D cartesian version of the 3D spherical framework.
!!  Due to symmetry, there is no need to evaluate and store a x-coordinate
!!  component, which is always equal to zero.
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleCen3Dspherical (idensvar)

  implicit none

  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCen3Dspherical
