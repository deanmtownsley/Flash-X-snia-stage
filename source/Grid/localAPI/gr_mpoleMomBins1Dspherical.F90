!!****if* source/Grid/localAPI/gr_mpoleMomBins1Dspherical
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
!!  gr_mpoleMomBins1Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMomBins1Dspherical (integer (in) :: maxQtype)
!!
!! DESCRIPTION
!!
!!  Compute the multipole moments of the density distribution in 1D spherical
!!  geometry. On output, the gr_mpoleMomentR and gr_mpoleMomentI arrays contain
!!  the mass moments over the regular and irregular solid harmonics. The moments
!!  are evaluated by an outer (threaded) loop over all radial bin types, ensuring
!!  separate evaluation of different radial bin types for all threads.
!!
!! ARGUMENTS
!!
!!  maxQtype : the total number of different radial bin types
!!
!!***

subroutine gr_mpoleMomBins1Dspherical (maxQtype)

  implicit none
  
  integer, intent (in) :: maxQtype

  return
end subroutine gr_mpoleMomBins1Dspherical
