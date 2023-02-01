!!****if* source/Grid/localAPI/gr_mpoleMom2Dspherical
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
!!  gr_mpoleMom2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom2Dspherical (integer (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Prepares for evaluation of the moments in 2D spherical geometry. In this
!!  routine, all the necessary arrays are prepared to enable evaluation of
!!  the moments in radial bin order. Each of the moments are grouped together
!!  according to their radial bins. This will ensure optimum unit stride values
!!  when accessing the big moment arrays and makes threading trivial.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom2Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMom2Dspherical
