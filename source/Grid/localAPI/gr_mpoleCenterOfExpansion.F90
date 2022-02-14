!!****if* source/Grid/localAPI/gr_mpoleCenterOfExpansion
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

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleCenterOfExpansion
