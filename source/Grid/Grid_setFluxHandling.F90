!!****f* source/Grid/Grid_setFluxHandling
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
!!  Grid_setFluxHandling
!!
!! SYNOPSIS
!!
!!  call Grid_setFluxHandling(character*(*){IN) :: handling,
!!                        OPTIONAL,integer(OUT) :: status)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   handling - a string indicating the handling requested;
!!              either 'consv_fluxes' or 'consv_flux_densities'.
!!   status - if present, will be set to 0 for success and 1 for failure.
!!
!!
!!***

subroutine Grid_setFluxHandling(handling, status)

  implicit none

  character(len=*),intent(IN) :: handling
  integer,intent(OUT),OPTIONAL :: status

  ! This stub pretends to always succeed:
  if (present(status)) status = 0
  

end subroutine Grid_setFluxHandling
