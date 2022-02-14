!!****f* source/Grid/Grid_setInterpValsGcell
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
!!
!! NAME
!!
!!  Grid_setInterpValsGcell
!!
!!
!! SYNOPSIS
!!
!!  call Grid_setInterpValsGcell(logical(IN) :: setval)
!!
!!
!! DESCRIPTION
!!
!! Sets interpolation values for guardcell filling within INS fractional
!! step method.
!!
!!
!! ARGUMENTS
!!
!! setval = .true.  set values 
!!          .false. restore all interpolation-restriction values to quadratic
!!
!! NOTES
!!
!!  The current functionality is specifically for use by the IncompNS unit.
!!
!!  Only implemented for PARAMESH 4.
!!
!! DEV: This should be generalized, to be useful for other code than IncompNS!
!!***

subroutine Grid_setInterpValsGcell(setval)

  implicit none

  logical, intent(IN) :: setval

  return

end subroutine Grid_setInterpValsGcell

