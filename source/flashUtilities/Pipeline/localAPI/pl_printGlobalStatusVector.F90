!!****if* source/flashUtilities/Pipeline/localAPI/pl_printGlobalStatusVector
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
!!  pl_printGlobalStatusVector
!!
!! SYNOPSIS
!!
!!  call pl_printGlobalStatusVector (integer, intent (in) :: fileUnit)
!!
!! DESCRIPTION
!!
!!  Prints the global status vector to the file associated with the passed file unit number.
!!
!! ARGUMENTS
!!
!!  fileUnit : the file unit number
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine pl_printGlobalStatusVector (fileUnit)

  implicit none

  integer, intent (in) :: fileUnit

  return
end subroutine pl_printGlobalStatusVector
