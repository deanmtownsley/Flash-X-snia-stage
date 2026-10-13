!!****f* source/Driver/Driver_abort_GPU
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
!!  Driver_abort_GPU
!!
!! SYNOPSIS
!!
!!  Driver_abort_GPU(character(len=*)(IN) :: errorMessage)
!!
!! DESCRIPTION
!!
!!  Record error message and flag the error code
!!  during computations on GPU.
!!  The error code will be checked later on host device,
!!  by calling Driver_checkGPUErrors subroutine.
!!
!! ARGUMENTS
!!
!!  errorMessage :    A string to write to the logfile (presumably
!!                    indicating what went wrong).
!!
!! NOTES
!!
!!  This function's implementation never returns control to the caller.
!!
!! SEE ALSO
!!
!!  Driver_checkGPUErrors
!!
!!***

subroutine Driver_abort_GPU(errorMessage)

  implicit none

  character(len=*), intent(in) :: errorMessage

  return
end subroutine Driver_abort_GPU
