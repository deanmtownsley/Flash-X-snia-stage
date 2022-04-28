!!****f* source/Driver/Driver_abort
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
!!  Driver_abort
!!
!! SYNOPSIS
!!
!!  Driver_abort(character(len=*)(IN) :: errorMessage)
!!
!! DESCRIPTION
!!
!!  Write an error message to the logfile and abort FLASH.
!!  Attempts to shut down all processes (using MPI_Abort()).
!!  If you wish to call Driver_abort from a 'c' routine
!!  use the API routine Driver_abortC
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
!!
!!***

subroutine Driver_abort (errorMessage)
  
  implicit none

  character(len=*), intent(in) :: errorMessage

  return
end subroutine Driver_abort
