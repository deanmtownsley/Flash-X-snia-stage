!!****f* source/Driver/DriverMain/GPU/Driver_abort_GPU
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

   use Driver_data, ONLY : dr_gpuErrorCode, dr_gpuErrorMsg, dr_gpuErrorMsgStart

   implicit none

   character(len=*), intent(in) :: errorMessage

   integer :: i, n
   integer :: truncate_len, available_buffer
   logical :: buffer_overflow

   !$acc routine

   i = dr_gpuErrorMsgStart
   n = len(errorMessage) + 2
   available_buffer = len(dr_gpuErrorMsg) - i + 1

   if (available_buffer > 0) then
      if (n <= available_buffer) then
         ! there is enough buffer to store the full message
         dr_gpuErrorMsg(i:i+n-1) = errorMessage // ", "
         dr_gpuErrorMsgStart = dr_gpuErrorMsgStart + n
      else
         ! not enough buffer available; truncate the message and add a ellipsis
         truncate_len = available_buffer - 3     ! reserve space for "..."
         if (truncate_len > 0) then
            dr_gpuErrorMsg(i:i+truncate_len-1) = errorMessage(1:truncate_len) // "..."
         else if (available_buffer >= 3) then
            ! just add an ellipsis
            dr_gpuErrorMsg(i:i+2) = "..."
         end if
         dr_gpuErrorMsgStart = len(dr_gpuErrorMsg) + 1
      end if
   else
      buffer_overflow = .true.
   end if

   dr_gpuErrorCode = 1

   ! if the message buffer is overflow, negate the error code
   if (buffer_overflow) dr_gpuErrorCode = -1 * dr_gpuErrorCode

end subroutine Driver_abort_GPU
