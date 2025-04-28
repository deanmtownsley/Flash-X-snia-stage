!!****f* source/Driver/Driver_checkGPUErrorCode
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
!!  Driver_checkGPUErrorCode
!!
!! SYNOPSIS
!!
!!  Driver_checkGPUErrorCode()
!!
!! DESCRIPTION
!!
!!  Check error code recorded on the GPU
!!  and abort the simulation if necessary.
!!
!! ARGUMENTS
!!
!!  None
!!
!! SEE ALSO
!!
!!  Driver_abort_GPU
!!
!!***

subroutine Driver_checkGPUErrorCode()

   use Driver_data, ONLY : dr_gpuErrorCode, dr_gpuErrorMsg
   use Driver_interface, ONLY: Driver_abort

   implicit none

   !$acc update host(dr_gpuErrorCode, dr_gpuErrorMsg)

   if (dr_gpuErrorCode /= 0) then
      if (dr_gpuErrorCode < 0) then
         ! message buffer overflowed
         call Driver_abort(dr_gpuErrorMsg // "...and more")
      else
         call Driver_abort(dr_gpuErrorMsg)
      end if
   end if

end subroutine Driver_checkGPUErrorCode
