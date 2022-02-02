!!****f* source/Driver/Driver_logMemoryUsage
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
!!  Driver_logMemoryUsage
!!
!! SYNOPSIS
!!
!!  Driver_logMemoryUsage(character(len=*)(IN) :: callsite)
!!
!! DESCRIPTION
!!
!!  Logs memory usage
!!
!!
!! ARGUMENTS
!!
!!  callsite -    A string storing from where we call this subroutine.
!!
!! NOTES
!!
!!  This routine is a no-op if the preprocessor symbol FLASH_USE_MEMORYUSAGE is undefined.
!!***

subroutine Driver_logMemoryUsage (callsite)
  implicit none
  character(len=*), intent(IN) :: callsite
end subroutine Driver_logMemoryUsage
