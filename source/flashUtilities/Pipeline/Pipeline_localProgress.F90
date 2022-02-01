!!****f* source/flashUtilities/Pipeline/Pipeline_localProgress
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
!!  Pipeline_localProgress
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localProgress ()
!!
!! DESCRIPTION
!!
!!  This routine can be considered the motor of the pipeline. It makes sure the items
!!  are transported through the channels on the local processor. As such the routine
!!  should be called often to ensure best possible communication progress between the
!!  channels. Progress on a processor can stall, if the item buffer is full.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  All MPI ranks must call this subroutine to avoid possible deadlock. This subroutine
!!  must remain non-blocking.
!!
!!***

subroutine Pipeline_localProgress ()

  implicit none

  return
end subroutine Pipeline_localProgress
