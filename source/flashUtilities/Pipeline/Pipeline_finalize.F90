!!****f* source/flashUtilities/Pipeline/Pipeline_finalize
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
!!  Pipeline_finalize
!!
!! SYNOPSIS
!! 
!!  call Pipeline_finalize ()
!!
!! DESCRIPTION
!!
!!  Finalizes the Pipeline unit. It deactivates and destroys any pipeline
!!  that is currently in use. The deactivation is done in a synchronous way
!!  to allow all processors to finalize the Pipeline unit at the same time.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  none.
!!
!!***

subroutine Pipeline_finalize ()

  implicit none

  return
end subroutine Pipeline_finalize
