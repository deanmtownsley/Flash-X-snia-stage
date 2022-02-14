!!****f* source/flashUtilities/Pipeline/Pipeline_localDestroy
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
!!  Pipeline_localDestroy
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localDestroy ()
!!
!! DESCRIPTION
!!
!!  Destroys the local pipeline structure. Destruction means that all arrays
!!  (the plumbing of the pipeline) are deallocated. Also the local log file
!!  will be closed (if present). 
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This routine should onle be called after the local pipeline has been
!!  deactivated.
!!
!!***

subroutine Pipeline_localDestroy ()

  implicit none

  return
end subroutine Pipeline_localDestroy
