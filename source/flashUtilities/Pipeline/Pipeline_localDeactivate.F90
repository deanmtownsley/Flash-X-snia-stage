!!****f* source/flashUtilities/Pipeline/Pipeline_localDeactivate
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
!!  Pipeline_localDeactivate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localDeactivate (logical, optional, intent (in) :: doAsyncReturn)
!!
!! DESCRIPTION
!!
!!  Deactivates the local pipeline environment. Deactivation means smooth shutdown
!!  of the local pipeline communitation environment. There are two ways the local
!!  pipeline communicator can shut down: 1) all processors return at will or 2) all
!!  processors return at the same time.
!!
!! ARGUMENTS
!!
!!  doAsyncReturn : if set true, each processor returns at will
!!
!! NOTES
!!
!!  This routine should be called at the end of a pipeline application. Deactivation
!!  is essential before destroying the pipeline. One should never destruct a pipeline
!!  before deactivating it.
!!
!!***

subroutine Pipeline_localDeactivate (doAsyncReturn)

  implicit none

  logical, optional, intent (in) :: doAsyncReturn

  return
end subroutine Pipeline_localDeactivate
