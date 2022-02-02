!!****f* source/flashUtilities/Pipeline/Pipeline_localActivate
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
!!  Pipeline_localActivate
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localActivate ()
!!
!! DESCRIPTION
!!
!!  Activates locally the pipeline created. Activation consists in posting speculative
!!  receive messages, awaiting for items to be received on the local processor.
!!  If no channels were created, no communication will happen on the local processor,
!!  but the processor is still considered to be active inside the pipeline.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  It is always a good idea to check beforehand, if the created global pipeline
!!  structure is ok by calling the routine 'Pipeline_globalCheckStructure'.
!!  Otherwise sends might be posted on channels which have no receive channel
!!  counterpart on another processor.
!!
!!***

subroutine Pipeline_localActivate ()

  implicit none

  return
end subroutine Pipeline_localActivate
