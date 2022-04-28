!!****if* source/flashUtilities/Pipeline/localAPI/pl_localPipelineSetup
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
!!  pl_localPipelineSetup
!!
!! SYNOPSIS
!! 
!!  call pl_localPipelineSetup ()
!!
!! DESCRIPTION
!!
!!  This routine sets up the local pipeline structure. It determines the number of
!!  channels and the channel processor list on the current processor. The default
!!  structure is based on the current grid structure and this is what is done below.
!!  If the user whishes a pipeline based on a different structure, (s)he needs to
!!  export this routine to her/his personal application (simmulation) unit.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  The pipeline structure is set once: 1) pl_numChannels and 2) pl_procList are defined.
!!  pl_procList is an array and is allocated here.
!!
!!***

subroutine pl_localPipelineSetup ()

  implicit none

  return
end subroutine pl_localPipelineSetup
