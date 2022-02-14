!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_init
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
!!  Pipeline_init
!!
!! SYNOPSIS
!! 
!!  call Pipeline_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the Pipeline unit. Anything that is common to all pipelines
!!  during an application should go in here.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  none
!!
!!***

subroutine Pipeline_init ()

  use Pipeline_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,            ONLY : Driver_abort

#include "Simulation.h"
#include "constants.h"

  implicit none
!
!
!     ...Get the needed external data.
!
!
  call Driver_getComm        (GLOBAL_COMM,  pl_commGlobal)
  call Driver_getComm        (  MESH_COMM,  pl_commMesh  )

  call RuntimeParameters_get ("basenm",     pl_baseName  )
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_init
