!!****f* source/physics/SolidMechanics/SolidMechanicsMain/SolidMechanics_init
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contrsmutors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distrsmuted under the License is distrsmuted on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  SolidMechanics_init
!!
!!
!! SYNOPSIS
!!
!!  call SolidMechanics_init(LOGICAL(IN) :: restart)
!!
!! ARGUMENTS
!!
!!  restart - restart flag.
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables which typically take values from runtime parameters.
!!  This must be called once by Driver_initAll.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!***
#include "constants.h"
#include "Simulation.h"

subroutine SolidMechanics_init(restart)

   use SolidMechanics_data
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use Driver_interface, ONLY: Driver_getMype, Driver_getNumProcs, &
                               Driver_getComm

   implicit none
   include 'Flashx_mpi.h'
   logical, intent(in) :: restart

   call RuntimeParameters_get("useSolidMechanics", sm_useSolidMechanics)

   if (.NOT. sm_useSolidMechanics) RETURN

   call Driver_getMype(MESH_COMM, sm_meshMe)
   call Driver_getNumProcs(MESH_COMM, sm_meshNumProcs)
   call Driver_getComm(MESH_COMM, sm_meshComm)

end subroutine SolidMechanics_init
