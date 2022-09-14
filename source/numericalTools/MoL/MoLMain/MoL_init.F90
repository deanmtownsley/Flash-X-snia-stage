!!****if* source/numericalTools/MoL/MoLMain/MoL_init
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
!!  NAME
!!
!!      MoL_init
!!
!!  SYNOPSIS
!!
!!      call MoL_init()
!!
!!  DESCRIPTION
!!
!!      Initialize the method of lines unit
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine MoL_init()
   use ml_interface, only: ml_init
   use MoL_data, only: ml_nscratch, ml_nscratch_total, &
                       ml_abortOnWarn, ml_verbosity, ml_mpiRank
   use MoL_variables, only: ml_nvars

   use Driver_interface, only: Driver_getMype
   use RuntimeParameters_interface, only: RuntimeParameters_get

#include "constants.h"

   implicit none

   call Driver_getMype(MESH_COMM, ml_mpiRank)

   call RuntimeParameters_get("MoL_verbosity", ml_verbosity)
   call RuntimeParameters_get("MoL_abortOnWarn", ml_abortOnWarn)

   ml_nscratch = 0
   ml_nvars = 0

   ! Specific integrator setup
   call ml_init()

   ! +2 for MOL_INITIAL & MOL_RHS
   ml_nscratch_total = 2 + ml_nscratch
end subroutine MoL_init
