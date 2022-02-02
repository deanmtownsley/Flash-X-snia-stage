!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/Driver_finalizeAll
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
!!  Driver_finalizeAll
!!
!! SYNOPSIS
!!
!!  Driver_finalizeAll ()
!!
!! DESCRIPTION
!!
!!  Stripped down version from original for testing single units.
!!
!!***
subroutine Driver_finalizeAll ()

  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Simulation_interface,        ONLY : Simulation_finalize
  use Timers_interface,            ONLY : Timers_finalize

  implicit none

#include "mpif.h"

  integer :: ierr
 
  call RuntimeParameters_finalize ()
  call Simulation_finalize ()
  call Timers_finalize ()
  call MPI_Finalize (ierr)

  return
end subroutine Driver_finalizeAll
