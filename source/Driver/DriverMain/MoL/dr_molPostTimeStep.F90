!!****if* source/Driver/DriverMain/MoL/dr_molPostTimeStep
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!      dr_molPostTimeStep
!!
!!  SYNOPSIS
!!
!!      call dr_molPostTimeStep(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Perform any post-timestep work.  This will be called after MoL_advance,
!!      unit-specific evolution calls.  Re-fluxing and any diagnostic calculations
!!      should be performed here.
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine dr_molPostTimeStep(t)
   use Hydro_interface, only: Hydro_molPostTimeStep
   use RadTrans_interface, only: RadTrans_molPostTimeStep
   use Simulation_interface, only: Simulation_molPostTimeStep

   implicit none

   real, intent(in) :: t

   call Hydro_molPostTimeStep(t)
   call RadTrans_molPostTimeStep(t)
   call Simulation_molPostTimeStep(t)
end subroutine dr_molPostTimeStep
