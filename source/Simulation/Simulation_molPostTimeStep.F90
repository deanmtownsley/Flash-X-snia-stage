!!****f* source/Simulation/Simulation_molPostTimeStep
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
!!      Simulation_molPostTimeStep
!!
!!  SYNOPSIS
!!
!!      call Simulation_molPostTimeStep(real, intent(in) :: t)
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
subroutine Simulation_molPostTimeStep(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Simulation_molPostTimeStep
