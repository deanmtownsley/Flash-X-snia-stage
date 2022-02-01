!!****if* source/Driver/DriverMain/Driver_finalizeAll
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
!!  Driver_finalizeAll
!!
!! SYNOPSIS
!!  Driver_finalizeAll()
!!
!! DESCRIPTION
!!
!!  Calls all the unit finalize routines
!!  which may need
!!  memory deallocated etc before the run end.
!!  Order does matter.
!!
!!***


subroutine Driver_finalizeAll()

  use Eos_interface, ONLY : Eos_finalize
  use Driver_interface, ONLY : Driver_finalizeSourceTerms
  use RuntimeParameters_interface, ONLY : RuntimeParameters_finalize
  use Multispecies_interface, ONLY : Multispecies_finalize
  use Particles_interface, ONLY : Particles_finalize
  use Grid_interface, ONLY : Grid_finalize
  use Hydro_interface, ONLY : Hydro_finalize
  use Driver_data, ONLY: dr_globalMe, dr_restart
  use Simulation_interface, ONLY : Simulation_finalize
  use IO_interface, ONLY : IO_finalize
  use Timers_interface, ONLY: Timers_finalize
  use Gravity_interface, ONLY: Gravity_finalize
  use IncompNS_interface, ONLY: IncompNS_finalize
  use Multiphase_interface, ONLY: Multiphase_finalize
  use HeatAD_interface, ONLY: HeatAD_finalize
  use Stencils_interface, ONLY: Stencils_finalize
  use ImBound_interface, ONLY: ImBound_finalize
implicit none
#include "mpif.h"

  integer :: ierr
 
  
!!$  call Profiler_finalize()
!!$  
  call RuntimeParameters_finalize()

  call Multispecies_finalize()

  call Driver_finalizeSourceTerms( dr_restart )

  call Grid_finalize()            ! Grid package
 
  call Particles_finalize()       ! Particles
  
  call Hydro_finalize()           ! Hydrodynamics
  
  call Eos_finalize()             ! Equation of State

  call Gravity_finalize()         ! Gravity

  call IncompNS_finalize()        ! IncompNS
 
  call Multiphase_finalize()      ! Multiphase

  call HeatAD_finalize()          ! Heat Advection Diffusion

  call Stencils_finalize()       ! Stencils units

  call ImBound_finalize()         ! Immersed Boundary

  call IO_finalize()

  call Simulation_finalize()

  call Timers_finalize()

  call MPI_Finalize(ierr)

  return
end subroutine Driver_finalizeAll








