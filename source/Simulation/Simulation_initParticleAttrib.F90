!!****f* source/Simulation/Simulation_initParticleAttrib
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
!!  Simulation_initParticleAttrib
!!
!! SYNOPSIS
!!
!!  call Simulation_initParticleAttrib(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  A stub for initializing additional particle attributes if necessary.
!!
!!  This interface is called during initialization after
!!   o  particles positions,
!!   o  particle velocities and other properties needed to advance
!!      particle trajectories,
!!   o  mass properties for active particles,
!!   o  and particle tags
!!  have all been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!  
!!
!!  By default this does nothing. A typical use would be to initialize
!!  particle properties that were defined in a simulation directory.
!!  It is likely that an implementation only needs to take action
!!  if the restart flag is false.
!!
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Particles_initPositions
!!
!!***

subroutine Simulation_initParticleAttrib(restart)
  implicit none
  logical,intent(in) :: restart

end subroutine Simulation_initParticleAttrib
