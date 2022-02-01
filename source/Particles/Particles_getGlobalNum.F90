!!****f* source/Particles/Particles_getGlobalNum
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
!!  Particles_getGlobalNum
!!
!! SYNOPSIS
!!
!!  Particles_getGlobalNum(integer(OUT)  :: globalNumParticles)
!!                
!! DESCRIPTION 
!!
!!  Returns the global (total) number of particles in the 
!!  simulation after computing it by calling an MPI reduction routine.
!!
!!  Needed for input/output across all processors.
!!
!! ARGUMENTS
!!
!!  globalNumParticles :   integer        number of particles across all processors
!!
!! SIDE EFFECTS
!!
!!  The default implementation writes to the log file a line with the current global
!!  number of particles when first called, and then whenever the number has changed
!!  since the previous time.
!!
!!***

subroutine Particles_getGlobalNum(globalNumParticles)

implicit none
  integer, intent(out)  :: globalNumParticles

  !return 0 particles for stub level
  globalNumParticles = 0

  return
end subroutine Particles_getGlobalNum
