!!****f* source/Particles/Particles_manageLost
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
!!    Particles_manageLost
!!
!! SYNOPSIS
!!    Particles_manageLost( integer(in) :: mode )
!!
!! DESCRIPTION
!!
!!    This routine is devised to facilitate keeping around particles
!!    even after they leave the physical domain. A runtime parameter
!!    "pt_keepLostParticles", when true triggers this behavior.
!!    These particles need special handling because they should not advance
!!    in time, or be mapped to or from the grid. We do this by assigning a value
!!    "LOST" (defined in constants.h) to the associated block number, and moving
!!    them to the end of the particles data structure. They are not
!!    counted in pt_numLocal, which is the local particles count, when
!!    processing (time advance, mapping or force computations), but
!!    need to be temporarily included for IO. This is the routine that
!!    allows the temporary inclusion by appending "LOST" particles
!!    at the end of the active particles list before writing, and moving
!!    them back at the end of the data structure after writing. 
!!
!! ARGUMENTS
!!
!!  mode -- indicates whether the "LOST" particles are appended at the end
!!          of active list of moved away from there
!!
!! EXAMPLE -- if pt_maxPerProc=1000 (the size of the particles data structure on
!!            a processor), count of particles still in the domain is 250,
!!            and the number of particles that have left the domain is 20, then
!!            for processing pt_numLocal=250, and the indices 981:1000 in the
!!            particles data structure {particles(:,981:1000)} hold the lost
!!            particles with pt_numLost=20. When this routine is invoked
!!            with pt_keepLostParticles=.true. then the following operations are
!!            performed :
!!            if mode==PART_EXPAND
!!                particles(:,251:270)=particles(:,981:1000); pt_numLocal=270
!!            if mode==PART_COLLAPSE
!!                particles(:,981:1000)=particles(:,251:270); pt_numLocal=250
!!            When pt_keepLostParticles=.false. the routine returns without
!!            performing any operations.
!! 
!!***

subroutine Particles_manageLost(mode)
  implicit none
  integer, intent(IN) :: mode
end subroutine Particles_manageLost
