!!****f* source/Particles/Particles_init
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
!!    Particles_init
!!
!! SYNOPSIS
!!    Particles_init( logical(in) :: restart )
!!
!! DESCRIPTION
!!
!!    General initialization routine for the particle module.
!!
!! ARGUMENTS
!!
!!    restart:   indicates if run is starting from scratch or restarting
!!               from checkpoint file
!!
!! PARAMETERS
!!
!!    useParticles   BOOLEAN [TRUE]  Should particles be used in this simulation?
!!    pt_maxPerProc  INTEGER [100]   Maximum number of particles per processor. Allocates array space
!!                                   Particles are distributed per PROCESSOR rather than per BLOCK
!!    pt_dtFactor    REAL    [0.5]   Factor to make sure that time step is small enough that particles
!!    pt_dtChangeTolerance REAL [0.4] For uncorrected Estimated Midpoint propagation scheme:
!!                                    Do Euler step if change in time step is greater than this
!!                                    percentage.  Set to 0 to always do Euler, set to a huge
!!                                    number to always use estimated midpoint velocities
!!    pt_small       REAL    [1.0E-10] Used for general comparisons of real values 
!!                                   For example, IF (abs(real1 - real2) .lt. pt_small) THEN
!!                                   don't move farther than one block in each step
!!
!!***
  
subroutine Particles_init( restart)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abort

  implicit none

  logical, INTENT(in) :: restart

  logical, save :: testUseParticles

  !! It is a failure to invoke the stub when useParticles is set TRUE.

  call RuntimeParameters_get ("useParticles", testUseParticles)
  if (testUseParticles) then
     call Driver_abort("Particles unit seems not to be compiled in, and the Particles_init stub does not &
          &allow the value of useParticles to be TRUE.")
  end if

end subroutine Particles_init

