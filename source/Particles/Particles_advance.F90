!!****f* source/Particles/Particles_advance
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
!!  Particles_advance
!!
!! SYNOPSIS
!!
!!  Particles_advance(real(in) :: dtOld,
!!                    real(in) :: dtNew)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the Particle Unit, called from 
!!  within a single block.
!!  Particle velocities are interpolated from the grid.
!!
!!  The two different time steps are needed if a multi-step method changes
!!    timestep sizes between steps. (e.g. Predictor-Corrector or Midpoint)
!!
!!  At the end of moving particles internally within a block, 
!!     Particles_advance calls the routine 
!!     Grid_moveParticles, which redistributes particles that have gone
!!     outside the current block.
!!
!! ARGUMENTS
!!
!!    dtOld:    real            time step for the last velocity scheme (not used in some routines)
!!    dtNew:    real            time step for the current round of sweeps
!!
!! PARAMETERS
!!  
!! NOTES
!!
!!  This routine is called TWICE within each time step.  It is called once after
!!     each split-direction Hydro sweep.
!!
!!***

!===============================================================================

subroutine Particles_advance (dtOld,dtNew)
  
  implicit none
  
  real, INTENT(in) :: dtOld, dtNew
  
  
  return
  
end subroutine Particles_advance

!===============================================================================

