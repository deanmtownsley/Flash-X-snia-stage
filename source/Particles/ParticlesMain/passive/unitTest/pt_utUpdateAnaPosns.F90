!!****if* source/Particles/ParticlesMain/passive/unitTest/pt_utUpdateAnaPosns
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
!!  pt_utUpdateAnaPosns
!!
!! SYNOPSIS
!!
!!  pt_utUpdateAnaPosns(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!
!!  Updates particles' POS{X,Y,Z}_PART_PROP and VEL{X,Y,Z}_PART_PROP
!!  properties.
!!
!!  This version just sets the new coordinates and velocities
!!  based on an analytically known solution.
!!
!! ARGUMENTS
!!
!!   dtOld -- not used in this first-order scheme
!!   dtNew -- current time increment
!!   t     -- time for which solution is sought
!!  
!!***

!===============================================================================

subroutine pt_utUpdateAnaPosns (dtOld,dtNew,t)
    
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, useParticles
  use Grid_interface, ONLY : Grid_moveParticles, Grid_fillGuardCells, &
                             Grid_mapMeshToParticles
  use pt_interface, ONLY : pt_utAnaGetNewPosComponents
  
  implicit none

#include "constants.h"  
#include "Simulation.h"
  
  real, INTENT(in)  :: dtOld, dtNew, t
  integer       :: i,nstep
  real          :: jumpx,jumpy,jumpz

!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
!!$      if (.not.useParticles ) return


  ! Update the particle positions.
  do i = 1, pt_numLocal
     call pt_utAnaGetNewPosComponents(particles, pt_maxPerProc, &
          POSANALX_PART_PROP, POSANALY_PART_PROP, POSANALZ_PART_PROP, t, i)
  enddo
     
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utUpdateAnaPosns


