!!****if* source/Particles/ParticlesMain/passive/unitTest/pt_utFakeParticlesAdvance
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
!!  pt_utFakeParticlesAdvance
!!
!! SYNOPSIS
!!
!!  pt_utFakeParticlesAdvance(real(in) :: dtOld,
!!                           real(in) :: dtNew,
!!                           real(in) :: t,
!!                           integer(in) :: ind)
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
!!   ind   -- index of this particle type
!!  
!! PARAMETERS
!!
!!
!!***

!===============================================================================

subroutine pt_utFakeParticlesAdvance (dtOld,dtNew,t, ind)
    
  use Particles_data, ONLY: particles, pt_numLocal, pt_maxPerProc, useParticles, &
       pt_gcMaskForAdvance, pt_posAttrib, pt_velNumAttrib, pt_velAttrib,&
       pt_typeInfo, pt_indexCount, pt_indexList
  use Grid_interface, ONLY : Grid_moveParticles, Grid_fillGuardCells, &
                             Grid_mapMeshToParticles
  use pt_interface, ONLY: pt_utAnaGetNewPosComponents
  
  implicit none

#include "constants.h"  
#include "Simulation.h"
#include "Particles.h"
  
  real, INTENT(in)  :: dtOld, dtNew, t
  integer, INTENT(in) :: ind
  integer       :: i,nstep
  real          :: jumpx,jumpy,jumpz
  integer       :: mapType
  integer,dimension(1) :: index_list
  integer       :: index_count=1
  logical       :: regrid=.false.

!!------------------------------------------------------------------------------
  
  ! Don't do anything if runtime parameter isn't set
!!$      if (.not.useParticles ) return
  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

  ! Prepare guardcell data needed for particle interpolation.
  !
  ! Experimentation (with the old way of advancing Grid particles)
  ! has shown that at least 2 layers of guardcells need to be filled
  ! with updated data for vel[xyz] and density, in order to get the
  ! same results as for a full guardcell fill, when using native grid interpolation. - KW

!!$  call Grid_fillGuardCells( CENTER, ALLDIR, &
!!$       maskSize=NUNK_VARS, mask = pt_gcMaskForAdvance)

     
  ! Update the particle positions.
  do i = 1, pt_numLocal
     call pt_utAnaGetNewPosComponents(particles, pt_maxPerProc, &
          POSX_PART_PROP,POSY_PART_PROP, POSZ_PART_PROP, t, i)
  enddo
     
  ! Obtain the updated particle velocities at the new positions.

  call Grid_mapMeshToParticles(particles,&
       NPART_PROPS, BLK_PART_PROP,pt_numLocal,&
       pt_posAttrib,pt_velNumAttrib,pt_velAttrib,mapType)

  
  ! Put the particles in the appropriate blocks if they've moved off
  call Grid_moveParticles(particles,NPART_PROPS, pt_maxPerProc, pt_numLocal,&
       pt_indexList, pt_indexCount, regrid) 
  
  return
!!------------------------------------------------------------------------------
  
end subroutine pt_utFakeParticlesAdvance


