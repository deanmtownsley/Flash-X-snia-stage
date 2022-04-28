!!****if* source/Particles/ParticlesMain/Particles_updateGridVar
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
!!    Particles_updateGridVar
!!
!! SYNOPSIS
!!    Particles_updateGridVar(integer(IN)    :: partProp,
!!                            integer(IN)    :: varGrid,
!!                            integer(IN),optional    :: mode)
!!
!! DESCRIPTION
!!
!!    Updates the given grid variable with data from the given particle
!!    property.
!!
!! ARGUMENTS
!!               partProp:  Index of particle attribute to interpolate onto 
!!                          the mesh
!!               varGrid:   Index of gridded variable to receive interpolated
!!                          quantity
!!               mode:      (Optional) If zero (default), zero varGrid first;
!!                          if present and nonzero, do not zero varGrid first
!!                          but add data from particles to existing grid data.
!!
!! PARAMETERS
!! 
!!***

#include "Simulation.h"
#include "Particles.h"
subroutine Particles_updateGridVar(partProp, varGrid, mode)
  use Grid_interface, only: Grid_mapParticlesToMesh, Grid_sortParticles
  use Particles_data, only: particles, pt_numLocal, pt_maxPerProc, pt_typeInfo
  implicit none

  integer, INTENT(in) :: partProp, varGrid
  integer, INTENT(in), optional :: mode

  logical, parameter :: doSort = .TRUE.
  integer,dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  integer :: mode1
  integer :: ty
  integer :: p0, pn
  logical :: skip
  
  mode1 = 0
  if(present(mode)) mode1 = mode

  if (doSort) then
  ! We need to call Grid_sortParticles to sort by particle type, just for
  ! the case that we have more than one particle type, and to fill the
  ! particlesPerBlk array.
  !  Sort by type, and then within each sort by block.
#ifdef TYPE_PART_PROP
     call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
          pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
#else
     call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
          pt_maxPerProc,particlesPerBlk,BLK_PART_PROP)
#endif
  ! Now update the pt_typeInfo data structure.
  ! We need to do this here, otherwise we cannot be sure that the pt_typeInfo
  ! elements used below are valid.
     call pt_updateTypeDS(particlesPerBlk)
  end if

  do ty=1, NPART_TYPES
    skip = .false.
#if defined(MASS_PART_PROP) && defined(PDEN_VAR) && defined(SINK_PART_TYPE)
    skip = partProp==MASS_PART_PROP .and. varGrid==PDEN_VAR .and. ty==SINK_PART_TYPE
#endif
    if(.not. skip) then
      p0 = pt_typeInfo(PART_TYPE_BEGIN,ty)
      pn = pt_typeInfo(PART_LOCAL,ty)
      call Grid_mapParticlesToMesh( &
        particles(:,p0:p0+pn-1), NPART_PROPS, pn, pt_maxPerProc, &
        partProp, varGrid, mode1, ptInfo=p0-1)
      mode1 = 1 ! stop zeroing out
    end if
  end do
end subroutine Particles_updateGridVar
