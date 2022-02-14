!!****if* source/Particles/ParticlesInitialization/Particles_initPositions
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
!!    Particles_initPositions
!!
!! SYNOPSIS
!!
!!   call Particles_initPositions( logical(inout) :: partPosInitialized,
!!                                 logical(out)   :: updateRefine)
!!
!!
!! DESCRIPTION
!!
!!    Initialize particle locations. This routine calls pt_initPositions
!!    which a Particles unit's local routine to initialize the positions
!!    on leaf blocks. The routine also creates tags for all the particles
!!    This routine will initialize based on Lattice or with Density 
!!    distribution depending upon which of the two is selected. 
!!
!! ARGUMENTS
!!
!!  partPosInitialized : boolean indicating whether particles positions were 
!!            successfully initialized.
!!            On entry, a value of TRUE is taken to mean that position
!!            initialization has already been completed previously,
!!            so the routine returns immediately (leaving partPosInitialized
!!            TRUE).
!!            If particles are disabled (as per runtime parameter useParticles),
!!            this implementation also returns immediately with
!!            partPosInitialized set to TRUE.
!!            Otherwise, partPosInitialized will be to TRUE if all particles
!!            have been placed in the domain successfully. A return value
!!            of FALSE may indicate that only some particles have been placed
!!            in the domain, perhaps because of space limitations in the
!!            particles array in some MPI tasks; partially initialized
!!            particles data of this kind may still be useful during FLASH
!!            initialization, in particular for the the purpose of providing
!!            refinement criteria for the initial Grid construction if the
!!            runtime parameter refine_on_particle_count is TRUE, but a
!!            fully successful Particles_initPositions invocation is still
!!            required before the simulation is allowed to proceed with
!!            its main evolution loop.
!!
!!  updateRefine : is set to TRUE if the routine wishes to indicate that during
!!                 the initial iterative construction of an AMR Grid (see
!!                 Grid_initDomain and gr_expandDomain), the initialization of
!!                 particle positions need not be repeated for each iteration if
!!                 all particles have already been placed in the domain in a
!!                 previous iteration.  Under that condition, subsequent calls
!!                 to Particles_initPositions from the Grid construction loop
!!                 will have partPosInitialized=.TRUE.  so will return
!!                 immediately, and Particles_updateRefinement will be called in
!!                 each iteration when the number of blocks or the block
!!                 distribution may have changed, to make sure that the retained
!!                 particles get moved to the correct block (hence the name of
!!                 the dummy argument).
!!
!!                 This implementation always returns FALSE.
!!                 Alternative implementations may wish to return TRUE
!!                 instead if initialization is very expensive.
!!                 
!!
!! NOTES
!!
!!  
!!
!! SIDE EFFECTS
!!
!!  Updates particles data that is private to the Particles unit.
!!
!!  May modify an internal flag (pt_posInitialized) that keeps track of
!!  whether initialization of particle positions is complete.
!!
!! SEE ALSO
!!
!!  Driver_initAll
!!  Grid_initDomain
!!  Particles_initData
!!***

!!#define DEBUG_PARTICLES

subroutine Particles_initPositions (partPosInitialized,updateRefine)


  use Grid_interface, ONLY : Grid_getTileIterator, &
                             Grid_releaseTileIterator
  use Driver_interface, ONLY : Driver_abort
  use pt_interface, ONLY : pt_initPositions,pt_createTag, &
                           pt_initLocal,                  &
                           pt_initPositionsLattice,       &
                           pt_initPositionsWithDensity

  use Particles_data, ONLY : pt_posInitialized,pt_numLocal,useParticles,&
       pt_typeInfo, particles, pt_meshNumProcs, pt_meshMe

  use Grid_iterator,       ONLY : Grid_iterator_t
  use Grid_tile,           ONLY : Grid_tile_t

  implicit none
#include "constants.h"
#include "Simulation.h"
#include "Particles.h"

  logical, INTENT(INOUT) :: partPosInitialized
  logical, INTENT(OUT) :: updateRefine

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  
  integer       :: i, j, k, b
  integer       :: p
  integer       :: numLocalThisType, numNewLocalThisType, numLocalPreviousTypes
  integer       :: numPreviousLocal

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  integer :: blkCount
!!  integer,dimension(MAXBLOCKS) :: blkList
!----------------------------------------------------------------------

  if(.not.useParticles) then
     partPosInitialized = .true.
  end if
  if(partPosInitialized) return

  !CD: We need to move the call to pt_initLocal to this level. 
  !Otherwise, if it is in pt_initPositions, we get a deadlock
  !when the number of blocks are not the same on all processors.
  call pt_initLocal()

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  ! Distribute initial positions


  updateRefine = .FALSE.

  partPosInitialized=.true.

  numLocalPreviousTypes = 0
  if(.not.updateRefine) then
     pt_numLocal = 0
  else
     pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
  end if
  do i = 1,NPART_TYPES
     if(.not.updateRefine) pt_typeInfo(PART_LOCAL,i) = 0
     numLocalThisType = pt_typeInfo(PART_LOCAL,i)
     
     !!     b=0
     call Grid_getTileIterator(itor,LEAF,tiling=.FALSE.)
     numNewLocalThisType = 0
     numPreviousLocal = pt_numLocal
     do while(itor%isValid().and.partPosInitialized)
        call itor%currentTile(tileDesc)
        select case(pt_typeInfo(PART_INITMETHOD,i))
        case(LATTICE)
           call pt_initPositionsLattice(tileDesc,partPosInitialized)
        case(WITH_DENSITY, CELLMASS, REJECTION)
           call pt_initPositionsWithDensity(tileDesc,partPosInitialized)
        case(CUSTOM)
           call pt_initPositions(tileDesc,partPosInitialized)
        case default
           call Driver_abort("Particles_initPosition: no valid initialization method")
        end select
        numNewLocalThisType = pt_numLocal - numPreviousLocal
        pt_typeInfo(PART_LOCAL,i) = numNewLocalThisType + numLocalThisType
        call itor%next()
     enddo
     call Grid_releaseTileIterator(itor)
#ifdef TYPE_PART_PROP
     particles(TYPE_PART_PROP, &
          pt_numLocal-numNewLocalThisType+1:pt_numLocal) = pt_typeInfo(PART_TYPE,i) 
#endif
     numLocalThisType=pt_typeInfo(PART_LOCAL,i)
     numLocalPreviousTypes = numLocalPreviousTypes + numLocalThisType
  end do
  pt_numLocal=sum(pt_typeInfo(PART_LOCAL,1:NPART_TYPES))
  !!print*,'Particles_initPositions: pt_numLocal now is',pt_numLocal

  pt_posInitialized = partPosInitialized

#ifdef DEBUG_PARTICLES
  if (pt_meshMe == MASTER_PE .OR. pt_meshNumProcs .LE. 4) then
     print*,'Particles_initPositions on processor', pt_meshMe, 'done, pt_numLocal=',pt_numLocal
  end if
#endif

  call pt_createTag()
  
  return
  
end subroutine Particles_initPositions
