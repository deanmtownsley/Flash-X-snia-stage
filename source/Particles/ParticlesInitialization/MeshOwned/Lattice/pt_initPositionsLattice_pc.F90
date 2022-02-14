!!****if* source/Particles/ParticlesInitialization/Amrex/Lattice/pt_initPositionsLattice
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
!!    pt_initPositionsLattice
!!
!! SYNOPSIS
!!
!!    call pt_initPositionsLattice(integer(in)  :: tileDesc,
!!                                 logical(OUT) :: success)
!!
!! DESCRIPTION
!!    Initialize particle locations.  This version sets up particles
!!      which are evenly distributed in space along the axes.  Distribution
!!      in non-Cartesian coordinate systems is not regular.
!!
!! ARGUMENTS
!!
!!  tileDesc:        local tile (Grid_tile_t) containing particles to create
!!  success:        returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!    pt_numX:      number of particles along physical x-axis of domain
!!    pt_numY:      number of particles along physical y-axis of domain
!!    pt_numZ:      number of particles along physical z-axis of domain
!!    pt_initialXMin, pt_initialXMax:  physical domain to initialize with particles in X
!!    pt_initialYMin, pt_initialYMax:  physical domain to initialize with particles in Y
!!    pt_initialZMin, pt_initialZMax:  physical domain to initialize with particles in Z
!!    pt_initialRadius: Radius of region about center of the domain to
!!                      which particle initialization is confined.  Ignored
!!                      if negative.
!!
!!***


subroutine pt_initPositionsLattice_pc (ind, tileDesc,success)

  use Particles_interface, ONLY: Particles_mapFromMesh
  use Particles_data, ONLY:  pt_numLocal, particles, pt_maxPerProc, &
       pt_xmin, pt_ymin, pt_zmin,pt_xmax, pt_ymax, pt_zmax,&
       pt_containers, &
       pt_posAttrib,pt_velNumAttrib, pt_velAttrib,pt_typeInfo, pt_meshMe
  use Driver_interface, ONLY : Driver_abort
       
  use Particles_data, ONLY : pt_numX, pt_numY, pt_numZ,&
       pt_initialZMin, pt_initialZMax, &
       pt_initialXMin, pt_initialXMax, pt_initialYMin, pt_initialYMax, &
       pt_initialRadius
  use Particles_data, ONLY : pt_geometry

  use Grid_tile,        ONLY : Grid_tile_t
  use pt_amrexParticleMod, ONLY : pt_amrexParticle_t, &
                                                            amrex_get_next_particle_id, amrex_get_cpu,&
                                                            amrex_set_particle_id, amrex_set_particle_cpu
   use amrex_amr_module, ONLY : amrex_get_amrcore

  implicit none
#include "constants.h"
#include "Simulation.h"
#include "Particles.h"

  integer, intent(in) :: ind
  type(Grid_tile_t)    :: tileDesc
  logical,intent(OUT) :: success

  integer       :: i, j, k, dir
  integer       :: p
  integer       :: blockType,mapType
  logical       :: IsInBlock,IsInSphere
  real          :: xpos, ypos, zpos, rpos, bxl, byl, bzl, bxu, byu, bzu
  real,dimension(MDIM) :: delta
  real, dimension(pt_velNumAttrib) :: partAttribVec
  real, contiguous, pointer :: solnData(:,:,:,:)

! NOTE dxParticle is particle spacing, not grid spacing
  real, dimension(MDIM) :: dxParticle = 0.0
  real, dimension(2,MDIM):: boundBox
  integer :: part_props=NPART_PROPS
  type(pt_amrexParticle_t) :: thisParticle
  integer :: grd_index,tile_index,level_amrex
!----------------------------------------------------------------------

  !       Initialization now done in Particles_init.
  
  !        Particle slot number
  
  p = pt_numLocal
  
  mapType=pt_typeInfo(PART_MAPMETHOD,ind)

  !Block information used by ParticleContainer functions
  level_amrex = tileDesc%level-1
  grd_index=tileDesc%grid_index
  tile_index=tileDesc%tile_index
  
  ! Get grid geometry for this tile 

  call tileDesc%boundBox(boundBox)
  bxl = boundBox(LOW,IAXIS)
  bxu = boundBox(HIGH,IAXIS)
  
  if (NDIM >= 2) then
     byl = boundBox(LOW,JAXIS)
     byu = boundBox(HIGH,JAXIS)
  endif
  if (NDIM == 3) then
     bzl = boundBox(LOW,KAXIS)
     bzu = boundBox(HIGH,KAXIS)
  endif
  
  call tileDesc%deltas(delta)
  nullify(solnData)
  call tileDesc%getDataPtr(solnData, CENTER)

  ! determine particle spacing as dxParticle
  dxParticle = 0.0   ! initialize to zero for all dimensions
  dxParticle(IAXIS) = (pt_initialXMax - pt_initialXMin) / pt_numX
  if (NDIM > 1) dxParticle(JAXIS) = (pt_initialYMax - pt_initialYMin) / pt_numY
  if (NDIM == 3) dxParticle(KAXIS) = (pt_initialZMax - pt_initialZMin) / pt_numZ
  
  
  !! initialization in case of lower dimensionality
  zpos = 0.0
  ypos = 0.0
  xpos = 0.0
  
  loop_x:  do i = 1, pt_numX
     xpos = (i-0.5)*dxParticle(IAXIS) + pt_initialXMin
     IsInBlock = (xpos >= bxl) .and. (xpos < bxu)
     if (.not. IsInBlock) cycle loop_x   !! skip rest of statements if not in block
     
     loop_y:  do j = 1, pt_numY
        if (NDIM >= 2) then
           ypos = (j-0.5)*dxParticle(JAXIS) + pt_initialYMin
           IsInBlock = (ypos >= byl) .and. (ypos < byu)
           if (.not. IsInBlock) cycle loop_y
        endif
        
        loop_z:  do k = 1, pt_numZ
           if (NDIM == 3) then
              zpos = (k-0.5)*dxParticle(KAXIS) + pt_initialZMin  ! location of particle in z
              IsInBlock = (zpos >= bzl) .and. (zpos < bzu)
           endif

!! Restriction to sphere only implemented for CARTESIAN, CYLINDRICAL, or SPHERICAL
!! coordinates.
           
           if ( pt_geometry == CARTESIAN ) then

             rpos = sqrt(xpos**2+ypos**2+zpos**2)

           else if ( pt_geometry == CYLINDRICAL ) then

             rpos = sqrt(xpos**2+ypos**2)

           else if ( pt_geometry == SPHERICAL ) then

             rpos = xpos

           else
           
             rpos = -1.0

           endif
           
           IsInSphere = ( pt_initialRadius <= 0.0 ) .or. ( rpos <= pt_initialRadius )

           if (IsInBlock .and. IsInSphere) then
              p = p + 1
             !! particle is defined, set up data structure
                thisParticle%pos(1) = xpos
#if NDIM > 1
                    thisParticle%pos(2) = ypos
#endif
#if NDIM == 3
                    thisParticle%pos(3) = zpos
#endif

                    !       Initialize velocity properties for the new particles
                    call Particles_mapFromMesh(mapType,pt_velNumAttrib,pt_velAttrib,&
                         (/xpos,ypos,zpos/),&
                         boundBox,delta, tileDesc%limits,solnData,partAttribVec)

                    thisParticle%vel = 0.d0
                    do dir=1,min(NDIM,pt_velNumAttrib)
                       thisParticle%vel(dir) = partAttribVec(dir)
                    end do

!!$                ! The naive code on the following 2 lines does not work any more!
!!$                thisParticle%id  = amrex_get_next_particle_id()
!!$                thisParticle%cpu = amrex_get_cpu()  !DevNote :: or = pt_meshMe

                ! This way of setting id/cpu is necessary since about AMReX 20.09.
                call amrex_set_particle_id(amrex_get_next_particle_id(), thisParticle)
                call amrex_set_particle_cpu(amrex_get_cpu(), thisParticle)

                call pt_containers(ind)%add_particle(level_amrex, grd_index, tile_index, thisParticle)  
              
           endif   !! end of IsInBlock .and. IsInSphere is true
           
        enddo loop_z
     enddo loop_y
  enddo loop_x
  !       Setting the particle database local number of particles
  pt_numLocal = p  
  
  call tileDesc%releaseDataPtr(solnData, CENTER)

  success=.true.

  return

!----------------------------------------------------------------------
  
end subroutine pt_initPositionsLattice_pc


