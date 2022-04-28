!!****if* source/Particles/ParticlesInitialization/WithDensity/CellMassBins/pt_initPositionsWithDensity
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
!!  pt_initPositionsWithDensity
!!
!! SYNOPSIS
!!
!!  call pt_initPositionsWithDensity(integer(IN)  :: blockID,
!!                                   logical(OUT) :: success)
!!
!! DESCRIPTION
!!
!!    Initialize particle locations.  This version sets up passive tracer
!!      particles that are distributed randomly according to the gas density     
!!
!! ARGUMENTS
!!
!!   blockID:  ID of block in current processor
!!
!!   success:       returns .TRUE. if positions for all particles
!!                  that should be assigned to this block have been
!!                  successfully initialized.
!!
!! PARAMETERS
!!
!!  pt_numParticlesWanted:   integer  Approximate number of tracer particles to use
!!                                throughout domain ??
!!
!!***
!========================================================================

subroutine pt_initPositionsWithDensity (tileDesc,success)

  use Particles_data, ONLY:  pt_meshMe, pt_numLocal, &
       pt_maxPerProc, particles, pt_velNumAttrib, pt_velAttrib,pt_typeInfo
  use Particles_data, ONLY : pt_numParticlesWanted, pt_totalMass

  use Grid_interface, ONLY : Grid_getCellVolumes, Grid_getCellCoords
!!$    Grid_getBlkIndexLimits, Grid_getSingleCellVol, Grid_getDeltas, &
!!$    Grid_getBlkBoundBox
  use Driver_interface, ONLY:  Driver_abort
  use Logfile_interface, ONLY:  Logfile_stamp
  use Particles_interface, ONLY : Particles_mapFromMesh
  use pt_interface, ONLY : pt_initLocal
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Particles.h"
  type(Grid_tile_t), INTENT(in) :: tileDesc
  logical, INTENT(out) :: success

  integer          :: newParticlesThisBlock
  real          :: blockMass
  real          :: dvol
  real          :: xpos, ypos, zpos, xvel, yvel, zvel, radius
  integer       :: p, b, bb, i, j, k, numParticlesThisBlock, tag

  real, dimension(:,:,:,:), pointer :: solnData

  real    :: urp  ! random number generator
  integer :: part_props=NPART_PROPS
  integer :: myType=1
  integer :: partID, ind
  integer :: blockList(MAXBLOCKS)
  integer :: iSize,jSize,kSize,numCells,icell
  real :: accu
  integer, dimension(MDIM) :: point  ! indices for a given location
  real, dimension(MDIM) :: delta, pos,partAttribVec
  integer, dimension(LOW:HIGH,MDIM):: blkLimits
  real,dimension(:),allocatable :: xLeftCoord,yLeftCoord,zLeftCoord
  real,dimension(:,:,:),allocatable :: cellVol
  integer,dimension(:,:),allocatable :: icell2ijk
  real,dimension(:),allocatable :: mass, cumuMass
  real,dimension(LOW:HIGH,MDIM) :: bndBox

  !----------------------------------------------------------------
  nullify(solnData)
  ! Access the mesh data for this block.
  
  call tileDesc%getDataPtr(solnData,CENTER)
  call tileDesc%deltas(delta)
  call tileDesc%boundBox(bndBox)
  blkLimits=tileDesc%limits

  iSize = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  jSize = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  kSize = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1
  numCells = iSize*jSize*kSize
  allocate(xLeftCoord(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)))
  allocate(yLeftCoord(blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)))
  allocate(zLeftCoord(blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))
  allocate(cellVol(blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
                   blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
                   blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)))

  !----------------------------------------------------------------------
  call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level,blkLimits(LOW,:),blkLimits(HIGH,:) , xLeftCoord)
  if (NDIM >= 2) then
  call Grid_getCellCoords(JAXIS, LEFT_EDGE, tileDesc%level,blkLimits(LOW,:),blkLimits(HIGH,:) , yLeftCoord)
  end if
  if (NDIM == 3) then
  call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level,blkLimits(LOW,:),blkLimits(HIGH,:) , zLeftCoord)
  end if

  ! Compute the amount of mass in this block and find the maximum density.

  blockMass = 0.
  call Grid_getCellVolumes(tileDesc%level,blkLimits(LOW,:),blkLimits(HIGH,:), cellvol)

! This duplicates the loop below, and could perhaps be consolidated with it.
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     point(3) = k
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        point(2) = j
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           point(1) = i

            dvol = cellVol(i,j,k)
           blockMass = blockMass + solnData(DENS_VAR,i,j,k)*dvol
        enddo
     enddo
  enddo

  numParticlesThisBlock = int(anint(pt_numParticlesWanted * (blockMass / pt_totalMass)))

  !  Check that this requested number isn't going to blow data allocation
  success=.true.
  if (pt_numLocal + numParticlesThisBlock > pt_maxPerProc) then
     call Logfile_stamp(numParticlesThisBlock,"This block mass would generate too many additional particles:")
     call tileDesc%releaseDataPtr(solnData, CENTER) !cleanup before premature return...
     deallocate(cellVol)
     deallocate(zLeftCoord)
     deallocate(yLeftCoord)
     deallocate(xLeftCoord)
     success=.false.
     pt_numLocal=0
     return
  end if
  
  allocate(icell2ijk(MDIM,numCells))
  allocate(mass(numCells))
  allocate(cumuMass(numCells))


  ! Initialize the particles for this block.



  mass = 0.0
  accu = 0.0
  icell = 0
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           icell = icell+1
           icell2ijk(IAXIS,icell) = i
           icell2ijk(JAXIS,icell) = j
           icell2ijk(KAXIS,icell) = k
           mass(icell) = solnData(DENS_VAR,i,j,k) * cellVol(i,j,k)
           accu = accu + mass(icell)
           cumuMass(icell) = accu
        end do
     end do
  end do
  if (icell .NE. numCells) then
     call Driver_abort('pt_initPositionsWithDensity: Garbage icell after accumulation loop!')
  end if
  deallocate(cellVol)

  do p = 1, numParticlesThisBlock


     !  NOTE that some whole complicated (and LBR doesn't understand it) procedure
     !  has been called in Particles_initWithDensity to initialize the random seed.

     call random_number(urp)
     urp = urp * blockMass
     call ut_hunt(cumuMass,numCells,urp,icell)
     icell = icell + 1
     if (icell > numCells) then
        call Driver_abort('pt_initPositionsWithDensity: Garbage after ut_hunt call!')
     end if


     i = icell2ijk(IAXIS,icell)
     j = icell2ijk(JAXIS,icell)
     k = icell2ijk(KAXIS,icell)

     ! Calculate X position within cell, needed in all geometries
     call random_number (harvest=urp)

     xpos = urp * delta(IAXIS)             
     xpos = xpos + xLeftCoord(i)

     if (NDIM >= 2) then
        call random_number (harvest=urp)
        ypos = urp * delta(JAXIS)
        ypos = ypos + yLeftCoord(j)
     else
        ypos = 0.
     endif

     if (NDIM == 3) then
        call random_number (harvest=urp)
        zpos = urp * delta(KAXIS)
        zpos = zpos + zLeftCoord(k)
     else
        zpos = 0.
     endif

     pos(IAXIS)=xpos
     pos(JAXIS)=ypos
     pos(KAXIS)=zpos

        
     ! now set up 
     partID=pt_numLocal+p

     particles(BLK_PART_PROP,partID) = real(tileDesc%id)
     particles(PROC_PART_PROP,partID) = real(pt_meshMe)
     particles(POSX_PART_PROP:POSZ_PART_PROP,partID) = pos(IAXIS:KAXIS)


     ! Now do velocities.  Originally this code (and FLASH2) had xvel=yvel=zvel=0
     ! Alan suggested that it should mimic what Lattice did, and initialize from the mesh
     call Particles_mapFromMesh(pt_typeInfo(PART_MAPMETHOD,myType),pt_velNumAttrib,pt_velAttrib,pos,&
          bndbox,delta, blkLimits,solnData,partAttribVec)
        
     do ind=1,pt_velNumAttrib
        particles(pt_velAttrib(PART_DS_IND,ind),partID) = partAttribVec(ind)
     end do
     
  enddo
  
  !---------------------------------------------------------------------
  ! Cleanup

  ! Release the mesh data pointer for this block.
  call tileDesc%releaseDataPtr(solnData, CENTER)

  ! Set the particle database local number of particles.

  pt_numLocal = pt_numLocal + numParticlesThisBlock

  ! cleanup
  deallocate(icell2ijk)
  deallocate(mass)
  deallocate(cumuMass)
  deallocate(xLeftCoord)
  deallocate(yLeftCoord)
  deallocate(zLeftCoord)

  return

end subroutine pt_initPositionsWithDensity


