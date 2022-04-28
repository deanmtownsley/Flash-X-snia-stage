!!****if* source/Particles/ParticlesInitialization/WithDensity/RejectionMethod/pt_initLocal
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
!!    pt_initLocal
!!
!! SYNOPSIS
!!
!!    pt_initLocal()
!!
!! DESCRIPTION
!!    Local initialization of  particle locations.  Specific initializations that are
!!      needed only with gas density.  Calculates the total volume and the average density
!!      across all blocks.
!!      Initializes random fields.
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pt_pRand:                integer  Something to do with initial random distribution
!!
!! NOTES
!!
!!   There is a nice description of the Fortran90 random number routines at
!!         http://www.nsc.liu.se/~boein/f77to90/a5.html#section21c
!!
!!***

!===============================================================================

subroutine pt_initLocal ()

  use Particles_data, ONLY:  pt_geometry,pt_meshMe, pt_meshNumProcs,pt_pRand, &
       pt_totalMass, pt_totalVolume, pt_averageDensity, pt_numParticlesWanted,&
       pt_meshComm

  use Driver_interface, ONLY : Driver_abort

  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator,&
                             Grid_getCellVolumes
  use Grid_iterator,       ONLY : Grid_iterator_t
  use Grid_tile,           ONLY : Grid_tile_t

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Flashx_mpi.h"



  ! for random number generator
  integer :: seed_size

  real          :: localMass, localVolume, localDensity
  integer       :: i, j, k, ierr
  integer       :: numPEs
  integer, dimension(MDIM) :: point
  integer, dimension(LOW:HIGH,MDIM):: lims
  real, dimension(:,:,:,:), pointer :: solnData
  real, allocatable, dimension(:,:,:) :: cellVol
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc


  !-------------------------------------------------------------------------------

  ! Runtime Parameters


  ! Currently we only support Cartesian and 2D axisymmetric geometries.
  nullify(solnData)
  if ( (pt_geometry /= CARTESIAN) .and. &
       (.not. ((pt_geometry == CYLINDRICAL) .and. (NDIM == 2))) ) &
       call Driver_abort ("pt_initLocal:  unsupported geometry for with density particle initialization!")

  ! In this routine, we determine the total volume and average density on the
  ! grid and save it.  Note that this will only work correctly if pt_initPositions
  ! has been called after the mesh has been set up OR the DENS_VAR variable contains
  ! accurate zone-average densities.

  ! In axisymmetric geometry we don't worry about factors of 2*pi or quadrants, since
  ! we are only concerned about the relative amount of mass in each block.

  localMass = 0.
  localVolume = 0.
  localDensity = 0.

  ! loop over all local leaf blocks
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     call tileDesc%getDataPtr(solnData, CENTER)
     lims=tileDesc%limits
     ! get dimension limits etc.
     allocate(cellVol(lims(LOW,IAXIS):lims(HIGH,IAXIS),&
                      lims(LOW,JAXIS):lims(HIGH,JAXIS),&
                      lims(LOW,KAXIS):lims(HIGH,KAXIS)))
     call Grid_getCellVolumes(tileDesc%level,lims(LOW,:),lims(HIGH,:),cellVol)
              
     do k = lims(LOW,KAXIS), lims(HIGH,KAXIS)
        do j = lims(LOW,JAXIS), lims(HIGH,JAXIS)
           do i = lims(LOW,IAXIS), lims(HIGH,IAXIS)
              localMass = localMass + solnData(DENS_VAR,i,j,k)*cellVol(i,j,k)
              localVolume = localVolume + cellVol(i,j,k)
           enddo
        enddo
     enddo

     !  release the pointer
     call tileDesc%releaseDataPtr(solnData,CENTER)
     deallocate(cellVol)
     call itor%next()
  enddo  !! of looping over all local leaf blocks
  
  call Grid_releaseTileIterator(itor)
  
  !! get the mass across all processors and store it in data variable pt_totalMass
  call mpi_allreduce(localMass, pt_totalMass, 1, FLASH_REAL, MPI_SUM, &
       pt_meshComm, ierr)
  call mpi_allreduce(localVolume, pt_totalVolume, 1, FLASH_REAL, MPI_SUM, &
       pt_meshComm, ierr)
  
  ! now calculate the average density
  pt_averageDensity = pt_totalMass / pt_totalVolume

  !-------------------------------------------------------------------------------

  ! randomize the initial particle positions

  !  returned value seed_size gives the number of integers the processor uses for the
  !    starting value
  call random_seed(SIZE=seed_size)
  
  !  generates a large (from pt_pRand) integer, in general different for each processor

  i = int(pt_pRand * pt_meshMe) + pt_meshNumProcs

  !  initializes the random number vector with a fixed seed (from i)
  call random_seed(PUT=(/(i, j = 1, seed_size)/))
  

  !We used to call random_number lots of times but this does not serve
  !any useful purpose for any of our FLASH simulations.


  !-------------------------------------------------------------------------------

  return

end subroutine pt_initLocal


