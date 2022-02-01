!!****if* source/Simulation/SimulationMain/YahilLattimerCollapse
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
!!  Eos_GetInnerMostDensity 
!!
!!
!! DESCRIPTION
!!  This routine determines if the density of inner-most-zone 
!!  reaches the maximum density, labled as "postBounce".
!!
!! ARGUMENTS
!!
!!
!!***

!!REORDER(4): solnData

subroutine Eos_getInnerMostDensity(postBounce)

#include "Simulation.h"
#include "constants.h"

  use Eos_data, ONLY : eos_meshMe, eos_meshComm
  use Grid_interface, ONLY :Grid_getTileIterator, &
      Grid_releaseTileIterator, &
      Grid_getMaxRefinement
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Driver_interface, ONLY : Driver_getSimTime
  use IO_interface, ONLY : IO_setScalar
  use Simulation_data, ONLY: sim_postBounce, sim_maxDens
  
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

  implicit none
  include "Flashx_mpi.h"

  logical, intent(OUT) :: postBounce

  integer :: blockCount

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits

  integer,dimension(MDIM)  :: dimSize
  real,allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer :: i,j,k
  integer :: ierr
  logical :: threadBlockList
  logical :: gcell = .true.

  character(len=100)  :: message

  real :: blkMaxDens, blkMaxEntr
  real :: localMaxDens, globalMaxDens
  real :: localMaxEntr, globalMaxEntr
  real :: radius
  real :: eos_centralDens
  real, dimension(2) :: localMax, globalMax

  integer :: imin,imax,jmin,jmax,kmin,kmax

  integer :: level, maxLev
  real :: time

  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc

  nullify(solnData)

  call Grid_getMaxRefinement(maxLev,mode=1)
  call Driver_getSimTime(time)

  blkMaxEntr = 0.
  blkMaxDens = 0.
  localMaxEntr = 0.
  localMaxDens = 0.

  do level = 1, maxLev

    call Grid_getTileIterator(itor, LEAF, level=level, tiling=.FALSE. )
    do while(itor%isValid())
       call itor%currentTile(tileDesc)
  
       blkLimits = tileDesc%limits(LOW:HIGH,1:MDIM)
  
       call tileDesc%getDataPtr(solnData, CENTER)
  
       imin = blkLimits(LOW,IAXIS)
       imax = blkLimits(HIGH,IAXIS)
       jmin = blkLimits(LOW,JAXIS)
       jmax = blkLimits(HIGH,JAXIS)
       kmin = blkLimits(LOW,KAXIS)
       kmax = blkLimits(HIGH,KAXIS)
  
       blkMaxDens = maxval(solnData(DENS_VAR,imin:imax,jmin:jmax,kmin:kmax))
  
       call tileDesc%releaseDataPtr(solnData, CENTER)
    
       if(blkMaxDens>=localMaxDens) localMaxDens = blkMaxDens
       if(blkMaxEntr>=localMaxEntr) localMaxEntr = blkMaxEntr
  
       call itor%next()  
    enddo
    call Grid_releaseTileIterator(itor)

  end do ! level

  localMax(1:2) = (/localMaxDens,localMaxEntr/)

  call MPI_Allreduce(localMax, globalMax, 2, FLASH_REAL, MPI_MAX, &
       eos_meshComm, ierr)

  globalMaxDens = globalMax(1)
  globalMaxEntr = globalMax(2)

  if (globalMaxDens > sim_maxDens) then
     sim_postBounce = .TRUE.
     if (eos_meshMe == MASTER_PE) then
        write(*,'(A,ES12.3,A,ES12.3)') &
        'Max Dens!  time=',time,' Max Dens=', globalMax(1)
        write(message,'(A,ES12.3,A,ES12.3)') &
        'Max Dens!  time=',time,' Max Dens=', globalMax(1)
        call Logfile_stampMessage(message)
     endif
  endif

  postBounce = sim_postBounce

  return
end subroutine Eos_getInnerMostDensity
