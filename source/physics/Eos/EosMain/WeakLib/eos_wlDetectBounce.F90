!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup physics_Eos
!!
!! @brief  This routine determines if collapse has proceeded to the point of 
!!         core bounce, as determined by the maximum density
!!
!!  @param postBounce  : flag that indicates if bounce has been detected
!!  @param bounceTime  : time at which bounce occurs
!!  @param centralDens : maximum density
!!  @param centralEntr : minimum density inside eos_shockEntrRad
!!


!!REORDER(4): solnData

subroutine eos_wlDetectBounce(postBounce,bounceTime,centralDens,centralEntr)
  !
  !==============================================================================
  !
#include "Simulation.h"
#include "constants.h"

  use Driver_interface, ONLY : Driver_getNStep, Driver_getSimTime
  use Grid_interface, ONLY : Grid_getCellCoords, &
     Grid_getTileIterator, Grid_releaseTileIterator, &
     Grid_coordTransfm
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use IO_interface, ONLY : IO_setScalar
  use Logfile_interface, ONLY : Logfile_stampMessage
  use eos_data, ONLY : eos_meshComm, eos_meshMe
  use eos_wlData, ONLY : eos_postBounce, eos_bounceTime, &
       eos_centralDens, eos_centralEntr, eos_nstep, &
       eos_bounceDens, eos_shockEntr, eos_shockEntrRad

#include "Flashx_mpi_implicitNone.fh"

  logical, intent(OUT) :: postBounce
  real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  integer, dimension(LOW:HIGH,MDIM) :: tileLimits, tileLimitsGC
  real, pointer :: solnData(:,:,:,:)

  integer,dimension(MDIM)  :: dimSize
  real, allocatable, dimension(:) :: xCenter, yCenter, zCenter

  integer :: i, j, k
  integer :: ierr

  character(len=100)  :: message

  real :: radius, theta, phi
  real :: localMaxDens
  real :: localMaxEntr, globalMaxEntr, globalMinEntr
  real :: localMinEntr
  real, dimension(2) :: localMax, globalMax

  integer :: nstep
  real :: time

  if (eos_postBounce) then
     postBounce = .TRUE.
     if (present(bounceTime)) bounceTime = eos_bounceTime
     if (present(centralDens)) centralDens = eos_centralDens
     if (present(centralEntr)) centralEntr = eos_centralEntr
     return !bounce already detected
  end if

  ! Calling routine would like us to check for bounce!
  ! First, verify that we haven't done this already for
  ! this time step.
  call Driver_getNStep(nstep)
  if (nstep == eos_nstep) then
     ! We've already checked for bounce
     postBounce = .FALSE.
     if (present(bounceTime)) bounceTime = 0.0
     if (present(centralDens)) centralDens = eos_centralDens
     if (present(centralEntr)) centralEntr = eos_centralEntr
     return
  end if

  ! We have yet to check for bounce this step
  eos_nstep = nstep

  ! Now proceed with check
  call Driver_getSimTime(time)

  localMaxDens = 0.0
  localMaxEntr = 0.0
  localMinEntr = 100.0

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)
     tileLimits = tileDesc%limits

     allocate(xCenter(tileLimits(LOW,IAXIS):tileLimits(HIGH,IAXIS)))
     allocate(yCenter(tileLimits(LOW,JAXIS):tileLimits(HIGH,JAXIS)))
     allocate(zCenter(tileLimits(LOW,KAXIS):tileLimits(HIGH,KAXIS)))

     call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, tileLimits(LOW,:), tileLimits(HIGH,:), xCenter)
     call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, tileLimits(LOW,:), tileLimits(HIGH,:), yCenter)
     call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, tileLimits(LOW,:), tileLimits(HIGH,:), zCenter)

     do k = tileLimits(LOW,KAXIS), tileLimits(HIGH,KAXIS)
        do j = tileLimits(LOW,JAXIS), tileLimits(HIGH,JAXIS)
           do i = tileLimits(LOW,IAXIS), tileLimits(HIGH,IAXIS)

              localMaxDens = max( localMaxDens, solnData(DENS_VAR,i,j,k) )

              call Grid_coordTransfm( xCenter(i), yCenter(j), zCenter(k), &
                 radius, theta, phi, geometryOut=SPHERICAL )
              if (radius < eos_shockEntrRad) then
                 localMaxEntr = max( localMaxEntr, solnData(ENTR_VAR,i,j,k) )
                 localMinEntr = min( localMinEntr, solnData(ENTR_VAR,i,j,k) )
              end if

           end do
        end do
     end do

     deallocate(xCenter)
     deallocate(yCenter)
     deallocate(zCenter)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  localMax(1:2) = [ localMaxDens, localMaxEntr ]

  call MPI_AllReduce(localMax, globalMax, 2, FLASH_REAL, MPI_MAX, &
       eos_meshComm, ierr)
  call MPI_AllReduce(localMinEntr, globalMinEntr, 1, FLASH_REAL, MPI_MIN, &
       eos_meshComm, ierr)

  eos_centralDens = globalMax(1)
  globalMaxEntr = globalMax(2)
  eos_centralEntr = globalMinEntr

  if (globalMaxEntr > eos_shockEntr .AND. eos_centralDens > eos_bounceDens) then
     ! Bounce
     eos_postBounce = .TRUE.
     eos_bounceTime = time
     call IO_setScalar("postBounce", eos_postBounce)
     call IO_setScalar("bounceTime", eos_bounceTime)
     if (eos_meshMe == MASTER_PE) then
        write(*,*) "Bounce!", time, eos_centralDens/1e14, eos_centralEntr
        write(message,*) "Bounce!", time, eos_centralDens/1e14, eos_centralEntr
        call Logfile_stampMessage(message)
     endif
  endif

  ! return the dummy variables
  postBounce = eos_postBounce
  if (present(bounceTime)) bounceTime = eos_bounceTime
  if (present(centralDens)) centralDens = eos_centralDens
  if (present(centralEntr)) centralEntr = eos_centralEntr

  return
end subroutine eos_wlDetectBounce
