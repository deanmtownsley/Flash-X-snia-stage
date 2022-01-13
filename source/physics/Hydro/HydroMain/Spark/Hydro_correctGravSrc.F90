!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_correctGravSrc
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!! NAME
!!
!!  Hydro_correctGravSrc
!!
!!
!! SYNOPSIS
!!                           
!!  Hydro_correctGravSrc(real(IN)    :: dt)
!!
!!
!! DESCRIPTION
!!
!!  Corrects gravitational source terms with n+1 potential to make source
!!  source terms second-order accurate
!!
!! ARGUMENTS
!!
!!  blockCount - the number of blocks in blockList
!!  blockList  - array holding local IDs of blocks on which to advance
!!  dt         - timestep
!!
!!***
!!Reorder(4):solnData
subroutine Hydro_correctGravSrc(dt)
  use Grid_interface, ONLY: Grid_getTileIterator,Grid_releaseTileIterator
  use Hydro_data, ONLY : hy_grav, hy_useTiling
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  
  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)   :: blockDesc

  real,    intent(in) :: dt

  integer :: n, i, j, k
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  real, pointer :: solnData(:,:,:,:)
  real, dimension(MDIM) :: momOld, momNew
  real :: ekin

  call Grid_getTileIterator(itor, LEAF, tiling=hy_useTiling)
  do while(itor%isValid())
     call itor%currentTile(blockDesc)

     blkLimits(:,:)   = blockDesc%limits
     blkLimitsGC(:,:) = blockDesc%blkLimitsGC
     call blockDesc%getDataPtr(solnData,CENTER)
     call hy_rk_getGravAccel(blockDesc,blkLimitsGC)

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              momOld = solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR:VELZ_VAR,i,j,k)
              momNew = momOld + 0.5*dt*solnData(DENS_VAR,i,j,k)*hy_grav(:,i,j,k)
              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) &
                   + 0.5*dt*dot_product(momNew, hy_grav(:,i,j,k))/solnData(DENS_VAR,i,j,k)
              solnData(VELX_VAR:VELZ_VAR,i,j,k) = momNew/solnData(DENS_VAR,i,j,k)
              ! ekin = 0.5*dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k), &
              !      solnData(VELX_VAR:VELZ_VAR,i,j,k))
              ! print *, (solnData(EINT_VAR,i,j,k) - (solnData(ENER_VAR,i,j,k)-ekin))/solnData(EINT_VAR,i,j,k)
           end do
        end do
     end do
     call blockDesc%releaseDataPtr(solnData,CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  return
end subroutine Hydro_correctGravSrc
