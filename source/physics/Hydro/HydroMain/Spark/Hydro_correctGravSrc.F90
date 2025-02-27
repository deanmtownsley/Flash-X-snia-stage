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

!> @ingroup HydroSpark
!!
!! @brief Corrects gravitational source terms with n+1 potential
!!        to make source terms second-order accurate
!!
!! @todo Is this file needed? There is no interface for this subroutine.
!<

!!Reorder(4):solnData

subroutine Hydro_correctGravSrc(dt)
  use Grid_interface, ONLY : Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_tile,         ONLY : Grid_tile_t
  use Grid_iterator,     ONLY : Grid_iterator_t
  use Hydro_data, ONLY : hya_grav, hy_useTiling, hy_del
  
  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"

  type(Grid_iterator_t) :: itor
  real,              pointer    :: Uin(:,:,:,:)
  type(Grid_tile_t)     :: tileDesc
  integer :: level
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits,blkLimitsGC,grownLimits
  real,dimension(MDIM) :: deltas
  real,    intent(in) :: dt
  real,dimension(:,:,:,:), pointer :: hy_grav

  integer :: n, i, j, k

  real, dimension(MDIM) :: momOld, momNew
  real :: ekin
  integer,dimension(LOW:HIGH,MDIM)::limits

  call Grid_getTileIterator(itor, LEAF, tiling=.false.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     blkLimits(:,:)=tileDesc%limits
     blkLimitsGC(:,:)=tileDesc%blkLimitsGC
     grownLimits(:,:)=tileDesc%grownLimits
     call tileDesc%deltas(deltas)
     level=tileDesc%level
     call tileDesc%getDataPtr(Uin, CENTER)
     limits=blkLimits
     call hy_rk_getGraveAccel(hy_del,limits,blkLimitsGC)
     hy_grav(1:MDIM,&
          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_grav
     
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
              momOld = Uin(DENS_VAR,i,j,k)*Uin(VELX_VAR:VELZ_VAR,i,j,k)
              momNew = momOld + 0.5*dt*Uin(DENS_VAR,i,j,k)*hy_grav(:,i,j,k)
              Uin(ENER_VAR,i,j,k) = Uin(ENER_VAR,i,j,k) &
                   + 0.5*dt*dot_product(momNew, hy_grav(:,i,j,k))/Uin(DENS_VAR,i,j,k)
              Uin(VELX_VAR:VELZ_VAR,i,j,k) = momNew/Uin(DENS_VAR,i,j,k)
              ! ekin = 0.5*dot_product(Uin(VELX_VAR:VELZ_VAR,i,j,k), &
              !      Uin(VELX_VAR:VELZ_VAR,i,j,k))
              ! print *, (Uin(EINT_VAR,i,j,k) - (Uin(ENER_VAR,i,j,k)-ekin))/Uin(EINT_VAR,i,j,k)
           end do
        end do
     end do
     nullify(hy_grav)
     call tileDesc%releaseDataPtr(Uin,CENTER)
     call itor%next()
  end do !!block loop
  call Grid_releaseTileIterator(itor)
  return
end subroutine Hydro_correctGravSrc
