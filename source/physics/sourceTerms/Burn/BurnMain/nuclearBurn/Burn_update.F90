!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_update
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
!!  Burn_update
!!
!!
!! SYNOPSIS
!!
!!   call Burn_update(
!!      real, intent(IN OUT)  ::  Uin(:,:,:,:),
!!      integer, intent(IN)   :: loGC(:),
!!      integer, intent(IN)   :: blkLimits(:,:),
!!      real, intent(IN)      :: dt
!!   )
!!
!! DESCRIPTION
!!
!!   Update the solution data using the outputs from the burner.
!!
!! ARGUMENTS
!!
!!   Uin       --
!!   loGC      --
!!   blkLimits --
!!   dt        --
!!
!! NOTES
!!
!!  This subroutine assumes that the nuclear energy generation rate is available
!!             at the Uin(ENUC_VAR,:,:,:) calculated from the burner.
!!
!!***

!!REORDER(4): Uin

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

subroutine Burn_update (Uin, loGC, blkLimits, dt)

  use Eos_interface, ONLY: Eos_multiDim
  use Burn_data, ONLY: bn_nuclearTempMin, &
                       bn_nuclearTempMax, &
                       bn_nuclearDensMin, &
                       bn_nuclearDensMax, &
                       bn_useShockBurn, &
                       bn_nuclearNI56Max

  implicit none

  !args
  integer, dimension(MDIM), intent(IN) ::  loGC
  real, dimension(1:,loGC(IAXIS):, loGC(JAXIS):, loGC(KAXIS):), intent(IN OUT) :: Uin
  integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
  real, intent(IN) :: dt

  ! locals
  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
  integer :: i, j, k
  real :: ek, enuc, ei
  real :: tmp, rho

  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           tmp  = Uin(TEMP_VAR,i,j,k)
           rho  = Uin(DENS_VAR,i,j,k)

           okBurnTemp = (tmp >= bn_nuclearTempMin .AND. tmp <= bn_nuclearTempMax)
           okBurnDens = (rho >= bn_nuclearDensMin .AND. rho <= bn_nuclearDensMax)

           okBurnShock = .true.
#ifdef SHOK_VAR
           okBurnShock = (Uin(SHOK_VAR,i,j,k) <= 0.0 .OR. (Uin(SHOK_VAR,i,j,k) > 0.0 .AND. bn_useShockBurn))
#endif
           if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

              if (NI56_SPEC /= NONEXISTENT) then
                 okBurnNickel = (Uin(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
              else    ! nickel is not even a species in this simulation, so we'll always burn
                 okBurnNickel = .TRUE.
              endif

              if (okBurnNickel) then
                 ek = 0.5e0*(Uin(VELX_VAR,i,j,k)**2 +  &
                             Uin(VELY_VAR,i,j,k)**2 +  &
                             Uin(VELZ_VAR,i,j,k)**2)
                 enuc = dt*Uin(ENUC_VAR,i,j,k)
                 ei = Uin(ENER_VAR,i,j,k) + enuc - ek
#ifdef EINT_VAR
                 Uin(EINT_VAR,i,j,k) = ei
#endif
                 Uin(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                 Uin(EELE_VAR,i,j,k) = Uin(EELE_VAR,i,j,k) + enuc
#endif
              end if
           end if

        end do
     end do
  end do

  call Eos_multidim(MODE_DENS_EI, blkLimits, loGC, Uin)

  return
end subroutine Burn_Update

