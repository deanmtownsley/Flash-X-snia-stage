!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_burner
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
!!  Burn_burner
!!
!!
!! SYNOPSIS
!!
!!   call Burn_burner(
!!   )
!!
!! DESCRIPTION
!!
!!   Apply burner to all blocks and save the results in the solution data.
!!
!! ARGUMENTS
!!
!!
!! NOTES
!!
!!
!!***

!!REORDER(4): Uin

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

subroutine Burn_burner(Uin, loGC, blkLimits, dt)

   use Burn_data, ONLY: bn_nuclearTempMin, &
                        bn_nuclearTempMax, &
                        bn_nuclearDensMin, &
                        bn_nuclearDensMax, &
                        bn_useShockBurn, &
                        bn_nuclearNI56Max
   use bn_interface, ONLY : bn_mapNetworkToSpecies, &
                            bn_burner
   implicit none

   ! input args
   integer, dimension(MDIM), intent(IN) ::  loGC
   real, dimension(1:, loGC(IAXIS):, loGC(JAXIS):, loGC(KAXIS):), intent(IN OUT) :: Uin
   integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
   real, intent(IN) :: dt

   ! local vars
   logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
   integer :: i, j, k, n, speciesMap, lo(MDIM), hi(MDIM)
   real, dimension(NSPECIES) :: xIn, xOut
   real :: tmp, rho, sdot


   lo(1:MDIM) = blkLimits(LOW,1:MDIM)
   hi(1:MDIM) = blkLimits(HIGH,1:MDIM)

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS)

            tmp  = Uin(TEMP_VAR,i,j,k)
            rho  = Uin(DENS_VAR,i,j,k)
            sdot = Uin(ENUC_VAR,i,j,k)

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

                  ! burnedZone = .TRUE.

                  ! Map the solution data into the order required by bn_burner
                  do n = 1, NSPECIES
                     call bn_mapNetworkToSpecies(n,speciesMap)
                     xIn(n) = Uin(speciesMap,i,j,k)
                  end do

                  ! Do the actual burn
                  call bn_burner(dt, tmp, rho, xIn, xOut, sdot)

                  !  Map it back NOTE someday make a nicer interface....
                  do n = 1, NSPECIES
                     call bn_mapNetworkToSpecies(n,speciesMap)
                     Uin(speciesMap,i,j,k) = xOut(n)
                  end do

               endif
            endif

            Uin(ENUC_VAR,i,j,k) = sdot

         end do
      end do
   end do

end subroutine Burn_burner
