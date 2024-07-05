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
!! @brief  set the flags to indicate which Eos to apply to each cell.
!!         some cells will compute both.
!!
!!
!! @param mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and
!!             MODE_DENS_TEMP as decribed above.
!!
!! @param vecLen   : number of points for each input variable
!!
!! @param eosData  : two dimension array for holding the input and output data
!!                   for applying Eos to a vector
!!
!!  @param eos_hybFlag : flag to select the one or both of the EoS in the hybrid implementation
!!  @param vecB  : the starting point in the vector where Eos is to be applied
!!  @param vecE  : the ending point in the vector where Eos is to be applied
!!


subroutine eos_hybridSetFlag(mode, vecLen, eosData, vecB, vecE, eos_hybFlag)

   use eos_hybridData, only: eos_hybTransitionDensLo, eos_hybTransitionDensHi

   implicit none

#include "constants.h"
#include "Eos.h"
#include "Simulation.h"

   integer, intent(in) :: mode, vecLen
   real, intent(in), dimension(vecLen,EOS_VARS) :: eosData
   integer, intent(in) :: vecB, vecE
   integer, intent(out) :: eos_hybFlag(vecLen)

   integer :: dens
   integer :: k

   ! These integers are indexes into the lowest location in UNK that contain the appropriate variable

   ! Initialize everything to unused flag
   eos_hybFlag = -1

   ! Set the flag for which EoS to use based on density
   do k = vecB, vecE

      ! Call nuclear eos above transition density
      if (eosData(k,EOS_DENS) > eos_hybTransitionDensHi) then

         eos_hybFlag(k) = EOS_WL

         ! Call both eos if between transition densities
      else if (eosData(k,EOS_DENS) <= eos_hybTransitionDensHi .and. &
               eosData(k,EOS_DENS) > eos_hybTransitionDensLo) then

         eos_hybFlag(k) = EOS_HYB

         ! Call helmholtz if below transition density
      else ! eosData(dens+k) <= eos_hybTransitionDensLo

         eos_hybFlag(k) = EOS_HLM

      end if

   end do

end subroutine eos_hybridSetFlag
