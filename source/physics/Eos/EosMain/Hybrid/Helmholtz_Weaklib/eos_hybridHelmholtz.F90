!!****if* source/physics/Eos/EosMain/hybrid/Helmholtz_WeakLib/eos_hybridHelmholtz
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  eos_hybridHelmholtz
!!
!! SYNOPSIS
!!
!!  call eos_hybridHelmholtz(integer(IN)  :: mode,
!!                         integer(IN)  :: vecLen,
!!                         real(INOUT)  :: eosData(vecLen*EOS_NUM),
!!               optional, real(IN)     :: massFrac(vecLen*NSPECIES),
!!               optional, logical(IN)  :: mask(EOS_VARS+1:EOS_NUM))
!!
!! DESCRIPTION
!!
!!  This procedure prepares and calls the Helmholtz Eos
!!
!! ARGUMENTS
!!
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points for each input variable
!!
!!  eosData  : This array is the data structure through which variable values are
!!             passed in and out of the Eos routine. The arrays is sized as
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on.
!!
!!  massFrac : If present, contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
!!             If you are working with electron abundance mass scalars, then you
!!             do not necessarily have to have mass fractions.
!!             Thus, for the Ye formulation this value is never used!
!!
!!  mask     : Mask is a logical array the size of EOS_DERIVS (number
!!              of partial derivatives that can be computed, defined in
!!              Eos.h), where each index represents a specific partial derivative
!!              that can be calculated by the Eos unit. A .true. value in mask
!!              results in the corresponding derivative being calculated and
!!              returned. It should preferably be dimensioned as
!!              mask(EOS_VARS+1:EOS_NUM) in the calling routine
!!              to exactly match the arguments declaration in Eos Unit.
!!             Note that the indexing of mask does not begin at 1, but rather at one past
!!             the number of variables.
!!
!!***
subroutine eos_hybridHelmholtz(mode, eosData, mask, massFrac)

   use eos_localInterface, only: eos_helmYe
   use eos_hybridData, only: ergPerKg_to_kbPerBaryon, kbPerBaryon_to_ergPerKg

   !! TODO: port this from BANG
   ! use Burn_interface, only: Burn_computeEbin

   implicit none

#include "constants.h"
#include "Eos.h"
#include "Simulation.h"

   integer, intent(in) :: mode
   real, dimension(EOS_NUM), target, intent(inout) :: eosData
   logical, dimension(EOS_VARS + 1:EOS_NUM), target, intent(in) :: mask
   real, dimension(NSPECIES), target, intent(in), optional :: massFrac

   real :: ebin
   logical :: pMassFrac

   integer :: k

   ebin = 0.0

   pMassFrac = present(massFrac)

   !! TODO: uncomment after porting Burn_computeEbin from BANG
   ! if (pMassFrac) then
   !    call Burn_computeEbin(massFrac, ebin)
   ! end if

   ! Shift internal energy back up by binding energy
   eosData(EOS_EINT) = eosData(EOS_EINT) + ebin

   ! Nuclear Eos gives entropy in k_B per baryon
   if (mode .eq. MODE_DENS_ENTR) &
      eosData(EOS_ENTR) = eosData(EOS_ENTR)*kbPerBaryon_to_ergPerKg

   if (pMassFrac) then
      call eos_helmYe(mode, 1, eosData, mask=mask, massFrac=massFrac)
   else
      call eos_helmYe(mode, 1, eosData, mask=mask)
   end if

   ! Take binding energy back off for NSE consistency
   eosData(EOS_EINT) = eosData(EOS_EINT) - ebin

   ! Back to nuclear units
   eosData(EOS_ENTR) = eosData(EOS_ENTR)*ergPerKg_to_kbPerBaryon

end subroutine eos_hybridHelmholtz
