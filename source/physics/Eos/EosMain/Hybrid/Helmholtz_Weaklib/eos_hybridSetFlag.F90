!!****if* source/physics/Eos/EosMain/hybrid/Helmholtz_WeakLib/eos_hybridSetFlag
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
!!  eos_hybridSetFlag
!!
!! SYNOPSIS
!!
!!  call eos_hybridSetFlag(integer(IN)  :: mode,
!!                         integer(IN)  :: vecLen,
!!                         real(INOUT)  :: eosData(vecLen*EOS_NUM),
!!               optional, real(IN)     :: massFrac(vecLen*NSPECIES),
!!               optional, logical(IN)  :: mask(EOS_VARS+1:EOS_NUM),
!!                         integer(OUT) :: eos_hybFlag(vecLen))
!!
!! DESCRIPTION
!!
!!  This routine implements the hybrid (Helmholtz+WeakLib) equation of state.
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
!!  eos_hybFlag : flag to select the one or both of the EoS in the hybrid implementation
!!
!! NOTES
!!
!!  AH: several arguments are not needed, but are passed in almost all EoS routines so I have
!!      included them here.
!!
!!***

subroutine eos_hybridSetFlag(mode, vecLen, eosData, eos_hybFlag)

   use eos_hybridData, only: eos_hybTransitionDensLo, eos_hybTransitionDensHi

   implicit none

#include "constants.h"
#include "Eos.h"
#include "Simulation.h"

   integer, intent(in) :: mode, vecLen
   real, intent(in), dimension(vecLen,EOS_VARS) :: eosData
   integer, intent(out) :: eos_hybFlag(vecLen)

   integer :: vecB, vecE

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
