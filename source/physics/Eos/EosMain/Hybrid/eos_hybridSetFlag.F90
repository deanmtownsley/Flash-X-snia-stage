!!****if* source/physics/Eos/EosMain/Hybrid/eos_hybridSetFlag
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

subroutine eos_hybridSetFlag(mode, vecLen, eosData, massFrac, mask, vecB, vecE, eos_hybFlag)

   implicit none

#include "constants.h"
#include "Eos.h"
#include "Simulation.h"

   integer, INTENT(in) :: mode, vecLen
   real, INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData
   logical, optional, INTENT(in), target, dimension(EOS_VARS + 1:EOS_NUM) :: mask
   real, optional, INTENT(in), dimension(NSPECIES*vecLen) :: massFrac
   integer, optional, INTENT(in) :: vecB, vecE
   integer, optional, INTENT(out) :: eos_hybFlag(vecLen)

   if (present(eos_hybFlag)) eos_hybFlag = 0

   return
end subroutine eos_hybridSetFlag
