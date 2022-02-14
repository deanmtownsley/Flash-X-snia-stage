!!****f* source/physics/Eos/Eos
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
!!  Eos
!!
!! SYNOPSIS
!!
!!       call Eos(integer(IN) :: mode,
!!                integer(IN) :: vecLen,
!!                real(INOUT) :: eosData(vecLen*EOS_NUM),
!!      optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional, logical(IN),target :: mask(EOS_VARS+1:EOS_NUM),
!!      optional, integer(IN) :: vecBegin,
!!      optional, integer(IN) :: vecEnd    )
!!
!! DESCRIPTION
!!
!!  This routine applies the equation of state to thermodynamic 
!!  quantities at one or more grid cells.  The number of cells is 
!!  determined by the argument veclen.  The routine expects data packaged
!!  for it in the 1d array, eosData.  The data in eosData is organized as:
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable, and so on. The number and order of
!!  variables in the array is determined by the constants defined in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and internal energy are generated
!!  as output; if mode=MODE_DENS_EI, density and internal energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  internal energy and temperature are generated as output. Note that
!!  internal energy is EINT_VAR, not ENER_VAR.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!
!!
!!  ARGUMENTS
!!
!!  mode :    Selects the mode of operation of the Eos unit.
!!             The valid values are MODE_DENS_EI, MODE_DENS_PRES and  
!!             MODE_DENS_TEMP as decribed above.
!!
!!  vecLen   : number of points (cells) for which the eosData array is sized.
!!             If vecBegin and vecEnd are not present, this is also the
!!             number of points (cells) for which EOS computation is to be done.
!!
!!  eosData  : This array is the data structure through which variable values are 
!!             passed in and out of the Eos routine. The arrays is sized as 
!!             EOS_NUM*vecLen. EOS_NUM, and individual input and output
!!             Eos variables are defined in Eos.h. The array is organizes such that
!!             the first 1:vecLen entries represent the first Eos variable, vecLen+1:
!!             2*vecLen represent the second Eos variable and so on. 
!!
!!  massFrac : Contains the mass fractions of the species included in
!!             the simulation. The array is sized as NSPECIES*vecLen.
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
!!             An implementation that does not need derivative quantities should
!!             set the mask equal to .false.
!!
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
!!
!!
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_wrapped  sets up the data structure.
!!
!!***

subroutine Eos(mode, vecLen, eosData, massFrac, mask, vecBegin,vecEnd, diagFlag)

!==============================================================================

  implicit none
#include "Eos.h"
#include "Simulation.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer, optional, INTENT(out)    :: diagFlag

  if (present(diagFlag)) diagFlag = 8

  return
end subroutine Eos
