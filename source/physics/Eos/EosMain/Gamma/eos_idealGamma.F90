!!****if* source/physics/Eos/EosMain/Gamma/eos_idealGamma
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
!!
!! NAME
!!
!!  eos_idealGamma
!!
!! SYNOPSIS
!!
!!  call eos_idealGamma(integer(IN) :: mode,
!!                      integer(IN) :: vecLen,
!!                      real(INOUT) :: eosData(vecLen*EOS_NUM),
!!            optional, integer(IN) :: vecBegin,
!!            optional, integer(IN) :: vecEnd,
!!            optional, real(IN)    :: massFrac(vecLen*NSPECIES),
!!      optional,target,logical(IN) :: mask(EOS_VARS+1:EOS_NUM)  )
!!
!! DESCRIPTION
!!
!!
!!  This routine applies the gamma law equation of state to thermodynamic 
!!  quantities at one or more grid points.  The number of points is 
!!  determined by the argument veclen.  Data should be packaged for this 
!!  routine in the 1d array, eosData.  The data in eosData is organized as: 
!!  1:vecLen points contain the first variable, vecLen+1:2*vecLen points 
!!  contain the second variable and so on. The number and 
!!  order of variables in the array is determined by the constants defined
!!  in Eos.h.
!!  
!!  The routine takes different quantities as givens depending on the
!!  value of the mode variable: if mode=MODE_DENS_TEMP, density and
!!  temperature are taken as given, and pressure and energy are generated
!!  as output; if mode=MODE_DENS_EI, density and energy are taken as
!!  givens, and pressure and temperature are generated as output.  If
!!  mode=MODE_DENS_PRES, density and pressure are taken as givens, and
!!  energy and temperature are generated as output.
!!  
!!  In addition to pressure, temperature, and internal energy, which are
!!  always thermodynamically consistent after this call, other quantities
!!  such as the various thermodynamic partial derivatives can be
!!  calculated based on the values in the argument, mask.  mask is a
!!  logical array with one entry per quantity, with the order determined
!!  by constants defined in Eos.h (the same as those for the eosData
!!  argument); .true. means return the quantity, .false. means don't.
!!  
!!  This version applies to a single fluid.
!!  
!!  
!! ARGUMENTS 
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
!!  vecBegin : Index of first cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is 1.
!!  vecEnd   : Index of last cell in eosData to handle.
!!             Can be used to limit operation to a subrange of cells, untested.
!!             If not present, the default is vecLen.
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
!!
!! PARAMETERS
!!
!!     gamma        :  Ratio of specific heats for the simulated gas
!!     eos_singleSpeciesA       :  Mass of nucleus for the simulated gas
!!     eos_singleSpeciesZ       :  Proton number for the simulated gas
!!
!!
!! NOTES
!!
!!  NSPECIES is defined in Simulation.h.
!!
!!  EOS_VARS and EOS_NUM  are defined in Eos.h.
!!  Calling funtions should included Eos.h, in order to get the definitions of
!!  Eos-specific constants to be able to populate the eosData and mask arrays.
!!  
!!  MODE_DENS_TEMP, MODE_DENS_EI, and MODE_DENS_PRES are defined in constants.h.
!!
!!  User code should not call this implementation routine directly, but
!!  should call Eos and make sure that the desired Gamma implementation
!!  is included in the simulation configuration.
!!  All code calling the Eos interface should include a 
!!    use Eos_interface 
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use Eos_interface, ONLY:  Eos
!!  All routines calling this routine directly should include a 
!!    use eos_localInterface
!!  statement, preferable with "ONLY" attribute, e.g.,
!!    use eos_localInterface, ONLY:  Eos
!!
!!  For Gamma and Multigamma routines, the entropy and entropy derivatives 
!!  calculations have not been confirmed to be correct.  Use with caution.
!!
!! SEE ALSO
!! 
!!  Eos.h    defines the variables used.
!!  Eos_multiDim  sets up the data structure.
!!
!!***
!!NOVARIANTS

!#ifdef DEBUG_ALL
#define DEBUG_EOS
!#endif

subroutine eos_idealGamma(mode, pres, temp, dens, gamc, eint, entr, abar, zbar, massFrac, derivs)

!==============================================================================
  use Eos_data, ONLY : eos_gasConstant, eos_gamma, &
       eos_singleSpeciesA, eos_singleSpeciesZ
  use eos_idealGammaData, ONLY: eos_gammam1
  use Driver_interface, ONLY : Driver_abort


  implicit none

#include "constants.h"
#include "Eos.h"
#include "Simulation.h"

  !     Arguments
  integer, INTENT(in) :: mode
  real,INTENT(inout) :: pres, temp, dens, gamc, eint, entr, abar, zbar
  real, optional, INTENT(in),dimension(NSPECIES)    :: massFrac
  real, optional, INTENT(out),dimension(EOS_VARS+1:EOS_NUM) :: derivs

  real ::  ggprod, ggprodinv, gam1inv
  integer ::  dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, pel, ne, eta
  integer :: i, ilo,ihi

  ggprod = eos_gammam1 * eos_gasConstant

!============================================================================


  
!!NOTE  But for NSPECIES = 1, abar is NOT just 1.
!!NOTE The Flash2 routines initialize a fake "fluid"
!!NOTE with properties A=1, Z=1, and Eb=1.  Therefore flash2 can
!!NOTE  always call the equivalent of Multispecies_getSumInV and it
!!NOTE returns a default number.   In FLASH3, the decision has been made
!!NOTE that Eos/Gamma cannot be run with more than one fluid.

  gamc = eos_gamma
  abar = eos_singleSpeciesA
  zbar = eos_singleSpeciesZ    


  ! density, temperature taken as input
  if (mode == MODE_DENS_TEMP) then

     pres = eos_gasConstant*dens * &
                                 temp / abar
     eint = ggprod * temp &
                                 / abar
     entr = (pres/dens +  &
             &  eint)/temp


  ! density, internal energy taken as input
  elseif (mode == MODE_DENS_EI) then

     ggprodinv = 1. / ggprod
     gam1inv   = 1. / eos_gammam1
     pres = dens * &
                                    eint * gam1inv
     temp = eint * ggprodinv * &
                                    abar
     entr = (pres/dens +  &
             &  eint)/temp


  ! density, pressure taken as input
  elseif (mode == MODE_DENS_PRES) then

     ggprodinv = 1. / ggprod
     gam1inv   = 1. / eos_gammam1
     eint = pres * eos_gammam1 / &
                                   dens
     temp = eint * ggprodinv * &
                                   abar
     entr = (pres/dens +  &
             &  eint)/temp

  ! unrecognized value for mode
  else 
     call Driver_abort("[Eos] Unrecognized input mode given to Eos")
  endif


  if(present(derivs)) then
     derivs(EOS_DPT) = eos_gasConstant*dens / abar
     derivs(EOS_DPD) = eos_gasConstant*temp / abar
     derivs(EOS_DET) = ggprod / abar
     derivs(EOS_DED) = 0.
    ! Entropy derivatives   
     derivs(EOS_DST) = ( (derivs(EOS_DPT)  / dens + derivs(EOS_DET)) -&
          &                      (pres/ dens + eint)/ &
          &                      temp ) / temp
     derivs(EOS_DSD) = &
               ( ((derivs(EOS_DPD) - pres/dens) / &
        &          dens) + derivs(EOS_DED)) / temp


     derivs(EOS_PEL) = 0.
     derivs(EOS_NE) = 0.
     derivs(EOS_ETA) = 0.
     derivs(EOS_CV) = derivs(EOS_DET)
     derivs(EOS_CP) = eos_gamma*derivs(EOS_CV)
     
#ifdef EOS_CVELE
     derivs(EOS_CV) = derivs(EOS_DET) * &
          zbar / (zbar + 1)
#endif
  end if


  return
end subroutine eos_idealGamma

!! FOR FUTURE  : This section is not in use in FLASH 3 yet. none
!! of the current setups use entropy. This will be taken care of 
!! in future releases

!!..no matter what the input mode compute the entropy
!!..ignore the -chemical_potential*number_density part for now
!!$  dens_inv = 1.0e0/derivs(dens+ilo:+ihi)
!!$  temp_inv = 1.0e0/derivs(temp+ilo:+ihi)
!!$  stot     = (pres*dens_inv + derivs(eint+ilo:+ihi))*temp_inv 
!!$  dstotdd  = (derivs(EOS_DPD)*dens_inv - pres*dens_inv*dens_inv + derivs(EOS_DED))*temp_inv
!!$  dstotdt  = (derivs(EOS_DPT)*dens_inv + derivs(EOS_DET))*temp_inv  - (pres*dens_inv + derivs(eint+ilo:+ihi)) * temp_inv*temp_inv 
!!$  




