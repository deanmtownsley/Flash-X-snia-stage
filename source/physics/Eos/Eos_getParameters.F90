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
!! @brief Eos_getParameters stub

!> @ingroup physics_Eos
!!
!! @brief return information of interest to other units
!!
!! @details
!! @anchor Eos_getParameters_stub
!!  This interface assumes that Eos initialization has already taken place when
!!  Eos_getParameters is called.
!!
!!  Since this interface uses optional arguments, all routines calling this routine must include
!!  a use Eos_interface statement, possibly with "ONLY" attribute listing Eos_getParameters
!!  explicitly, e.g.
!!      use Eos_interface, ONLY:  Eos_getParameters
!!
!!  
!!  @param eintSwitch -          value of the Eos unit's eintSwitch runtime parameter.
!!  @paraminputsAreUnchanged -  Indicates whether calls to Eos (or Eos_multiDim, etc.) can result
!!                         in modification of some state variables that should be, strictly
!!                         speaking, input only to the EOS in the given mode.
!!                         If this is true, then calls to Eos with MODE_DENS_PRES can modify
!!                         the pressure, and calls with MODE_DENS_EI can modify enery variables.
!!   @paraminputTempIsGuess -    Indicates whether the Eos implementation uses the temperature
!!                         provided on entry to a call as in initial gues in an iterative scheme.
!!   @param constantGammaC -      Indicates whether the gamc returned by Eos will always be constant.
!!   @param inputMassFracNeeded - Indicates whether the Eos implementation makes use of mass fractions.
!!   @param smalle -              value of the Eos unit's smallE runtime parameter.
!!   @param smallE1,smallE2,smallE3 - values of lower bounds for temperature components in multiTemp 
!!                         implementations.
!!
subroutine Eos_getParameters(eintSwitch,inputsAreUnchanged,inputTempIsGuess,constantGammaC,&
     inputMassFracNeeded,smalle,smallE1,smallE2,smallE3)

  implicit none

  real,OPTIONAL,intent(OUT) :: eintSwitch
  logical,OPTIONAL,intent(OUT) :: inputsAreUnchanged
  logical,OPTIONAL,intent(OUT) :: inputTempIsGuess
  logical,OPTIONAL,intent(OUT) :: constantGammaC
  logical,OPTIONAL,intent(OUT) :: inputMassFracNeeded
  real,OPTIONAL,intent(OUT) :: smalle
  real,OPTIONAL,intent(OUT) :: smallE1,smallE2,smallE3

  if (present(eintSwitch)) eintSwitch = 0.0
  if (present(inputsAreUnchanged)) inputsAreUnchanged = .FALSE.
  if (present(inputTempIsGuess)) inputTempIsGuess = .FALSE.
  if (present(constantGammaC)) constantGammaC = .FALSE.
  if (present(inputMassFracNeeded)) inputMassFracNeeded = .FALSE.
  if (present(smalle)) smalle = 0.0
  if (present(smallE1)) smallE1 = 0.0
  if (present(smallE2)) smallE2 = 0.0
  if (present(smallE3)) smallE3 = 0.0

  return
end subroutine Eos_getParameters
