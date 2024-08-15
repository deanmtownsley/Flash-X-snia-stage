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
!! @brief Helmholtz specific implementation for getting parameters
!! 
!! @subref{Eos_getParameters}
!!



subroutine Eos_getParameters(eintSwitch,inputsAreUnchanged,inputTempIsGuess,constantGammaC,&
     inputMassFracNeeded,smalle,smallE1,smallE2,smallE3)

  use Eos_data, ONLY : eos_smalle, eos_eintSwitch
  use Eos_data, ONLY : eos_smallEion, eos_smallEele, eos_smallErad
  use eos_helmData, ONLY : eos_forceConstantInput, eos_useMultiSpecies
  use Driver_interface, ONLY : Driver_abort

  implicit none

#include "Simulation.h"

  real,OPTIONAL,intent(OUT) :: eintSwitch
  logical,OPTIONAL,intent(OUT) :: inputsAreUnchanged
  logical,OPTIONAL,intent(OUT) :: inputTempIsGuess
  logical,OPTIONAL,intent(OUT) :: constantGammaC
  logical,OPTIONAL,intent(OUT) :: inputMassFracNeeded
  real,OPTIONAL,intent(OUT) :: smalle
  real,OPTIONAL,intent(OUT) :: smallE1,smallE2,smallE3

  ! This assumes that runtime parameters have already been gotten.
  if (present(eintSwitch)) eintSwitch = eos_eintSwitch
  if (present(inputsAreUnchanged)) inputsAreUnchanged = eos_forceConstantInput
  if (present(inputTempIsGuess)) inputTempIsGuess = .TRUE.
  if (present(constantGammaC)) constantGammaC = .FALSE.

  if (present(inputMassFracNeeded)) inputMassFracNeeded = eos_useMultiSpecies
  if (present(smalle)) smalle = eos_smalle
  if (present(smallE1)) smallE1 = eos_smallEion
  if (present(smallE2)) smallE2 = eos_smallEele
  if (present(smallE3)) smallE3 = eos_smallErad

  return
end subroutine Eos_getParameters
