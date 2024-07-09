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
!! @brief This routine does ideal gamma law specific initialization
!!
!!
!!

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"
subroutine eos_idealGammaInit()

  use Eos_data
  use eos_idealGammaData, ONLY : eos_gammam1
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abort
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
       Driver_getComm
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

  logical :: threadWithinBlockBuild
  
  ! Everybody should know this
  call Driver_getMype(MESH_COMM,eos_meshMe)
  call Driver_getNumProcs(MESH_COMM,eos_meshNumProcs)
  call Driver_getComm(MESH_COMM,eos_meshComm)

  call PhysicalConstants_get("ideal gas constant", eos_gasConstant)

  call RuntimeParameters_get("gamma", eos_gamma)
  call RuntimeParameters_get("eos_singleSpeciesA", eos_singleSpeciesA)
  call RuntimeParameters_get("eos_singleSpeciesZ", eos_singleSpeciesZ)
  call RuntimeParameters_get("eos_logLevel", eos_logLevel)
  call RuntimeParameters_get("smalle",eos_smalle)
  call RuntimeParameters_get("eintSwitch",eos_eintSwitch)
#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abort("[eos_idealGammaInit] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif


  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadEosWithinBlock", eos_threadWithinBlock)

  if (eos_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[eos_idealGammaInit]')
     eos_threadWithinBlock = .false.
  end if


  eos_type=EOS_GAM

  eos_gammam1 = 1.0/(eos_gamma-1.0)

  return
end subroutine eos_idealGammaInit
