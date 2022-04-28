!!****if* source/physics/Eos/EosMain/Eos_init
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
!!  Eos_init
!!
!! 
!! SYNOPSIS
!!
!!  call Eos_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants facilities. This is the version for simple Gamma law
!!  implementation with a single fluid.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos. 
!!   Particular implementations (Gamma,Helmholz,etc) of the unit
!!   define their own runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might overwrite these values with the 
!!   flash.par values for your specific run.  
!!
!!   gamma[Real]   -- Ratio of specific heats, default 1.6667
!!   eos_singleSpeciesA[Real]  -- Nucleon number for the gas, default 1.00
!!   eos_singleSpeciesZ[Real]  -- Proton number for the gas, default 1.00
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based variables GAMC_VAR and GAME_VAR in Simulation.h
!!
!!***
subroutine Eos_init()

  use Eos_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abort
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs, &
       Driver_getComm
  use eos_localInterface, ONLY : eos_initMgamma, eos_initHelmholtz, eos_initGamma
  use Logfile_interface, ONLY : Logfile_stamp

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"  

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
     call Driver_abort("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif


  call eos_fillMapLookup()

  call eos_initGamma()
  call eos_initMgamma()
  call eos_initHelmholtz()
  call eos_initWeaklib()


  call RuntimeParameters_get("threadWithinBlockBuild", threadWithinBlockBuild)
  call RuntimeParameters_get("threadEosWithinBlock", eos_threadWithinBlock)

  if (eos_threadWithinBlock .and. .not. threadWithinBlockBuild) then
     call Logfile_stamp('WARNING! Turning off within block threading '//&
          'because FLASH is not built appropriately','[Eos_init]')
     eos_threadWithinBlock = .false.
  end if

  return

end subroutine Eos_init
