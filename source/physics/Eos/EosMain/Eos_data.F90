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
!! @brief This is the data module that contains data that is common to all the
!!      implementations of Eos 
!! 
!!


module Eos_data
#include "Simulation.h"
#include "Eos.h"
#include "Eos_map.h"

  real, save :: eos_gasConstant 
  real, save :: eos_smalle
  real, save :: eos_gamma
  real, save :: eos_singleSpeciesA
  real, save :: eos_singleSpeciesZ
!  real, save :: eos_gammam1
  real, save :: eos_eintSwitch
  integer, save :: eos_type, eos_meshMe, eos_meshNumProcs, eos_globalMe
  integer, save :: eos_meshComm

  ! maximum number of iterations for the Newton loop to find T from e
  integer, save :: eos_maxNewton

  ! how accurately to do the Newton iteration to get T from e
  real, save :: eos_tol

  ! force the iterative solver to leave inputs alone (always true in MODE_DENS_TEMP)
  logical, save :: eos_forceConstantInput

  real, save :: eos_smallt
  real, save :: eos_largeT = 1.0e10 ! used by some modified implementations of Newton-Raphson
  real, save :: eos_smallRho

  integer, save :: eos_logLevel = 700

  ! Some stuff that is only used by multiTemp implementations
  integer, save :: eos_combinedTempRule = -1

  integer, save :: eos_entrEleScaleChoice = -1

  real, save :: eos_smallEion=0.0, eos_smallEele=0.0, eos_smallErad=0.0

#ifdef FIXEDBLOCKSIZE
  real,save,dimension(NSPECIES*MAXCELLS) :: eos_massFr
  real,save,dimension(EOS_NUM*MAXCELLS) :: eos_inOut
#else
  real,save, allocatable :: eos_inOut(:),eos_massFr(:)
#endif

  real,parameter :: eos_pradScaleFactor = 1.0 !hardwired for now
#ifdef FLLM_VAR
  integer,parameter :: eos_pradScaleVar = FLLM_VAR
#else
  integer,parameter :: eos_pradScaleVar = -1
#endif

  integer, save, dimension(1:EOSMAP_NUM_ROLES, 1:2, 1:5) :: eos_mapLookup
  logical, save :: eos_threadWithinBlock = .false.
  logical, save :: eos_doYe = .false.
end module Eos_data
