!!****f* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_init
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
!!  rt_init
!!
!! SYNOPSIS
!!
!!  call rt_init ()
!!
!! DESCRIPTION
!!
!!  Initialize Thornado
!!
!! ARGUMENTS
!!
!!***

subroutine rt_init()

#include "Simulation.h"
#include "constants.h"

  use rt_data
  use RadTrans_data, ONLY : rt_gcMask, rt_meshMe
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use ThornadoInitializationModule, ONLY : InitThornado
#ifdef FLASH_EOS_WEAKLIB
  use eos_wlData, ONLY : eos_pointer
#endif

  implicit none

  character(len=MAX_STRING_LENGTH) :: eos_file
  real :: eos_gamma
  logical :: Verbose

  call RuntimeParameters_get ("rt_writeTimers", rt_writeTimers)

  call RuntimeParameters_get ("rt_doExplicit", rt_doExplicit)
  call RuntimeParameters_get ("rt_doImplicit", rt_doImplicit)

  call RuntimeParameters_get ("rt_eL",    rt_eL)
  call RuntimeParameters_get ("rt_eR",    rt_eR)
  call RuntimeParameters_get ("rt_zoomE", rt_zoomE)
  call RuntimeParameters_get ("rt_bcE",   rt_bcE)

  call RuntimeParameters_get ("rt_use_emab", rt_use_emab)
  call RuntimeParameters_get ("rt_use_iso",  rt_use_iso)
  call RuntimeParameters_get ("rt_use_nes",  rt_use_nes)
  call RuntimeParameters_get ("rt_use_pair", rt_use_pair)

  call RuntimeParameters_get ("rt_emab_file", rt_emab_file)
  if ( .not. rt_use_emab ) rt_emab_file = ""
  call RuntimeParameters_get ("rt_iso_file", rt_iso_file)
  if ( .not. rt_use_iso ) rt_iso_file = ""
  call RuntimeParameters_get ("rt_nes_file", rt_nes_file)
  if ( .not. rt_use_nes ) rt_nes_file = ""
  call RuntimeParameters_get ("rt_pair_file", rt_pair_file)
  if ( .not. rt_use_pair ) rt_pair_file = ""

  call RuntimeParameters_get ("rt_positivityLimiter", rt_positivityLimiter)
  call RuntimeParameters_get ("rt_UpperBry1", rt_UpperBry1)
  rt_UpperBry1 = NEAREST(rt_UpperBry1,-1.0)

  rt_gcMask(THORNADO_BEGIN:THORNADO_END) = .TRUE.

  Verbose = ( rt_meshMe == MASTER_PE )
#ifdef FLASH_EOS_WEAKLIB
  call RuntimeParameters_get("eos_file", eos_file)
  call InitThornado( THORNADO_NNODES, NDIM, THORNADO_NE, &
     THORNADO_SWE, rt_eL, rt_eR, rt_zoomE, rt_bcE, &
     EquationOfStateTableName_Option = eos_file, &
     External_EOS = eos_pointer, &
     PositivityLimiter_Option = rt_positivityLimiter, &
     UpperBry1_Option = rt_UpperBry1, &
     OpacityTableName_EmAb_Option = rt_emab_file, &
     OpacityTableName_Iso_Option = rt_iso_file, &
     OpacityTableName_NES_Option = rt_nes_file, &
     OpacityTableName_Pair_Option = rt_pair_file, &
     Verbose_Option = Verbose )
#else
  call RuntimeParameters_get("gamma", eos_gamma)
  call InitThornado( THORNADO_NNODES, NDIM, THORNADO_NE, &
     THORNADO_SWE, rt_eL, rt_eR, rt_zoomE, rt_bcE, &
     Gamma_IDEAL_Option = eos_gamma, &
     PositivityLimiter_Option = rt_positivityLimiter, &
     UpperBry1_Option = rt_UpperBry1, &
     Verbose_Option = Verbose )
#endif

  return
end subroutine rt_init
