!!****if* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_data
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
!!  NAME 
!!
!!  rt_data
!!
!!  SYNOPSIS
!!   use rt_data
!!
!!  DESCRIPTION 
!!    Stores local data for Thornado
!!
!!***

#include "Simulation.h"

module rt_data

  implicit none
  
  logical, save :: rt_writeTimers

  logical, save :: rt_doExplicit, rt_doImplicit
  integer, save :: rt_bcE
  real,    save :: rt_eL, rt_eR, rt_zoomE

  logical, save :: rt_units

  logical, save :: rt_positivityLimiter
  real,    save :: rt_UpperBry1

  logical, save :: rt_troubledCellIndicator
  real,    save :: rt_cTCI

  logical, save :: rt_slopeLimiter
  logical, save :: rt_energyLimiter

  logical, save :: rt_use_emab, rt_use_iso, rt_use_nes, rt_use_pair, rt_use_brem
  character(len=80), save :: rt_emab_file, rt_iso_file, rt_nes_file, rt_pair_file, rt_brem_file

  real,    save :: rt_Op_MinD, rt_Op_MaxD

  logical, save :: rt_muShift

  real,    save :: rt_offGridFluxR(2*THORNADO_NMOMENTS)

  ! neutrino-matter solver parameters
  integer, save :: rt_M_outer, rt_M_inner
  integer, save :: rt_MaxIter_outer, rt_MaxIter_inner
  real,    save :: rt_Rtol_outer, rt_Rtol_inner
  logical, save :: rt_Include_LinCorr

  real,    save :: rt_wMatrRHS(5)
  logical, save :: rt_freezeOpacities

  integer, save :: rt_ivar(1:THORNADO_NNODESE, &
                           1-THORNADO_SWE:THORNADO_NE+THORNADO_SWE, &
                           1:THORNADO_NMOMENTS, &
                           1:THORNADO_NSPECIES)

  integer, save :: rt_irhs(1:THORNADO_NNODESE, &
                           1-THORNADO_SWE:THORNADO_NE+THORNADO_SWE, &
                           1:THORNADO_NMOMENTS, &
                           1:THORNADO_NSPECIES)

  real,    save :: rt_D_0, rt_Chi, rt_Sigma

  real, save :: rt_transitionDens
  
end module rt_data
