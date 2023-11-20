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
  use Grid_data, ONLY : gr_minRefine, gr_lrefineMax
  use MoL_interface, ONLY : MoL_registerVariable
  use RadTrans_data, ONLY : rt_gcMask, rt_meshMe, rt_meshNumProcs, rt_str_geometry
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use ThornadoInitializationModule, ONLY : InitThornado
#ifdef FLASH_EOS_WEAKLIB
  use eos_wlData, ONLY : eos_pointer
#endif

#ifdef FLASH_GRID_PARAMESH
  use physicaldata, ONLY : interp_mask_unk, interp_mask_unk_res
#endif

  implicit none

  character(len=MAX_STRING_LENGTH) :: eos_file
  real :: eos_gamma
  logical :: Verbose
  integer :: iS, iCR, iE, iNodeE, ivar
  character(len=4) :: unk_name

  call RuntimeParameters_get ("rt_writeTimers", rt_writeTimers)

  call RuntimeParameters_get ("rt_doExplicit", rt_doExplicit)
  call RuntimeParameters_get ("rt_doImplicit", rt_doImplicit)

  call RuntimeParameters_get ("rt_eL",    rt_eL)
  call RuntimeParameters_get ("rt_eR",    rt_eR)
  call RuntimeParameters_get ("rt_zoomE", rt_zoomE)
  call RuntimeParameters_get ("rt_bcE",   rt_bcE)

  call RuntimeParameters_get ("rt_units", rt_units)

  call RuntimeParameters_get ("rt_use_emab", rt_use_emab)
  call RuntimeParameters_get ("rt_use_iso",  rt_use_iso)
  call RuntimeParameters_get ("rt_use_nes",  rt_use_nes)
  call RuntimeParameters_get ("rt_use_pair", rt_use_pair)
  call RuntimeParameters_get ("rt_use_brem", rt_use_brem)

  call RuntimeParameters_get ("rt_emab_file", rt_emab_file)
  if ( .not. rt_use_emab ) rt_emab_file = ""
  call RuntimeParameters_get ("rt_iso_file", rt_iso_file)
  if ( .not. rt_use_iso ) rt_iso_file = ""
  call RuntimeParameters_get ("rt_nes_file", rt_nes_file)
  if ( .not. rt_use_nes ) rt_nes_file = ""
  call RuntimeParameters_get ("rt_pair_file", rt_pair_file)
  if ( .not. rt_use_pair ) rt_pair_file = ""
  call RuntimeParameters_get ("rt_brem_file", rt_brem_file)
  if ( .not. rt_use_brem ) rt_brem_file = ""

  call RuntimeParameters_get ("rt_positivityLimiter", rt_positivityLimiter)
  call RuntimeParameters_get ("rt_UpperBry1", rt_UpperBry1)
  call RuntimeParameters_get ("rt_troubledCellIndicator", rt_troubledCellIndicator)
  call RuntimeParameters_get ("rt_cTCI", rt_cTCI)
  call RuntimeParameters_get ("rt_slopeLimiter", rt_slopeLimiter)
  call RuntimeParameters_get ("rt_energyLimiter", rt_energyLimiter)
  rt_UpperBry1 = NEAREST(rt_UpperBry1,-1.0)

  call RuntimeParameters_get ("rt_M_outer", rt_M_outer)
  call RuntimeParameters_get ("rt_M_inner", rt_M_inner)
  call RuntimeParameters_get ("rt_MaxIter_outer", rt_MaxIter_outer)
  call RuntimeParameters_get ("rt_MaxIter_inner", rt_MaxIter_inner)
  call RuntimeParameters_get ("rt_Rtol_outer", rt_Rtol_outer)
  call RuntimeParameters_get ("rt_Rtol_inner", rt_Rtol_inner)
  call RuntimeParameters_get ("rt_Include_LinCorr", rt_Include_LinCorr)

  rt_gcMask(THORNADO_BEGIN:THORNADO_END) = .TRUE.

  rt_offGridFluxR = 0.0

  ! interpolation of DG variables only works for "native" interpolation mode
#if defined(FLASH_GRID_PARAMESH)
  ! use DG interpolation/averaging for prolongation/restriction
  interp_mask_unk    (THORNADO_BEGIN:THORNADO_END) = 40
  ! if curvilinear, these looks like they are overwritten by mpi_amr_1blk_restrict.F90 on lines 332
  interp_mask_unk_res(THORNADO_BEGIN:THORNADO_END) = 40
#else
  if ( gr_lrefineMax > gr_minRefine ) then
     call Driver_abort("lrefine_max > lrefine_min is only supported for PM4 w/ -gridinterpolation=native)")
   end if
#endif

  rt_wMatrRHS(1:2) = 1.0
  rt_wMatrRHS(3:5) = 0.0

  Verbose = ( rt_meshMe == MASTER_PE )
#ifdef FLASH_EOS_WEAKLIB
  call RuntimeParameters_get("eos_file", eos_file)
  call InitThornado( THORNADO_NNODES, NDIM, THORNADO_NE, &
     THORNADO_SWE, rt_eL, rt_eR, rt_zoomE, rt_bcE, THORNADO_NSPECIES, &
     EquationOfStateTableName_Option = eos_file, &
     External_EOS = eos_pointer, &
     PositivityLimiter_Option = rt_positivityLimiter, &
     UpperBry1_Option = rt_UpperBry1, &
     TroubledCellIndicator_Option = rt_troubledCellIndicator, &
     C_TCI_Option = rt_cTCI, &
     SlopeLimiter_Option = rt_slopeLimiter, &
     EnergyLimiter_Option = rt_energyLimiter, &
     OpacityTableName_EmAb_Option = rt_emab_file, &
     OpacityTableName_Iso_Option = rt_iso_file, &
     OpacityTableName_NES_Option = rt_nes_file, &
     OpacityTableName_Pair_Option = rt_pair_file, &
     OpacityTableName_Brem_Option = rt_brem_file, &
     M_outer_Option = rt_M_outer, &
     M_inner_Option = rt_M_inner, &
     MaxIter_outer_Option = rt_MaxIter_outer, &
     MaxIter_inner_Option = rt_MaxIter_inner, &
     Rtol_outer_Option = rt_Rtol_outer, &
     Rtol_inner_Option = rt_Rtol_inner, &
     Include_LinCorr_Option = rt_Include_LinCorr, &
     wMatrRHS_Option = rt_wMatrRHS, &
     Include_NES_Option = rt_use_nes, &
     Include_Pair_Option = rt_use_pair, &
     Include_Brem_Option = rt_use_brem, &
     ActivateUnits_Option = rt_units, &
     CoordinateSystem_Option = rt_str_geometry, &
     Verbose_Option = Verbose )
#else
  call RuntimeParameters_get("gamma", eos_gamma)
  call InitThornado( THORNADO_NNODES, NDIM, THORNADO_NE, &
     THORNADO_SWE, rt_eL, rt_eR, rt_zoomE, rt_bcE, THORNADO_NSPECIES, &
     EquationOfStateTableName_Option = '', &
     External_EOS = 0, &
     Gamma_IDEAL_Option = eos_gamma, &
     PositivityLimiter_Option = rt_positivityLimiter, &
     UpperBry1_Option = rt_UpperBry1, &
     SlopeLimiter_Option = rt_slopeLimiter, &
     EnergyLimiter_Option = rt_energyLimiter, &
     CoordinateSystem_Option = rt_str_geometry, &
     ActivateUnits_Option = rt_units, &
     Verbose_Option = Verbose )
#endif

  do iS = 1, THORNADO_NSPECIES
     do iCR = 1, THORNADO_NMOMENTS
        do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE
           do iNodeE = 1, THORNADO_NNODESE

              ivar = THORNADO_BEGIN &
                 + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                 + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                 + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE) &
                 + iNodeE - 1

              rt_ivar(iNodeE,iE,iCR,iS) = ivar

              call Simulation_mapIntToStr(ivar, unk_name, MAPBLOCK_UNK)
              call MoL_registerVariable(unk_name, ivar, rt_irhs(iNodeE,iE,iCR,iS))

           end do
        end do
     end do
  end do

  call RuntimeParameters_get( "rt_D_0"  , rt_D_0 )
  call RuntimeParameters_get( "rt_Chi"  , rt_Chi )
  call RuntimeParameters_get( "rt_Sigma", rt_Sigma )

  return
end subroutine rt_init
