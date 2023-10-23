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
!! @brief Eos_wlPotentials stub

!> @ingroup physics_Eos_EosMain_WeakLib
!!
!! @brief Obtain the energy shift assumed by WeakLib
!!
!! @details
!! @anchor Eos_wlPotentials_stub
!!
!! This procedure can be used to obtain the neutron and proton chemical
!! potentials, and optionally the electron and electron neutrino chemical
!! potentials.  This procedure assumes that the provided density, temperature
!! and electron fraction represent a thermodynamically consistent state.
!!
!! @param xDens   Density (g/cc)
!! @param xTemp   Temperature (K)
!! @param xYe     Electron fraction
!! @param xMu_n   Neutron chemical potential
!! @param xMu_p   Proton chemical potential
!! @param xMu_e   Electron chemical potential
subroutine Eos_wlPotentials(xDens, xTemp, xYe, xMu_n, xMu_p, xMu_e)
   use eos_wlData, only: EosNewTable

   use wlInterpolationModule, only: LogInterpolateSingleVariable_3D_Custom_Point

   implicit none

   real, intent(in) :: xDens, xTemp, xYe
   real, intent(out) :: xMu_n, xMu_p
   real, intent(out), optional :: xMu_e

   integer :: iD_T, iT_T, iY_T
   integer :: iMe_T, iMp_T, iMn_T

   real :: OS_Me, OS_Mp, OS_Mn

   integer :: error_flag

   iD_T = EosNewTable%TS%Indices%iRho
   iT_T = EosNewTable%TS%Indices%iT
   iY_T = EosNewTable%TS%Indices%iYe

   iMe_T = EosNewTable%DV%Indices%iElectronChemicalPotential
   iMp_T = EosNewTable%DV%Indices%iProtonChemicalPotential
   iMn_T = EosNewTable%DV%Indices%iNeutronChemicalPotential

end subroutine Eos_wlPotentials
