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
!> @ingroup physics_eos_EosMain_WeakLib
!!
!! @brief Obtain the energy shift assumed by WeakLib
!!
!! @details
!!
!! This procedure can be used to obtain the neutron and proton chemical
!! potentials, and optionally the electron and electron neutrino chemical
!! potentials.  This procedure assumes that the provided density, temperature
!! and electron fraction represent a thermodynamically consistent state.
!!
!! @param xDens   Density (g/cc)
!! @param xTemp   Temperature (K)
!! @param xYe     Electron fraction
!! @param xMu_n   Neutron chemical potential (erg)
!! @param xMu_p   Proton chemical potential (erg)
!! @param xMu_e   Electron chemical potential (erg)

subroutine eos_wlPotentials(xDens, xTemp, xYe, xMu_n, xMu_p, xMu_e)
   use eos_wlData, only: EosNewTable

   use wlInterpolationModule, only: LogInterpolateSingleVariable_3D_Custom_Point

   implicit none

   real, intent(in) :: xDens, xTemp, xYe
   real, intent(out) :: xMu_n, xMu_p
   real, intent(out), optional :: xMu_e

   integer :: iD_T, iT_T, iY_T
   integer :: iMn_T, iMp_T, iMe_T

   real :: OS_Me, OS_Mp, OS_Mn

   iD_T = EosNewTable%TS%Indices%iRho
   iT_T = EosNewTable%TS%Indices%iT
   iY_T = EosNewTable%TS%Indices%iYe

   iMn_T = EosNewTable%DV%Indices%iNeutronChemicalPotential
   iMp_T = EosNewTable%DV%Indices%iProtonChemicalPotential
   iMe_T = EosNewTable%DV%Indices%iElectronChemicalPotential

   OS_Mn = EosNewTable%DV%Offsets(iMn_T)
   OS_Mp = EosNewTable%DV%Offsets(iMp_T)
   OS_Me = EosNewTable%DV%Offsets(iMe_T)

   associate (D_T => EosNewTable%TS%States(iD_T)%Values, &
              T_T => EosNewTable%TS%States(iT_T)%Values, &
              Y_T => EosNewTable%TS%States(iY_T)%Values, &
              Me_T => EosNewTable%DV%Variables(iMe_T)%Values, &
              Mn_T => EosNewTable%DV%Variables(iMn_T)%Values, &
              Mp_T => EosNewTable%DV%Variables(iMp_T)%Values)

      call LogInterpolateSingleVariable_3D_Custom_Point(xDens, xTemp, xYe, &
                                                        D_T, T_T, Y_T, &
                                                        OS_Mn, Mn_T, xMu_n)
      call LogInterpolateSingleVariable_3D_Custom_Point(xDens, xTemp, xYe, &
                                                        D_T, T_T, Y_T, &
                                                        OS_Mp, Mp_T, xMu_p)

      if (present(xMu_e)) &
         call LogInterpolateSingleVariable_3D_Custom_Point(xDens, xTemp, xYe, &
                                                           D_T, T_T, Y_T, &
                                                           OS_Me, Me_T, xMu_e)
   end associate

end subroutine eos_wlPotentials
