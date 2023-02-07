!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief Spacetime_init implementation

!> @ingroup SpacetimeMain
!!
!! @brief Implements Spacetime_init
!!
!! @stubref{Spacetime_init}
subroutine Spacetime_init()
   ! Intentionally leaving off the use statement since we need all variables in the module
   use z4c_data

   use RuntimeParameters_interface, only: RuntimeParameters_get
   use MoL_interface, only: MoL_registerVariable

#include "Z4c.h"
#include "constants.h"

   implicit none

   call RuntimeParameters_get("z4c_cfl", z4c_cfl)

   call RuntimeParameters_get("z4c_eta", z4c_eta)
   call RuntimeParameters_get("z4c_kappa1", z4c_kappa1)
   call RuntimeParameters_get("z4c_kappa2", z4c_kappa2)
   call RuntimeParameters_get("z4c_KOSigma", z4c_KOSigma)

   call RuntimeParameters_get("z4c_harmonicLapse", z4c_harmonicLapse)
   call RuntimeParameters_get("z4c_stationaryShift", z4c_stationaryShift)

   call MoL_registerVariable("Z400", CHI_VAR, CHI_RHS)
   call MoL_registerVariable("Z401", GAMTILDE_LL_00_VAR, GAMTILDE_LL_00_RHS)
   call MoL_registerVariable("Z402", GAMTILDE_LL_01_VAR, GAMTILDE_LL_01_RHS)
   call MoL_registerVariable("Z403", GAMTILDE_LL_02_VAR, GAMTILDE_LL_02_RHS)
   call MoL_registerVariable("Z404", GAMTILDE_LL_11_VAR, GAMTILDE_LL_11_RHS)
   call MoL_registerVariable("Z405", GAMTILDE_LL_12_VAR, GAMTILDE_LL_12_RHS)
   call MoL_registerVariable("Z406", GAMTILDE_LL_22_VAR, GAMTILDE_LL_22_RHS)
   call MoL_registerVariable("Z407", KHAT_VAR, KHAT_RHS)
   call MoL_registerVariable("Z408", ATILDE_LL_00_VAR, ATILDE_LL_00_RHS)
   call MoL_registerVariable("Z409", ATILDE_LL_01_VAR, ATILDE_LL_01_RHS)
   call MoL_registerVariable("Z410", ATILDE_LL_02_VAR, ATILDE_LL_02_RHS)
   call MoL_registerVariable("Z411", ATILDE_LL_11_VAR, ATILDE_LL_11_RHS)
   call MoL_registerVariable("Z412", ATILDE_LL_12_VAR, ATILDE_LL_12_RHS)
   call MoL_registerVariable("Z413", ATILDE_LL_22_VAR, ATILDE_LL_22_RHS)
   call MoL_registerVariable("Z414", THETAFUNC_VAR, THETAFUNC_RHS)
   call MoL_registerVariable("Z415", GAMTILDE_U_0_VAR, GAMTILDE_U_0_RHS)
   call MoL_registerVariable("Z416", GAMTILDE_U_1_VAR, GAMTILDE_U_1_RHS)
   call MoL_registerVariable("Z417", GAMTILDE_U_2_VAR, GAMTILDE_U_2_RHS)
   call MoL_registerVariable("Z418", ALPHA_VAR, ALPHA_RHS)

   if (.not. z4c_stationaryShift) then
      call MoL_registerVariable("Z419", BETA_U_0_VAR, BETA_U_0_RHS)
      call MoL_registerVariable("Z420", BETA_U_1_VAR, BETA_U_1_RHS)
      call MoL_registerVariable("Z421", BETA_U_2_VAR, BETA_U_2_RHS)
   end if
end subroutine Spacetime_init
