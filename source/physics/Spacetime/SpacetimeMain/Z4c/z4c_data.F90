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
!! @brief Main data module for Z4c solver

!> @ingroup Z4c
!! Variables storing data required in the Z4c solver
module z4c_data

#include "Simulation.h"

   implicit none

   public

   real, save :: z4c_eta
   real, save :: z4c_kappa1
   real, save :: z4c_kappa2
   real, save :: z4c_KOSigma

   real, save :: z4c_cfl

   ! RHS indices from MoL - capitalized for consistency with *_VAR indices
   integer, save :: CHI_RHS, &
                    GAMTILDE_LL_00_RHS, &
                    GAMTILDE_LL_01_RHS, &
                    GAMTILDE_LL_02_RHS, &
                    GAMTILDE_LL_11_RHS, &
                    GAMTILDE_LL_12_RHS, &
                    GAMTILDE_LL_22_RHS, &
                    KHAT_RHS, &
                    ATILDE_LL_00_RHS, &
                    ATILDE_LL_01_RHS, &
                    ATILDE_LL_02_RHS, &
                    ATILDE_LL_11_RHS, &
                    ATILDE_LL_12_RHS, &
                    ATILDE_LL_22_RHS, &
                    THETAFUNC_RHS, &
                    GAMTILDE_U_0_RHS, &
                    GAMTILDE_U_1_RHS, &
                    GAMTILDE_U_2_RHS, &
                    ALPHA_RHS, &
                    BETA_U_0_RHS, &
                    BETA_U_1_RHS, &
                    BETA_U_2_RHS
end module z4c_data
