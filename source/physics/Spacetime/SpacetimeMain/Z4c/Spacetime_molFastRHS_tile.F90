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
!! @brief Spacetime_molFastRHS_tile implementation

!> @ingroup Z4c
!!
!! @brief Implements Spacetime_molFastRHS_tile for the Z4c solver
!!
!! @stubref{Spacetime_molFastRHS_tile}
subroutine Spacetime_molFastRHS_tile(tileDesc, t, activeRHS, dtWeight)
   use Grid_tile, only: Grid_tile_t

   use z4c_data, only: CHI_RHS, GAMTILDE_LL_00_RHS, GAMTILDE_LL_01_RHS, GAMTILDE_LL_02_RHS, &
                       GAMTILDE_LL_11_RHS, GAMTILDE_LL_12_RHS, GAMTILDE_LL_22_RHS, KHAT_RHS, &
                       ATILDE_LL_00_RHS, ATILDE_LL_01_RHS, ATILDE_LL_02_RHS, &
                       ATILDE_LL_11_RHS, ATILDE_LL_12_RHS, ATILDE_LL_22_RHS, &
                       THETAFUNC_RHS, GAMTILDE_U_0_RHS, GAMTILDE_U_1_RHS, GAMTILDE_U_2_RHS, &
                       ALPHA_RHS, BETA_U_0_RHS, BETA_U_1_RHS, BETA_U_2_RHS, &
                       eta => z4c_eta, kappa1 => z4c_kappa1, kappa2 => z4c_kappa2, KOSigma => z4c_KOSigma

   use MoL_interface, only: MoL_getDataPtr, MoL_releaseDataPtr

#include "Z4c.h"
#include "MoL.h"
#include "constants.h"

   implicit none

   type(Grid_tile_t), intent(in) :: tileDesc
   real, intent(in) :: t
   integer, intent(in) :: activeRHS
   real, intent(in) :: dtWeight

   real, dimension(:, :, :, :), pointer :: rhs, vars

   integer, dimension(LOW:HIGH, MDIM) :: lim, bcs
   integer :: i, j, k

   real :: del(MDIM), dx0, dx1, dx2, idx0, idx1, idx2

   real :: CHI
   real :: GAMTILDE_LL_00, GAMTILDE_LL_01, GAMTILDE_LL_02, GAMTILDE_LL_11, GAMTILDE_LL_12, GAMTILDE_LL_22
   real :: KHAT
   real :: ATILDE_LL_00, ATILDE_LL_01, ATILDE_LL_02, ATILDE_LL_11, ATILDE_LL_12, ATILDE_LL_22
   real :: THETAFUNC
   real :: GAMTILDE_U_0, GAMTILDE_U_1, GAMTILDE_U_2
   real :: ALPHA
   real :: BETA_U_0, BETA_U_1, BETA_U_2

   real :: dDCHI_L_0, dDCHI_L_1, dDCHI_L_2
   real :: dDDCHI_LL_00, dDDCHI_LL_01, dDDCHI_LL_02, &
           dDDCHI_LL_11, dDDCHI_LL_12, dDDCHI_LL_22

   real :: dKODCHI, AdvDBETACHI

   real :: dDGAMTILDE_LLL_000, dDGAMTILDE_LLL_001, dDGAMTILDE_LLL_002, &
           dDGAMTILDE_LLL_010, dDGAMTILDE_LLL_011, dDGAMTILDE_LLL_012, &
           dDGAMTILDE_LLL_020, dDGAMTILDE_LLL_021, dDGAMTILDE_LLL_022, &
           dDGAMTILDE_LLL_110, dDGAMTILDE_LLL_111, dDGAMTILDE_LLL_112, &
           dDGAMTILDE_LLL_120, dDGAMTILDE_LLL_121, dDGAMTILDE_LLL_122, &
           dDGAMTILDE_LLL_220, dDGAMTILDE_LLL_221, dDGAMTILDE_LLL_222

   real :: dDDGAMTILDE_LLLL_0000, dDDGAMTILDE_LLLL_0001, dDDGAMTILDE_LLLL_0002, &
           dDDGAMTILDE_LLLL_0011, dDDGAMTILDE_LLLL_0012, dDDGAMTILDE_LLLL_0022, &
           dDDGAMTILDE_LLLL_0100, dDDGAMTILDE_LLLL_0101, dDDGAMTILDE_LLLL_0102, &
           dDDGAMTILDE_LLLL_0111, dDDGAMTILDE_LLLL_0112, dDDGAMTILDE_LLLL_0122, &
           dDDGAMTILDE_LLLL_0200, dDDGAMTILDE_LLLL_0201, dDDGAMTILDE_LLLL_0202, &
           dDDGAMTILDE_LLLL_0211, dDDGAMTILDE_LLLL_0212, dDDGAMTILDE_LLLL_0222, &
           dDDGAMTILDE_LLLL_1100, dDDGAMTILDE_LLLL_1101, dDDGAMTILDE_LLLL_1102, &
           dDDGAMTILDE_LLLL_1111, dDDGAMTILDE_LLLL_1112, dDDGAMTILDE_LLLL_1122, &
           dDDGAMTILDE_LLLL_1200, dDDGAMTILDE_LLLL_1201, dDDGAMTILDE_LLLL_1202, &
           dDDGAMTILDE_LLLL_1211, dDDGAMTILDE_LLLL_1212, dDDGAMTILDE_LLLL_1222, &
           dDDGAMTILDE_LLLL_2200, dDDGAMTILDE_LLLL_2201, dDDGAMTILDE_LLLL_2202, &
           dDDGAMTILDE_LLLL_2211, dDDGAMTILDE_LLLL_2212, dDDGAMTILDE_LLLL_2222

   real :: dKODGAMTILDE_LL_00, dKODGAMTILDE_LL_01, dKODGAMTILDE_LL_02, &
           dKODGAMTILDE_LL_11, dKODGAMTILDE_LL_12, dKODGAMTILDE_LL_22

   real :: AdvDBETAGAMTILDE_LL_00, AdvDBETAGAMTILDE_LL_01, AdvDBETAGAMTILDE_LL_02, &
           AdvDBETAGAMTILDE_LL_11, AdvDBETAGAMTILDE_LL_12, AdvDBETAGAMTILDE_LL_22

   real :: dDKHAT_L_0, dDKHAT_L_1, dDKHAT_L_2
   real :: dDDKHAT_LL_00, dDDKHAT_LL_01, dDDKHAT_LL_02, &
           dDDKHAT_LL_11, dDDKHAT_LL_12, dDDKHAT_LL_22

   real :: dKODKHAT, AdvDBETAKHAT

   real :: dDATILDE_LLL_000, dDATILDE_LLL_001, dDATILDE_LLL_002, &
           dDATILDE_LLL_010, dDATILDE_LLL_011, dDATILDE_LLL_012, &
           dDATILDE_LLL_020, dDATILDE_LLL_021, dDATILDE_LLL_022, &
           dDATILDE_LLL_110, dDATILDE_LLL_111, dDATILDE_LLL_112, &
           dDATILDE_LLL_120, dDATILDE_LLL_121, dDATILDE_LLL_122, &
           dDATILDE_LLL_220, dDATILDE_LLL_221, dDATILDE_LLL_222

   real :: dDDATILDE_LLLL_0000, dDDATILDE_LLLL_0001, dDDATILDE_LLLL_0002, &
           dDDATILDE_LLLL_0011, dDDATILDE_LLLL_0012, dDDATILDE_LLLL_0022, &
           dDDATILDE_LLLL_0100, dDDATILDE_LLLL_0101, dDDATILDE_LLLL_0102, &
           dDDATILDE_LLLL_0111, dDDATILDE_LLLL_0112, dDDATILDE_LLLL_0122, &
           dDDATILDE_LLLL_0200, dDDATILDE_LLLL_0201, dDDATILDE_LLLL_0202, &
           dDDATILDE_LLLL_0211, dDDATILDE_LLLL_0212, dDDATILDE_LLLL_0222, &
           dDDATILDE_LLLL_1100, dDDATILDE_LLLL_1101, dDDATILDE_LLLL_1102, &
           dDDATILDE_LLLL_1111, dDDATILDE_LLLL_1112, dDDATILDE_LLLL_1122, &
           dDDATILDE_LLLL_1200, dDDATILDE_LLLL_1201, dDDATILDE_LLLL_1202, &
           dDDATILDE_LLLL_1211, dDDATILDE_LLLL_1212, dDDATILDE_LLLL_1222, &
           dDDATILDE_LLLL_2200, dDDATILDE_LLLL_2201, dDDATILDE_LLLL_2202, &
           dDDATILDE_LLLL_2211, dDDATILDE_LLLL_2212, dDDATILDE_LLLL_2222

   real :: dKODATILDE_LL_00, dKODATILDE_LL_01, dKODATILDE_LL_02, &
           dKODATILDE_LL_11, dKODATILDE_LL_12, dKODATILDE_LL_22

   real :: AdvDBETAATILDE_LL_00, AdvDBETAATILDE_LL_01, AdvDBETAATILDE_LL_02, &
           AdvDBETAATILDE_LL_11, AdvDBETAATILDE_LL_12, AdvDBETAATILDE_LL_22

   real :: dDTHETAFUNC_L_0, dDTHETAFUNC_L_1, dDTHETAFUNC_L_2
   real :: dDDTHETAFUNC_LL_00, dDDTHETAFUNC_LL_01, dDDTHETAFUNC_LL_02, &
           dDDTHETAFUNC_LL_11, dDDTHETAFUNC_LL_12, dDDTHETAFUNC_LL_22

   real :: dKODTHETAFUNC, AdvDBETATHETAFUNC

   real :: dDGAMTILDE_UL_00, dDGAMTILDE_UL_01, dDGAMTILDE_UL_02, &
           dDGAMTILDE_UL_10, dDGAMTILDE_UL_11, dDGAMTILDE_UL_12, &
           dDGAMTILDE_UL_20, dDGAMTILDE_UL_21, dDGAMTILDE_UL_22

   real :: dDDGAMTILDE_ULL_000, dDDGAMTILDE_ULL_001, dDDGAMTILDE_ULL_002, &
           dDDGAMTILDE_ULL_011, dDDGAMTILDE_ULL_012, dDDGAMTILDE_ULL_022, &
           dDDGAMTILDE_ULL_100, dDDGAMTILDE_ULL_101, dDDGAMTILDE_ULL_102, &
           dDDGAMTILDE_ULL_111, dDDGAMTILDE_ULL_112, dDDGAMTILDE_ULL_122, &
           dDDGAMTILDE_ULL_200, dDDGAMTILDE_ULL_201, dDDGAMTILDE_ULL_202, &
           dDDGAMTILDE_ULL_211, dDDGAMTILDE_ULL_212, dDDGAMTILDE_ULL_222

   real :: dKODGAMTILDE_U_0, dKODGAMTILDE_U_1, dKODGAMTILDE_U_2

   real :: AdvDBETAGAMTILDE_U_0, AdvDBETAGAMTILDE_U_1, AdvDBETAGAMTILDE_U_2

   real :: dDALPHA_L_0, dDALPHA_L_1, dDALPHA_L_2
   real :: dDDALPHA_LL_00, dDDALPHA_LL_01, dDDALPHA_LL_02, &
           dDDALPHA_LL_11, dDDALPHA_LL_12, dDDALPHA_LL_22

   real :: dKODALPHA, AdvDBETAALPHA

   real :: dDBETA_UL_00, dDBETA_UL_01, dDBETA_UL_02, &
           dDBETA_UL_10, dDBETA_UL_11, dDBETA_UL_12, &
           dDBETA_UL_20, dDBETA_UL_21, dDBETA_UL_22

   real :: dDDBETA_ULL_000, dDDBETA_ULL_001, dDDBETA_ULL_002, &
           dDDBETA_ULL_011, dDDBETA_ULL_012, dDDBETA_ULL_022, &
           dDDBETA_ULL_100, dDDBETA_ULL_101, dDDBETA_ULL_102, &
           dDDBETA_ULL_111, dDDBETA_ULL_112, dDDBETA_ULL_122, &
           dDDBETA_ULL_200, dDDBETA_ULL_201, dDDBETA_ULL_202, &
           dDDBETA_ULL_211, dDDBETA_ULL_212, dDDBETA_ULL_222

   real :: dKODBETA_U_0, dKODBETA_U_1, dKODBETA_U_2
   real :: AdvDBETABETA_U_0, AdvDBETABETA_U_1, AdvDBETABETA_U_2

   real :: invDetGAMTILDE_LL
   real :: invGAMTILDE_UU_00, invGAMTILDE_UU_01, invGAMTILDE_UU_02, &
           invGAMTILDE_UU_11, invGAMTILDE_UU_12, invGAMTILDE_UU_22

   real :: Gamtilde_ULL_000, Gamtilde_ULL_001, Gamtilde_ULL_002, &
           Gamtilde_ULL_011, Gamtilde_ULL_012, Gamtilde_ULL_022, &
           Gamtilde_ULL_100, Gamtilde_ULL_101, Gamtilde_ULL_102, &
           Gamtilde_ULL_111, Gamtilde_ULL_112, Gamtilde_ULL_122, &
           Gamtilde_ULL_200, Gamtilde_ULL_201, Gamtilde_ULL_202, &
           Gamtilde_ULL_211, Gamtilde_ULL_212, Gamtilde_ULL_222

   real :: Gamtilde_LLL_000, Gamtilde_LLL_001, Gamtilde_LLL_002, &
           Gamtilde_LLL_011, Gamtilde_LLL_012, Gamtilde_LLL_022, &
           Gamtilde_LLL_100, Gamtilde_LLL_101, Gamtilde_LLL_102, &
           Gamtilde_LLL_111, Gamtilde_LLL_112, Gamtilde_LLL_122, &
           Gamtilde_LLL_200, Gamtilde_LLL_201, Gamtilde_LLL_202, &
           Gamtilde_LLL_211, Gamtilde_LLL_212, Gamtilde_LLL_222

   real :: invCHI
   real :: Gam_ULL_000, Gam_ULL_001, Gam_ULL_002, &
           Gam_ULL_011, Gam_ULL_012, Gam_ULL_022, &
           Gam_ULL_100, Gam_ULL_101, Gam_ULL_102, &
           Gam_ULL_111, Gam_ULL_112, Gam_ULL_122, &
           Gam_ULL_200, Gam_ULL_201, Gam_ULL_202, &
           Gam_ULL_211, Gam_ULL_212, Gam_ULL_222

   real :: CovDDalpha_LL_00, CovDDalpha_LL_01, CovDDalpha_LL_02, &
           CovDDalpha_LL_11, CovDDalpha_LL_12, CovDDalpha_LL_22

   real :: gam_LL_00, gam_LL_01, gam_LL_02, &
           gam_LL_11, gam_LL_12, gam_LL_22

   real :: invgam_UU_00, invgam_UU_01, invgam_UU_02, &
           invgam_UU_11, invgam_UU_12, invgam_UU_22

   real :: Atilde_UU_00, Atilde_UU_01, Atilde_UU_02, &
           Atilde_UU_11, Atilde_UU_12, Atilde_UU_22

   real :: Atilde_UL_00, Atilde_UL_01, Atilde_UL_02, &
           Atilde_UL_10, Atilde_UL_11, Atilde_UL_12, &
           Atilde_UL_20, Atilde_UL_21, Atilde_UL_22

   real :: thirdCovDDalpha
   real :: CovDDalphaTF_LL_00, CovDDalphaTF_LL_01, CovDDalphaTF_LL_02, &
           CovDDalphaTF_LL_11, CovDDalphaTF_LL_12, CovDDalphaTF_LL_22

   real :: GamtildeD_U_0, GamtildeD_U_1, GamtildeD_U_2

   real :: Rtilde_LL_00, Rtilde_LL_01, Rtilde_LL_02, &
           Rtilde_LL_11, Rtilde_LL_12, Rtilde_LL_22

   real :: CovDtildechi_L_0, CovDtildechi_L_1, CovDtildechi_L_2

   real :: CovDtildeDtildechi_LL_00, CovDtildeDtildechi_LL_01, CovDtildeDtildechi_LL_02, &
           CovDtildeDtildechi_LL_11, CovDtildeDtildechi_LL_12, CovDtildeDtildechi_LL_22

   real :: Rchi_LL_00, Rchi_LL_01, Rchi_LL_02, &
           Rchi_LL_11, Rchi_LL_12, Rchi_LL_22

   real :: R_LL_00, R_LL_01, R_LL_02, &
           R_LL_11, R_LL_12, R_LL_22

   real :: Rsclr

   real :: RTF_LL_00, RTF_LL_01, RTF_LL_02, &
           RTF_LL_11, RTF_LL_12, RTF_LL_22

   real :: CovDtildebeta_UL_00, CovDtildebeta_UL_01, CovDtildebeta_UL_02, &
           CovDtildebeta_UL_10, CovDtildebeta_UL_11, CovDtildebeta_UL_12, &
           CovDtildebeta_UL_20, CovDtildebeta_UL_21, CovDtildebeta_UL_22

   real :: divTildebeta

   real :: mul, mus

   ! For finite differences (stencils given as dimensions)
   real, dimension(-2:2), parameter :: dDi = [1d0/12d0, -2d0/3d0, 0d0, 2d0/3d0, -1d0/12d0]
   real, dimension(-1:3), parameter :: dupDi = [-1d0/4d0, -5d0/6d0, 3d0/2d0, -1d0/2d0, 1d0/12d0]
   real, dimension(-3:1), parameter :: ddnDi = [-1d0/12d0, 1d0/2d0, -3d0/2d0, 5d0/6d0, 1d0/4d0]
   real, dimension(-2:2), parameter :: dDDii = [-1d0/12d0, 4d0/3d0, -5d0/2d0, 4d0/3d0, -1d0/12d0]
   real, dimension(-2:2, -2:2), parameter :: dDDij = reshape([1d0/144d0, -1d0/18d0, 0d0, 1d0/18d0, -1d0/144d0, &
                                                              -1d0/18d0, 4d0/9d0, 0d0, -4d0/9d0, 1d0/18d0, &
                                                              0d0, 0d0, 0d0, 0d0, 0d0, &
                                                              1d0/18d0, -4d0/9d0, 0d0, 4d0/9d0, -1d0/18d0, &
                                                              -1d0/144d0, 1d0/18d0, 0d0, -1d0/18d0, 1d0/144d0], [5, 5])
   real, dimension(-3:3), parameter :: dKODi = [1d0/64d0, -3d0/32d0, 15d0/64d0, -5d0/16d0, 15d0/64d0, -3d0/32d0, 1d0/64d0]

   nullify (rhs); nullify (vars)

   lim = tileDesc%limits

   call tileDesc%deltas(del)
   dx0 = del(IAXIS)
   dx1 = del(JAXIS)
   dx2 = del(KAXIS)

   idx0 = 1d0/dx0
   idx1 = 1d0/dx1
   idx2 = 1d0/dx2

   call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)
   call MoL_getDataPtr(tileDesc, rhs, activeRHS)

   do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
      do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
         do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
            CHI = vars(CHI_VAR, i, j, k)

            GAMTILDE_LL_00 = vars(GAMTILDE_LL_00_VAR, i, j, k)
            GAMTILDE_LL_01 = vars(GAMTILDE_LL_01_VAR, i, j, k)
            GAMTILDE_LL_02 = vars(GAMTILDE_LL_02_VAR, i, j, k)
            GAMTILDE_LL_11 = vars(GAMTILDE_LL_11_VAR, i, j, k)
            GAMTILDE_LL_12 = vars(GAMTILDE_LL_12_VAR, i, j, k)
            GAMTILDE_LL_22 = vars(GAMTILDE_LL_22_VAR, i, j, k)

            KHAT = vars(KHAT_VAR, i, j, k)

            ATILDE_LL_00 = vars(ATILDE_LL_00_VAR, i, j, k)
            ATILDE_LL_01 = vars(ATILDE_LL_01_VAR, i, j, k)
            ATILDE_LL_02 = vars(ATILDE_LL_02_VAR, i, j, k)
            ATILDE_LL_11 = vars(ATILDE_LL_11_VAR, i, j, k)
            ATILDE_LL_12 = vars(ATILDE_LL_12_VAR, i, j, k)
            ATILDE_LL_22 = vars(ATILDE_LL_22_VAR, i, j, k)

            THETAFUNC = vars(THETAFUNC_VAR, i, j, k)

            GAMTILDE_U_0 = vars(GAMTILDE_U_0_VAR, i, j, k)
            GAMTILDE_U_1 = vars(GAMTILDE_U_1_VAR, i, j, k)
            GAMTILDE_U_2 = vars(GAMTILDE_U_2_VAR, i, j, k)

            ALPHA = vars(ALPHA_VAR, i, j, k)

            BETA_U_0 = vars(BETA_U_0_VAR, i, j, k)
            BETA_U_1 = vars(BETA_U_1_VAR, i, j, k)
            BETA_U_2 = vars(BETA_U_2_VAR, i, j, k)

            dDCHI_L_0 = sum(dDi*vars(CHI_VAR, i - 2:i + 2, j, k))*idx0
            dDCHI_L_1 = sum(dDi*vars(CHI_VAR, i, j - 2:j + 2, k))*idx1
            dDCHI_L_2 = sum(dDi*vars(CHI_VAR, i, j, k - 2:k + 2))*idx2

            dDDCHI_LL_00 = sum(dDDii*vars(CHI_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDCHI_LL_01 = sum(dDDij*vars(CHI_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDCHI_LL_02 = sum(dDDij*vars(CHI_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDCHI_LL_11 = sum(dDDii*vars(CHI_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDCHI_LL_12 = sum(dDDij*vars(CHI_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDCHI_LL_22 = sum(dDDii*vars(CHI_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODCHI = sum(dKODi*vars(CHI_VAR, i - 3:i + 3, j, k))*idx0 + &
                      sum(dKODi*vars(CHI_VAR, i, j - 3:j + 3, k))*idx1 + &
                      sum(dKODi*vars(CHI_VAR, i, j, k - 3:k + 3))*idx2

            dDGAMTILDE_LLL_000 = sum(dDi*vars(GAMTILDE_LL_00_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_LLL_010 = sum(dDi*vars(GAMTILDE_LL_01_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_LLL_020 = sum(dDi*vars(GAMTILDE_LL_02_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_LLL_110 = sum(dDi*vars(GAMTILDE_LL_11_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_LLL_120 = sum(dDi*vars(GAMTILDE_LL_12_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_LLL_220 = sum(dDi*vars(GAMTILDE_LL_22_VAR, i - 2:i + 2, j, k))*idx0

            dDGAMTILDE_LLL_001 = sum(dDi*vars(GAMTILDE_LL_00_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_LLL_011 = sum(dDi*vars(GAMTILDE_LL_01_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_LLL_021 = sum(dDi*vars(GAMTILDE_LL_02_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_LLL_111 = sum(dDi*vars(GAMTILDE_LL_11_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_LLL_121 = sum(dDi*vars(GAMTILDE_LL_12_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_LLL_221 = sum(dDi*vars(GAMTILDE_LL_22_VAR, i, j - 2:j + 2, k))*idx1

            dDGAMTILDE_LLL_002 = sum(dDi*vars(GAMTILDE_LL_00_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_LLL_012 = sum(dDi*vars(GAMTILDE_LL_01_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_LLL_022 = sum(dDi*vars(GAMTILDE_LL_02_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_LLL_112 = sum(dDi*vars(GAMTILDE_LL_11_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_LLL_122 = sum(dDi*vars(GAMTILDE_LL_12_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_LLL_222 = sum(dDi*vars(GAMTILDE_LL_22_VAR, i, j, k - 2:k + 2))*idx2

            dDDGAMTILDE_LLLL_0000 = sum(dDDii*vars(GAMTILDE_LL_00_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_LLLL_0100 = sum(dDDii*vars(GAMTILDE_LL_01_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_LLLL_0200 = sum(dDDii*vars(GAMTILDE_LL_02_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_LLLL_1100 = sum(dDDii*vars(GAMTILDE_LL_11_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_LLLL_1200 = sum(dDDii*vars(GAMTILDE_LL_12_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_LLLL_2200 = sum(dDDii*vars(GAMTILDE_LL_22_VAR, i - 2:i + 2, j, k))*idx0**2

            dDDGAMTILDE_LLLL_0001 = sum(dDDij*vars(GAMTILDE_LL_00_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_LLLL_0101 = sum(dDDij*vars(GAMTILDE_LL_01_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_LLLL_0201 = sum(dDDij*vars(GAMTILDE_LL_02_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_LLLL_1101 = sum(dDDij*vars(GAMTILDE_LL_11_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_LLLL_1201 = sum(dDDij*vars(GAMTILDE_LL_12_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_LLLL_2201 = sum(dDDij*vars(GAMTILDE_LL_22_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1

            dDDGAMTILDE_LLLL_0002 = sum(dDDij*vars(GAMTILDE_LL_00_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_LLLL_0102 = sum(dDDij*vars(GAMTILDE_LL_01_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_LLLL_0202 = sum(dDDij*vars(GAMTILDE_LL_02_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_LLLL_1102 = sum(dDDij*vars(GAMTILDE_LL_11_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_LLLL_1202 = sum(dDDij*vars(GAMTILDE_LL_12_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_LLLL_2202 = sum(dDDij*vars(GAMTILDE_LL_22_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2

            dDDGAMTILDE_LLLL_0011 = sum(dDDii*vars(GAMTILDE_LL_00_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_LLLL_0111 = sum(dDDii*vars(GAMTILDE_LL_01_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_LLLL_0211 = sum(dDDii*vars(GAMTILDE_LL_02_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_LLLL_1111 = sum(dDDii*vars(GAMTILDE_LL_11_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_LLLL_1211 = sum(dDDii*vars(GAMTILDE_LL_12_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_LLLL_2211 = sum(dDDii*vars(GAMTILDE_LL_22_VAR, i, j - 2:j + 2, k))*idx1**2

            dDDGAMTILDE_LLLL_0012 = sum(dDDij*vars(GAMTILDE_LL_00_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_LLLL_0112 = sum(dDDij*vars(GAMTILDE_LL_01_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_LLLL_0212 = sum(dDDij*vars(GAMTILDE_LL_02_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_LLLL_1112 = sum(dDDij*vars(GAMTILDE_LL_11_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_LLLL_1212 = sum(dDDij*vars(GAMTILDE_LL_12_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_LLLL_2212 = sum(dDDij*vars(GAMTILDE_LL_22_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2

            dDDGAMTILDE_LLLL_0022 = sum(dDDii*vars(GAMTILDE_LL_00_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_LLLL_0122 = sum(dDDii*vars(GAMTILDE_LL_01_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_LLLL_0222 = sum(dDDii*vars(GAMTILDE_LL_02_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_LLLL_1122 = sum(dDDii*vars(GAMTILDE_LL_11_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_LLLL_1222 = sum(dDDii*vars(GAMTILDE_LL_12_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_LLLL_2222 = sum(dDDii*vars(GAMTILDE_LL_22_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODGAMTILDE_LL_00 = sum(dKODi*vars(GAMTILDE_LL_00_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_00_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_00_VAR, i, j, k - 3:k + 3))*idx2
            dKODGAMTILDE_LL_01 = sum(dKODi*vars(GAMTILDE_LL_01_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_01_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_01_VAR, i, j, k - 3:k + 3))*idx2
            dKODGAMTILDE_LL_02 = sum(dKODi*vars(GAMTILDE_LL_02_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_02_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_02_VAR, i, j, k - 3:k + 3))*idx2
            dKODGAMTILDE_LL_11 = sum(dKODi*vars(GAMTILDE_LL_11_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_11_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_11_VAR, i, j, k - 3:k + 3))*idx2
            dKODGAMTILDE_LL_12 = sum(dKODi*vars(GAMTILDE_LL_12_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_12_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_12_VAR, i, j, k - 3:k + 3))*idx2
            dKODGAMTILDE_LL_22 = sum(dKODi*vars(GAMTILDE_LL_22_VAR, i - 3:i + 3, j, k))*idx0 + &
                                 sum(dKODi*vars(GAMTILDE_LL_22_VAR, i, j - 3:j + 3, k))*idx1 + &
                                 sum(dKODi*vars(GAMTILDE_LL_22_VAR, i, j, k - 3:k + 3))*idx2

            dDKHAT_L_0 = sum(dDi*vars(KHAT_VAR, i - 2:i + 2, j, k))*idx0
            dDKHAT_L_1 = sum(dDi*vars(KHAT_VAR, i, j - 2:j + 2, k))*idx1
            dDKHAT_L_2 = sum(dDi*vars(KHAT_VAR, i, j, k - 2:k + 2))*idx2

            dDDKHAT_LL_00 = sum(dDDii*vars(KHAT_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDKHAT_LL_01 = sum(dDDij*vars(KHAT_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDKHAT_LL_02 = sum(dDDij*vars(KHAT_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDKHAT_LL_11 = sum(dDDii*vars(KHAT_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDKHAT_LL_12 = sum(dDDij*vars(KHAT_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDKHAT_LL_22 = sum(dDDii*vars(KHAT_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODKHAT = sum(dKODi*vars(KHAT_VAR, i - 3:i + 3, j, k))*idx0 + &
                       sum(dKODi*vars(KHAT_VAR, i, j - 3:j + 3, k))*idx1 + &
                       sum(dKODi*vars(KHAT_VAR, i, j, k - 3:k + 3))*idx2

            dDATILDE_LLL_000 = sum(dDi*vars(ATILDE_LL_00_VAR, i - 2:i + 2, j, k))*idx0
            dDATILDE_LLL_010 = sum(dDi*vars(ATILDE_LL_01_VAR, i - 2:i + 2, j, k))*idx0
            dDATILDE_LLL_020 = sum(dDi*vars(ATILDE_LL_02_VAR, i - 2:i + 2, j, k))*idx0
            dDATILDE_LLL_110 = sum(dDi*vars(ATILDE_LL_11_VAR, i - 2:i + 2, j, k))*idx0
            dDATILDE_LLL_120 = sum(dDi*vars(ATILDE_LL_12_VAR, i - 2:i + 2, j, k))*idx0
            dDATILDE_LLL_220 = sum(dDi*vars(ATILDE_LL_22_VAR, i - 2:i + 2, j, k))*idx0

            dDATILDE_LLL_001 = sum(dDi*vars(ATILDE_LL_00_VAR, i, j - 2:j + 2, k))*idx1
            dDATILDE_LLL_011 = sum(dDi*vars(ATILDE_LL_01_VAR, i, j - 2:j + 2, k))*idx1
            dDATILDE_LLL_021 = sum(dDi*vars(ATILDE_LL_02_VAR, i, j - 2:j + 2, k))*idx1
            dDATILDE_LLL_111 = sum(dDi*vars(ATILDE_LL_11_VAR, i, j - 2:j + 2, k))*idx1
            dDATILDE_LLL_121 = sum(dDi*vars(ATILDE_LL_12_VAR, i, j - 2:j + 2, k))*idx1
            dDATILDE_LLL_221 = sum(dDi*vars(ATILDE_LL_22_VAR, i, j - 2:j + 2, k))*idx1

            dDATILDE_LLL_002 = sum(dDi*vars(ATILDE_LL_00_VAR, i, j, k - 2:k + 2))*idx2
            dDATILDE_LLL_012 = sum(dDi*vars(ATILDE_LL_01_VAR, i, j, k - 2:k + 2))*idx2
            dDATILDE_LLL_022 = sum(dDi*vars(ATILDE_LL_02_VAR, i, j, k - 2:k + 2))*idx2
            dDATILDE_LLL_112 = sum(dDi*vars(ATILDE_LL_11_VAR, i, j, k - 2:k + 2))*idx2
            dDATILDE_LLL_122 = sum(dDi*vars(ATILDE_LL_12_VAR, i, j, k - 2:k + 2))*idx2
            dDATILDE_LLL_222 = sum(dDi*vars(ATILDE_LL_22_VAR, i, j, k - 2:k + 2))*idx2

            dDDATILDE_LLLL_0000 = sum(dDDii*vars(ATILDE_LL_00_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDATILDE_LLLL_0100 = sum(dDDii*vars(ATILDE_LL_01_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDATILDE_LLLL_0200 = sum(dDDii*vars(ATILDE_LL_02_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDATILDE_LLLL_1100 = sum(dDDii*vars(ATILDE_LL_11_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDATILDE_LLLL_1200 = sum(dDDii*vars(ATILDE_LL_12_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDATILDE_LLLL_2200 = sum(dDDii*vars(ATILDE_LL_22_VAR, i - 2:i + 2, j, k))*idx0**2

            dDDATILDE_LLLL_0001 = sum(dDDij*vars(ATILDE_LL_00_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDATILDE_LLLL_0101 = sum(dDDij*vars(ATILDE_LL_01_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDATILDE_LLLL_0201 = sum(dDDij*vars(ATILDE_LL_02_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDATILDE_LLLL_1101 = sum(dDDij*vars(ATILDE_LL_11_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDATILDE_LLLL_1201 = sum(dDDij*vars(ATILDE_LL_12_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDATILDE_LLLL_2201 = sum(dDDij*vars(ATILDE_LL_22_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1

            dDDATILDE_LLLL_0002 = sum(dDDij*vars(ATILDE_LL_00_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDATILDE_LLLL_0102 = sum(dDDij*vars(ATILDE_LL_01_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDATILDE_LLLL_0202 = sum(dDDij*vars(ATILDE_LL_02_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDATILDE_LLLL_1102 = sum(dDDij*vars(ATILDE_LL_11_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDATILDE_LLLL_1202 = sum(dDDij*vars(ATILDE_LL_12_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDATILDE_LLLL_2202 = sum(dDDij*vars(ATILDE_LL_22_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2

            dDDATILDE_LLLL_0011 = sum(dDDii*vars(ATILDE_LL_00_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDATILDE_LLLL_0111 = sum(dDDii*vars(ATILDE_LL_01_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDATILDE_LLLL_0211 = sum(dDDii*vars(ATILDE_LL_02_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDATILDE_LLLL_1111 = sum(dDDii*vars(ATILDE_LL_11_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDATILDE_LLLL_1211 = sum(dDDii*vars(ATILDE_LL_12_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDATILDE_LLLL_2211 = sum(dDDii*vars(ATILDE_LL_22_VAR, i, j - 2:j + 2, k))*idx1**2

            dDDATILDE_LLLL_0012 = sum(dDDij*vars(ATILDE_LL_00_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDATILDE_LLLL_0112 = sum(dDDij*vars(ATILDE_LL_01_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDATILDE_LLLL_0212 = sum(dDDij*vars(ATILDE_LL_02_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDATILDE_LLLL_1112 = sum(dDDij*vars(ATILDE_LL_11_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDATILDE_LLLL_1212 = sum(dDDij*vars(ATILDE_LL_12_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDATILDE_LLLL_2212 = sum(dDDij*vars(ATILDE_LL_22_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2

            dDDATILDE_LLLL_0022 = sum(dDDii*vars(ATILDE_LL_00_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDATILDE_LLLL_0122 = sum(dDDii*vars(ATILDE_LL_01_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDATILDE_LLLL_0222 = sum(dDDii*vars(ATILDE_LL_02_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDATILDE_LLLL_1122 = sum(dDDii*vars(ATILDE_LL_11_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDATILDE_LLLL_1222 = sum(dDDii*vars(ATILDE_LL_12_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDATILDE_LLLL_2222 = sum(dDDii*vars(ATILDE_LL_22_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODATILDE_LL_00 = sum(dKODi*vars(ATILDE_LL_00_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_00_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_00_VAR, i, j, k - 3:k + 3))*idx2
            dKODATILDE_LL_01 = sum(dKODi*vars(ATILDE_LL_01_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_01_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_01_VAR, i, j, k - 3:k + 3))*idx2
            dKODATILDE_LL_02 = sum(dKODi*vars(ATILDE_LL_02_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_02_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_02_VAR, i, j, k - 3:k + 3))*idx2
            dKODATILDE_LL_11 = sum(dKODi*vars(ATILDE_LL_11_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_11_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_11_VAR, i, j, k - 3:k + 3))*idx2
            dKODATILDE_LL_12 = sum(dKODi*vars(ATILDE_LL_12_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_12_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_12_VAR, i, j, k - 3:k + 3))*idx2
            dKODATILDE_LL_22 = sum(dKODi*vars(ATILDE_LL_22_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(ATILDE_LL_22_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(ATILDE_LL_22_VAR, i, j, k - 3:k + 3))*idx2

            dDTHETAFUNC_L_0 = sum(dDi*vars(THETAFUNC_VAR, i - 2:i + 2, j, k))*idx0
            dDTHETAFUNC_L_1 = sum(dDi*vars(THETAFUNC_VAR, i, j - 2:j + 2, k))*idx1
            dDTHETAFUNC_L_2 = sum(dDi*vars(THETAFUNC_VAR, i, j, k - 2:k + 2))*idx2

            dDDTHETAFUNC_LL_00 = sum(dDDii*vars(THETAFUNC_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDTHETAFUNC_LL_01 = sum(dDDij*vars(THETAFUNC_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDTHETAFUNC_LL_02 = sum(dDDij*vars(THETAFUNC_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDTHETAFUNC_LL_11 = sum(dDDii*vars(THETAFUNC_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDTHETAFUNC_LL_12 = sum(dDDij*vars(THETAFUNC_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDTHETAFUNC_LL_22 = sum(dDDii*vars(THETAFUNC_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODTHETAFUNC = sum(dKODi*vars(THETAFUNC_VAR, i - 3:i + 3, j, k))*idx0 + &
                            sum(dKODi*vars(THETAFUNC_VAR, i, j - 3:j + 3, k))*idx1 + &
                            sum(dKODi*vars(THETAFUNC_VAR, i, j, k - 3:k + 3))*idx2

            dDGAMTILDE_UL_00 = sum(dDi*vars(GAMTILDE_U_0_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_UL_10 = sum(dDi*vars(GAMTILDE_U_1_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_UL_20 = sum(dDi*vars(GAMTILDE_U_2_VAR, i - 2:i + 2, j, k))*idx0

            dDGAMTILDE_UL_01 = sum(dDi*vars(GAMTILDE_U_0_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_UL_11 = sum(dDi*vars(GAMTILDE_U_1_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_UL_21 = sum(dDi*vars(GAMTILDE_U_2_VAR, i, j - 2:j + 2, k))*idx1

            dDGAMTILDE_UL_02 = sum(dDi*vars(GAMTILDE_U_0_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_UL_12 = sum(dDi*vars(GAMTILDE_U_1_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_UL_22 = sum(dDi*vars(GAMTILDE_U_2_VAR, i, j, k - 2:k + 2))*idx2

            dDDGAMTILDE_ULL_000 = sum(dDDii*vars(GAMTILDE_U_0_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_ULL_100 = sum(dDDii*vars(GAMTILDE_U_1_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDGAMTILDE_ULL_200 = sum(dDDii*vars(GAMTILDE_U_2_VAR, i - 2:i + 2, j, k))*idx0**2

            dDDGAMTILDE_ULL_001 = sum(dDDij*vars(GAMTILDE_U_0_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_ULL_101 = sum(dDDij*vars(GAMTILDE_U_1_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDGAMTILDE_ULL_201 = sum(dDDij*vars(GAMTILDE_U_2_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1

            dDDGAMTILDE_ULL_002 = sum(dDDij*vars(GAMTILDE_U_0_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_ULL_102 = sum(dDDij*vars(GAMTILDE_U_1_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDGAMTILDE_ULL_202 = sum(dDDij*vars(GAMTILDE_U_2_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2

            dDDGAMTILDE_ULL_011 = sum(dDDii*vars(GAMTILDE_U_0_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_ULL_111 = sum(dDDii*vars(GAMTILDE_U_1_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDGAMTILDE_ULL_211 = sum(dDDii*vars(GAMTILDE_U_2_VAR, i, j - 2:j + 2, k))*idx1**2

            dDDGAMTILDE_ULL_012 = sum(dDDij*vars(GAMTILDE_U_0_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_ULL_112 = sum(dDDij*vars(GAMTILDE_U_1_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDGAMTILDE_ULL_212 = sum(dDDij*vars(GAMTILDE_U_2_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2

            dDDGAMTILDE_ULL_022 = sum(dDDii*vars(GAMTILDE_U_0_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_ULL_122 = sum(dDDii*vars(GAMTILDE_U_1_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDGAMTILDE_ULL_222 = sum(dDDii*vars(GAMTILDE_U_2_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODGAMTILDE_U_0 = sum(dKODi*vars(GAMTILDE_U_0_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(GAMTILDE_U_0_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(GAMTILDE_U_0_VAR, i, j, k - 3:k + 3))*idx2

            dKODGAMTILDE_U_1 = sum(dKODi*vars(GAMTILDE_U_1_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(GAMTILDE_U_1_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(GAMTILDE_U_1_VAR, i, j, k - 3:k + 3))*idx2

            dKODGAMTILDE_U_2 = sum(dKODi*vars(GAMTILDE_U_2_VAR, i - 3:i + 3, j, k))*idx0 + &
                               sum(dKODi*vars(GAMTILDE_U_2_VAR, i, j - 3:j + 3, k))*idx1 + &
                               sum(dKODi*vars(GAMTILDE_U_2_VAR, i, j, k - 3:k + 3))*idx2

            dDALPHA_L_0 = sum(dDi*vars(ALPHA_VAR, i - 2:i + 2, j, k))*idx0
            dDALPHA_L_1 = sum(dDi*vars(ALPHA_VAR, i, j - 2:j + 2, k))*idx1
            dDALPHA_L_2 = sum(dDi*vars(ALPHA_VAR, i, j, k - 2:k + 2))*idx2

            dDDALPHA_LL_00 = sum(dDDii*vars(ALPHA_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDALPHA_LL_01 = sum(dDDij*vars(ALPHA_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDALPHA_LL_02 = sum(dDDij*vars(ALPHA_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDALPHA_LL_11 = sum(dDDii*vars(ALPHA_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDALPHA_LL_12 = sum(dDDij*vars(ALPHA_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDALPHA_LL_22 = sum(dDDii*vars(ALPHA_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODALPHA = sum(dKODi*vars(ALPHA_VAR, i - 3:i + 3, j, k))*idx0 + &
                        sum(dKODi*vars(ALPHA_VAR, i, j - 3:j + 3, k))*idx1 + &
                        sum(dKODi*vars(ALPHA_VAR, i, j, k - 3:k + 3))*idx2

            dDBETA_UL_00 = sum(dDi*vars(BETA_U_0_VAR, i - 2:i + 2, j, k))*idx0
            dDBETA_UL_10 = sum(dDi*vars(BETA_U_1_VAR, i - 2:i + 2, j, k))*idx0
            dDBETA_UL_20 = sum(dDi*vars(BETA_U_2_VAR, i - 2:i + 2, j, k))*idx0

            dDBETA_UL_01 = sum(dDi*vars(BETA_U_0_VAR, i, j - 2:j + 2, k))*idx1
            dDBETA_UL_11 = sum(dDi*vars(BETA_U_1_VAR, i, j - 2:j + 2, k))*idx1
            dDBETA_UL_21 = sum(dDi*vars(BETA_U_2_VAR, i, j - 2:j + 2, k))*idx1

            dDBETA_UL_02 = sum(dDi*vars(BETA_U_0_VAR, i, j, k - 2:k + 2))*idx2
            dDBETA_UL_12 = sum(dDi*vars(BETA_U_1_VAR, i, j, k - 2:k + 2))*idx2
            dDBETA_UL_22 = sum(dDi*vars(BETA_U_2_VAR, i, j, k - 2:k + 2))*idx2

            dDDBETA_ULL_000 = sum(dDDii*vars(BETA_U_0_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDBETA_ULL_100 = sum(dDDii*vars(BETA_U_1_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDBETA_ULL_200 = sum(dDDii*vars(BETA_U_2_VAR, i - 2:i + 2, j, k))*idx0**2

            dDDBETA_ULL_001 = sum(dDDij*vars(BETA_U_0_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDBETA_ULL_101 = sum(dDDij*vars(BETA_U_1_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDBETA_ULL_201 = sum(dDDij*vars(BETA_U_2_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1

            dDDBETA_ULL_002 = sum(dDDij*vars(BETA_U_0_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDBETA_ULL_102 = sum(dDDij*vars(BETA_U_1_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDBETA_ULL_202 = sum(dDDij*vars(BETA_U_2_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2

            dDDBETA_ULL_011 = sum(dDDii*vars(BETA_U_0_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDBETA_ULL_111 = sum(dDDii*vars(BETA_U_1_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDBETA_ULL_211 = sum(dDDii*vars(BETA_U_2_VAR, i, j - 2:j + 2, k))*idx1**2

            dDDBETA_ULL_012 = sum(dDDij*vars(BETA_U_0_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDBETA_ULL_112 = sum(dDDij*vars(BETA_U_1_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDBETA_ULL_212 = sum(dDDij*vars(BETA_U_2_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2

            dDDBETA_ULL_022 = sum(dDDii*vars(BETA_U_0_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDBETA_ULL_122 = sum(dDDii*vars(BETA_U_1_VAR, i, j, k - 2:k + 2))*idx2**2
            dDDBETA_ULL_222 = sum(dDDii*vars(BETA_U_2_VAR, i, j, k - 2:k + 2))*idx2**2

            dKODBETA_U_0 = sum(dKODi*vars(BETA_U_0_VAR, i - 3:i + 3, j, k))*idx0 + &
                           sum(dKODi*vars(BETA_U_0_VAR, i, j - 3:j + 3, k))*idx1 + &
                           sum(dKODi*vars(BETA_U_0_VAR, i, j, k - 3:k + 3))*idx2

            dKODBETA_U_1 = sum(dKODi*vars(BETA_U_1_VAR, i - 3:i + 3, j, k))*idx0 + &
                           sum(dKODi*vars(BETA_U_1_VAR, i, j - 3:j + 3, k))*idx1 + &
                           sum(dKODi*vars(BETA_U_1_VAR, i, j, k - 3:k + 3))*idx2

            dKODBETA_U_2 = sum(dKODi*vars(BETA_U_2_VAR, i - 3:i + 3, j, k))*idx0 + &
                           sum(dKODi*vars(BETA_U_2_VAR, i, j - 3:j + 3, k))*idx1 + &
                           sum(dKODi*vars(BETA_U_2_VAR, i, j, k - 3:k + 3))*idx2

            AdvDBETACHI = 0d0
            AdvDBETAKHAT = 0d0
            AdvDBETATHETAFUNC = 0d0
            AdvDBETAALPHA = 0d0
            AdvDBETAGAMTILDE_U_0 = 0d0; AdvDBETAGAMTILDE_U_1 = 0d0; AdvDBETAGAMTILDE_U_2 = 0d0
            AdvDBETABETA_U_0 = 0d0; AdvDBETABETA_U_1 = 0d0; AdvDBETABETA_U_2 = 0d0
            AdvDBETAGAMTILDE_LL_00 = 0d0; AdvDBETAGAMTILDE_LL_01 = 0d0; AdvDBETAGAMTILDE_LL_02 = 0d0
            AdvDBETAGAMTILDE_LL_11 = 0d0; AdvDBETAGAMTILDE_LL_12 = 0d0; AdvDBETAGAMTILDE_LL_22 = 0d0
            AdvDBETAATILDE_LL_00 = 0d0; AdvDBETAATILDE_LL_01 = 0d0; AdvDBETAATILDE_LL_02 = 0d0
            AdvDBETAATILDE_LL_11 = 0d0; AdvDBETAATILDE_LL_12 = 0d0; AdvDBETAATILDE_LL_22 = 0d0

            if (BETA_U_0 .lt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_0*sum(ddnDi*vars(CHI_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_0*sum(ddnDi*vars(KHAT_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_0*sum(ddnDi*vars(THETAFUNC_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_0*sum(ddnDi*vars(ALPHA_VAR, i - 3:i + 1, j, k))*idx0

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_0*sum(ddnDi*vars(GAMTILDE_U_0_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_0*sum(ddnDi*vars(GAMTILDE_U_1_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_0*sum(ddnDi*vars(GAMTILDE_U_2_VAR, i - 3:i + 1, j, k))*idx0

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_0*sum(ddnDi*vars(BETA_U_0_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_0*sum(ddnDi*vars(BETA_U_1_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_0*sum(ddnDi*vars(BETA_U_2_VAR, i - 3:i + 1, j, k))*idx0

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_00_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_01_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_02_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_11_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_12_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_0*sum(ddnDi*vars(GAMTILDE_LL_22_VAR, i - 3:i + 1, j, k))*idx0

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_00_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_01_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_02_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_11_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_12_VAR, i - 3:i + 1, j, k))*idx0
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_22_VAR, i - 3:i + 1, j, k))*idx0
            else if (BETA_U_0 .gt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_0*sum(dupDi*vars(CHI_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_0*sum(dupDi*vars(KHAT_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_0*sum(dupDi*vars(THETAFUNC_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_0*sum(dupDi*vars(ALPHA_VAR, i - 1:i + 3, j, k))*idx0

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_0*sum(dupDi*vars(GAMTILDE_U_0_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_0*sum(dupDi*vars(GAMTILDE_U_1_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_0*sum(dupDi*vars(GAMTILDE_U_2_VAR, i - 1:i + 3, j, k))*idx0

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_0*sum(dupDi*vars(BETA_U_0_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_0*sum(dupDi*vars(BETA_U_1_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_0*sum(dupDi*vars(BETA_U_2_VAR, i - 1:i + 3, j, k))*idx0

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_00_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_01_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_02_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_11_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_12_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_0*sum(dupDi*vars(GAMTILDE_LL_22_VAR, i - 1:i + 3, j, k))*idx0

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_00_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_01_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_02_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_11_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_12_VAR, i - 1:i + 3, j, k))*idx0
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_22_VAR, i - 1:i + 3, j, k))*idx0
            end if

            if (BETA_U_1 .lt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_1*sum(ddnDi*vars(CHI_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_1*sum(ddnDi*vars(KHAT_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_1*sum(ddnDi*vars(THETAFUNC_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_1*sum(ddnDi*vars(ALPHA_VAR, i, j - 3:j + 1, k))*idx1

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_1*sum(ddnDi*vars(GAMTILDE_U_0_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_1*sum(ddnDi*vars(GAMTILDE_U_1_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_1*sum(ddnDi*vars(GAMTILDE_U_2_VAR, i, j - 3:j + 1, k))*idx1

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_1*sum(ddnDi*vars(BETA_U_0_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_1*sum(ddnDi*vars(BETA_U_1_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_1*sum(ddnDi*vars(BETA_U_2_VAR, i, j - 3:j + 1, k))*idx1

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_00_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_01_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_02_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_11_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_12_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_1*sum(ddnDi*vars(GAMTILDE_LL_22_VAR, i, j - 3:j + 1, k))*idx1

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_00_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_01_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_02_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_11_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_12_VAR, i, j - 3:j + 1, k))*idx1
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_0*sum(ddnDi*vars(ATILDE_LL_22_VAR, i, j - 3:j + 1, k))*idx1
            else if (BETA_U_1 .gt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_1*sum(dupDi*vars(CHI_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_1*sum(dupDi*vars(KHAT_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_1*sum(dupDi*vars(THETAFUNC_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_1*sum(dupDi*vars(ALPHA_VAR, i, j - 1:j + 3, k))*idx1

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_1*sum(dupDi*vars(GAMTILDE_U_0_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_1*sum(dupDi*vars(GAMTILDE_U_1_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_1*sum(dupDi*vars(GAMTILDE_U_2_VAR, i, j - 1:j + 3, k))*idx1

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_1*sum(dupDi*vars(BETA_U_0_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_1*sum(dupDi*vars(BETA_U_1_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_1*sum(dupDi*vars(BETA_U_2_VAR, i, j - 1:j + 3, k))*idx1

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_00_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_01_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_02_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_11_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_12_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_1*sum(dupDi*vars(GAMTILDE_LL_22_VAR, i, j - 1:j + 3, k))*idx1

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_00_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_01_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_02_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_11_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_12_VAR, i, j - 1:j + 3, k))*idx1
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_0*sum(dupDi*vars(ATILDE_LL_22_VAR, i, j - 1:j + 3, k))*idx1
            end if

            if (BETA_U_2 .lt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_2*sum(ddnDi*vars(CHI_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_2*sum(ddnDi*vars(KHAT_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_2*sum(ddnDi*vars(THETAFUNC_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_2*sum(ddnDi*vars(ALPHA_VAR, i, j, k - 3:k + 1))*idx2

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_2*sum(ddnDi*vars(GAMTILDE_U_0_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_2*sum(ddnDi*vars(GAMTILDE_U_1_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_2*sum(ddnDi*vars(GAMTILDE_U_2_VAR, i, j, k - 3:k + 1))*idx2

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_2*sum(ddnDi*vars(BETA_U_0_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_2*sum(ddnDi*vars(BETA_U_1_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_2*sum(ddnDi*vars(BETA_U_2_VAR, i, j, k - 3:k + 1))*idx2

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_00_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_01_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_02_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_11_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_12_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_2*sum(ddnDi*vars(GAMTILDE_LL_22_VAR, i, j, k - 3:k + 1))*idx2

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_00_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_01_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_02_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_11_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_12_VAR, i, j, k - 3:k + 1))*idx2
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_2*sum(ddnDi*vars(ATILDE_LL_22_VAR, i, j, k - 3:k + 1))*idx2
            else if (BETA_U_2 .gt. 0d0) then
               AdvDBETACHI = AdvDBETACHI + BETA_U_2*sum(dupDi*vars(CHI_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAKHAT = AdvDBETAKHAT + BETA_U_2*sum(dupDi*vars(KHAT_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETATHETAFUNC = AdvDBETATHETAFUNC + BETA_U_2*sum(dupDi*vars(THETAFUNC_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAALPHA = AdvDBETAALPHA + BETA_U_2*sum(dupDi*vars(ALPHA_VAR, i, j, k - 1:k + 3))*idx2

               AdvDBETAGAMTILDE_U_0 = AdvDBETAGAMTILDE_U_0 &
                                      + BETA_U_2*sum(dupDi*vars(GAMTILDE_U_0_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_U_1 = AdvDBETAGAMTILDE_U_1 &
                                      + BETA_U_2*sum(dupDi*vars(GAMTILDE_U_1_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_U_2 = AdvDBETAGAMTILDE_U_2 &
                                      + BETA_U_2*sum(dupDi*vars(GAMTILDE_U_2_VAR, i, j, k - 1:k + 3))*idx2

               AdvDBETABETA_U_0 = AdvDBETABETA_U_0 &
                                  + BETA_U_2*sum(dupDi*vars(BETA_U_0_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETABETA_U_1 = AdvDBETABETA_U_1 &
                                  + BETA_U_2*sum(dupDi*vars(BETA_U_1_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETABETA_U_2 = AdvDBETABETA_U_2 &
                                  + BETA_U_2*sum(dupDi*vars(BETA_U_2_VAR, i, j, k - 1:k + 3))*idx2

               AdvDBETAGAMTILDE_LL_00 = AdvDBETAGAMTILDE_LL_00 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_00_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_LL_01 = AdvDBETAGAMTILDE_LL_01 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_01_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_LL_02 = AdvDBETAGAMTILDE_LL_02 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_02_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_LL_11 = AdvDBETAGAMTILDE_LL_11 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_11_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_LL_12 = AdvDBETAGAMTILDE_LL_12 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_12_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAGAMTILDE_LL_22 = AdvDBETAGAMTILDE_LL_22 &
                                        + BETA_U_2*sum(dupDi*vars(GAMTILDE_LL_22_VAR, i, j, k - 1:k + 3))*idx2

               AdvDBETAATILDE_LL_00 = AdvDBETAATILDE_LL_00 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_00_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAATILDE_LL_01 = AdvDBETAATILDE_LL_01 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_01_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAATILDE_LL_02 = AdvDBETAATILDE_LL_02 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_02_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAATILDE_LL_11 = AdvDBETAATILDE_LL_11 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_11_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAATILDE_LL_12 = AdvDBETAATILDE_LL_12 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_12_VAR, i, j, k - 1:k + 3))*idx2
               AdvDBETAATILDE_LL_22 = AdvDBETAATILDE_LL_22 &
                                      + BETA_U_2*sum(dupDi*vars(ATILDE_LL_22_VAR, i, j, k - 1:k + 3))*idx2
            end if

            invDetGAMTILDE_LL = GAMTILDE_LL_00*GAMTILDE_LL_11*GAMTILDE_LL_22 &
                                - GAMTILDE_LL_00*GAMTILDE_LL_12**2 &
                                - GAMTILDE_LL_01**2*GAMTILDE_LL_22 &
                                - GAMTILDE_LL_02**2*GAMTILDE_LL_11 &
                                + 2d0*GAMTILDE_LL_01*GAMTILDE_LL_02*GAMTILDE_LL_12

            invGAMTILDE_UU_00 = (GAMTILDE_LL_11*GAMTILDE_LL_22 - GAMTILDE_LL_12**2)*invDetGAMTILDE_LL
            invGAMTILDE_UU_01 = (-GAMTILDE_LL_01*GAMTILDE_LL_22 + GAMTILDE_LL_02*GAMTILDE_LL_12)*invDetGAMTILDE_LL
            invGAMTILDE_UU_02 = (GAMTILDE_LL_01*GAMTILDE_LL_12 - GAMTILDE_LL_02*GAMTILDE_LL_11)*invDetGAMTILDE_LL
            invGAMTILDE_UU_11 = (GAMTILDE_LL_00*GAMTILDE_LL_22 - GAMTILDE_LL_02**2)*invDetGAMTILDE_LL
            invGAMTILDE_UU_12 = (-GAMTILDE_LL_00*GAMTILDE_LL_12 + GAMTILDE_LL_01*GAMTILDE_LL_02)*invDetGAMTILDE_LL
            invGAMTILDE_UU_22 = (GAMTILDE_LL_00*GAMTILDE_LL_11 - GAMTILDE_LL_01**2)*invDetGAMTILDE_LL

            Gamtilde_ULL_000 = 0.5d0*(dDGAMTILDE_LLL_000*invGAMTILDE_UU_00 + &
                                      invGAMTILDE_UU_01*(-dDGAMTILDE_LLL_001 + 2*dDGAMTILDE_LLL_010) + &
                                      invGAMTILDE_UU_02*(-dDGAMTILDE_LLL_002 + 2*dDGAMTILDE_LLL_020))
            Gamtilde_ULL_001 = 0.5d0*(dDGAMTILDE_LLL_001*invGAMTILDE_UU_00 + &
                                      dDGAMTILDE_LLL_110*invGAMTILDE_UU_01 + &
                                      invGAMTILDE_UU_02*(-dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_002 = 0.5d0*(dDGAMTILDE_LLL_002*invGAMTILDE_UU_00 + &
                                      dDGAMTILDE_LLL_220*invGAMTILDE_UU_02 + &
                                      invGAMTILDE_UU_01*(dDGAMTILDE_LLL_012 - dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_011 = 0.5d0*(dDGAMTILDE_LLL_111*invGAMTILDE_UU_01 + &
                                      invGAMTILDE_UU_00*(2*dDGAMTILDE_LLL_011 - dDGAMTILDE_LLL_110) + &
                                      invGAMTILDE_UU_02*(-dDGAMTILDE_LLL_112 + 2*dDGAMTILDE_LLL_121))
            Gamtilde_ULL_012 = 0.5d0*(dDGAMTILDE_LLL_112*invGAMTILDE_UU_01 + &
                                      dDGAMTILDE_LLL_221*invGAMTILDE_UU_02 + &
                                      invGAMTILDE_UU_00*(dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 - &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_022 = 0.5d0*(dDGAMTILDE_LLL_222*invGAMTILDE_UU_02 + &
                                      invGAMTILDE_UU_00*(2*dDGAMTILDE_LLL_022 - dDGAMTILDE_LLL_220) + &
                                      invGAMTILDE_UU_01*(2*dDGAMTILDE_LLL_122 - dDGAMTILDE_LLL_221))
            Gamtilde_ULL_100 = 0.5d0*(dDGAMTILDE_LLL_000*invGAMTILDE_UU_01 + &
                                      invGAMTILDE_UU_11*(-dDGAMTILDE_LLL_001 + 2*dDGAMTILDE_LLL_010) + &
                                      invGAMTILDE_UU_12*(-dDGAMTILDE_LLL_002 + 2*dDGAMTILDE_LLL_020))
            Gamtilde_ULL_101 = 0.5d0*(dDGAMTILDE_LLL_001*invGAMTILDE_UU_01 + &
                                      dDGAMTILDE_LLL_110*invGAMTILDE_UU_11 + &
                                      invGAMTILDE_UU_12*(-dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_102 = 0.5d0*(dDGAMTILDE_LLL_002*invGAMTILDE_UU_01 + &
                                      dDGAMTILDE_LLL_220*invGAMTILDE_UU_12 + &
                                      invGAMTILDE_UU_11*(dDGAMTILDE_LLL_012 - dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_111 = 0.5d0*(dDGAMTILDE_LLL_111*invGAMTILDE_UU_11 + &
                                      invGAMTILDE_UU_01*(2*dDGAMTILDE_LLL_011 - dDGAMTILDE_LLL_110) + &
                                      invGAMTILDE_UU_12*(-dDGAMTILDE_LLL_112 + 2*dDGAMTILDE_LLL_121))
            Gamtilde_ULL_112 = 0.5d0*(dDGAMTILDE_LLL_112*invGAMTILDE_UU_11 + &
                                      dDGAMTILDE_LLL_221*invGAMTILDE_UU_12 + &
                                      invGAMTILDE_UU_01*(dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 - &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_122 = 0.5d0*(dDGAMTILDE_LLL_222*invGAMTILDE_UU_12 + &
                                      invGAMTILDE_UU_01*(2*dDGAMTILDE_LLL_022 - dDGAMTILDE_LLL_220) + &
                                      invGAMTILDE_UU_11*(2*dDGAMTILDE_LLL_122 - dDGAMTILDE_LLL_221))
            Gamtilde_ULL_200 = 0.5d0*(dDGAMTILDE_LLL_000*invGAMTILDE_UU_02 + &
                                      invGAMTILDE_UU_12*(-dDGAMTILDE_LLL_001 + 2*dDGAMTILDE_LLL_010) + &
                                      invGAMTILDE_UU_22*(-dDGAMTILDE_LLL_002 + 2*dDGAMTILDE_LLL_020))
            Gamtilde_ULL_201 = 0.5d0*(dDGAMTILDE_LLL_001*invGAMTILDE_UU_02 + &
                                      dDGAMTILDE_LLL_110*invGAMTILDE_UU_12 + &
                                      invGAMTILDE_UU_22*(-dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_202 = 0.5d0*(dDGAMTILDE_LLL_002*invGAMTILDE_UU_02 + &
                                      dDGAMTILDE_LLL_220*invGAMTILDE_UU_22 + &
                                      invGAMTILDE_UU_12*(dDGAMTILDE_LLL_012 - dDGAMTILDE_LLL_021 + &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_211 = 0.5d0*(dDGAMTILDE_LLL_111*invGAMTILDE_UU_12 + &
                                      invGAMTILDE_UU_02*(2*dDGAMTILDE_LLL_011 - dDGAMTILDE_LLL_110) + &
                                      invGAMTILDE_UU_22*(-dDGAMTILDE_LLL_112 + 2*dDGAMTILDE_LLL_121))
            Gamtilde_ULL_212 = 0.5d0*(dDGAMTILDE_LLL_112*invGAMTILDE_UU_12 + &
                                      dDGAMTILDE_LLL_221*invGAMTILDE_UU_22 + &
                                      invGAMTILDE_UU_02*(dDGAMTILDE_LLL_012 + dDGAMTILDE_LLL_021 - &
                                                         dDGAMTILDE_LLL_120))
            Gamtilde_ULL_222 = 0.5d0*(dDGAMTILDE_LLL_222*invGAMTILDE_UU_22 + &
                                      invGAMTILDE_UU_02*(2*dDGAMTILDE_LLL_022 - dDGAMTILDE_LLL_220) + &
                                      invGAMTILDE_UU_12*(2*dDGAMTILDE_LLL_122 - dDGAMTILDE_LLL_221))

            Gamtilde_LLL_000 = GAMTILDE_LL_00*Gamtilde_ULL_000 + GAMTILDE_LL_01*Gamtilde_ULL_100 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_200
            Gamtilde_LLL_001 = GAMTILDE_LL_00*Gamtilde_ULL_001 + GAMTILDE_LL_01*Gamtilde_ULL_101 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_201
            Gamtilde_LLL_002 = GAMTILDE_LL_00*Gamtilde_ULL_002 + GAMTILDE_LL_01*Gamtilde_ULL_102 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_202
            Gamtilde_LLL_011 = GAMTILDE_LL_00*Gamtilde_ULL_011 + GAMTILDE_LL_01*Gamtilde_ULL_111 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_211
            Gamtilde_LLL_012 = GAMTILDE_LL_00*Gamtilde_ULL_012 + GAMTILDE_LL_01*Gamtilde_ULL_112 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_212
            Gamtilde_LLL_022 = GAMTILDE_LL_00*Gamtilde_ULL_022 + GAMTILDE_LL_01*Gamtilde_ULL_122 &
                               + GAMTILDE_LL_02*Gamtilde_ULL_222
            Gamtilde_LLL_100 = GAMTILDE_LL_01*Gamtilde_ULL_000 + GAMTILDE_LL_11*Gamtilde_ULL_100 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_200
            Gamtilde_LLL_101 = GAMTILDE_LL_01*Gamtilde_ULL_001 + GAMTILDE_LL_11*Gamtilde_ULL_101 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_201
            Gamtilde_LLL_102 = GAMTILDE_LL_01*Gamtilde_ULL_002 + GAMTILDE_LL_11*Gamtilde_ULL_102 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_202
            Gamtilde_LLL_111 = GAMTILDE_LL_01*Gamtilde_ULL_011 + GAMTILDE_LL_11*Gamtilde_ULL_111 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_211
            Gamtilde_LLL_112 = GAMTILDE_LL_01*Gamtilde_ULL_012 + GAMTILDE_LL_11*Gamtilde_ULL_112 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_212
            Gamtilde_LLL_122 = GAMTILDE_LL_01*Gamtilde_ULL_022 + GAMTILDE_LL_11*Gamtilde_ULL_122 &
                               + GAMTILDE_LL_12*Gamtilde_ULL_222
            Gamtilde_LLL_200 = GAMTILDE_LL_02*Gamtilde_ULL_000 + GAMTILDE_LL_12*Gamtilde_ULL_100 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_200
            Gamtilde_LLL_201 = GAMTILDE_LL_02*Gamtilde_ULL_001 + GAMTILDE_LL_12*Gamtilde_ULL_101 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_201
            Gamtilde_LLL_202 = GAMTILDE_LL_02*Gamtilde_ULL_002 + GAMTILDE_LL_12*Gamtilde_ULL_102 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_202
            Gamtilde_LLL_211 = GAMTILDE_LL_02*Gamtilde_ULL_011 + GAMTILDE_LL_12*Gamtilde_ULL_111 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_211
            Gamtilde_LLL_212 = GAMTILDE_LL_02*Gamtilde_ULL_012 + GAMTILDE_LL_12*Gamtilde_ULL_112 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_212
            Gamtilde_LLL_222 = GAMTILDE_LL_02*Gamtilde_ULL_022 + GAMTILDE_LL_12*Gamtilde_ULL_122 &
                               + GAMTILDE_LL_22*Gamtilde_ULL_222

            invCHI = 1d0/CHI
            Gam_ULL_000 = Gamtilde_ULL_000 + 0.5d0*(GAMTILDE_LL_00*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_00*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_00*dDCHI_L_2*invGAMTILDE_UU_02 - 2*dDCHI_L_0)*invCHI
            Gam_ULL_001 = Gamtilde_ULL_001 + 0.5d0*(GAMTILDE_LL_01*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_01*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_01*dDCHI_L_2*invGAMTILDE_UU_02 - dDCHI_L_1)*invCHI
            Gam_ULL_002 = Gamtilde_ULL_002 + 0.5d0*(GAMTILDE_LL_02*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_02*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_02*dDCHI_L_2*invGAMTILDE_UU_02 - dDCHI_L_2)*invCHI
            Gam_ULL_011 = Gamtilde_ULL_011 + 0.5d0*(GAMTILDE_LL_11*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_11*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_11*dDCHI_L_2*invGAMTILDE_UU_02)*invCHI
            Gam_ULL_012 = Gamtilde_ULL_012 + 0.5d0*(GAMTILDE_LL_12*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_12*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_12*dDCHI_L_2*invGAMTILDE_UU_02)*invCHI
            Gam_ULL_022 = Gamtilde_ULL_022 + 0.5d0*(GAMTILDE_LL_22*dDCHI_L_0*invGAMTILDE_UU_00 + &
                                                    GAMTILDE_LL_22*dDCHI_L_1*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_22*dDCHI_L_2*invGAMTILDE_UU_02)*invCHI
            Gam_ULL_100 = Gamtilde_ULL_100 + 0.5d0*(GAMTILDE_LL_00*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_00*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_00*dDCHI_L_2*invGAMTILDE_UU_12)*invCHI
            Gam_ULL_101 = Gamtilde_ULL_101 + 0.5d0*(GAMTILDE_LL_01*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_01*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_01*dDCHI_L_2*invGAMTILDE_UU_12 - dDCHI_L_0)*invCHI
            Gam_ULL_102 = Gamtilde_ULL_102 + 0.5d0*(GAMTILDE_LL_02*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_02*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_02*dDCHI_L_2*invGAMTILDE_UU_12)*invCHI
            Gam_ULL_111 = Gamtilde_ULL_111 + 0.5d0*(GAMTILDE_LL_11*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_11*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_11*dDCHI_L_2*invGAMTILDE_UU_12 - 2*dDCHI_L_1)*invChi
            Gam_ULL_112 = Gamtilde_ULL_112 + 0.5d0*(GAMTILDE_LL_12*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_12*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_12*dDCHI_L_2*invGAMTILDE_UU_12 - dDCHI_L_2)*invCHI
            Gam_ULL_122 = Gamtilde_ULL_122 + 0.5d0*(GAMTILDE_LL_22*dDCHI_L_0*invGAMTILDE_UU_01 + &
                                                    GAMTILDE_LL_22*dDCHI_L_1*invGAMTILDE_UU_11 + &
                                                    GAMTILDE_LL_22*dDCHI_L_2*invGAMTILDE_UU_12)*invCHI
            Gam_ULL_200 = Gamtilde_ULL_200 + 0.5d0*(GAMTILDE_LL_00*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_00*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_00*dDCHI_L_2*invGAMTILDE_UU_22)*invCHI
            Gam_ULL_201 = Gamtilde_ULL_201 + 0.5d0*(GAMTILDE_LL_01*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_01*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_01*dDCHI_L_2*invGAMTILDE_UU_22)*invCHI
            Gam_ULL_202 = Gamtilde_ULL_202 + 0.5d0*(GAMTILDE_LL_02*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_02*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_02*dDCHI_L_2*invGAMTILDE_UU_22 - dDCHI_L_0)*invCHI
            Gam_ULL_211 = Gamtilde_ULL_211 + 0.5d0*(GAMTILDE_LL_11*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_11*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_11*dDCHI_L_2*invGAMTILDE_UU_22)*invCHI
            Gam_ULL_212 = Gamtilde_ULL_212 + 0.5d0*(GAMTILDE_LL_12*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_12*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_12*dDCHI_L_2*invGAMTILDE_UU_22 - dDCHI_L_1)*invCHI
            Gam_ULL_222 = Gamtilde_ULL_222 + 0.5d0*(GAMTILDE_LL_22*dDCHI_L_0*invGAMTILDE_UU_02 + &
                                                    GAMTILDE_LL_22*dDCHI_L_1*invGAMTILDE_UU_12 + &
                                                    GAMTILDE_LL_22*dDCHI_L_2*invGAMTILDE_UU_22 - 2*dDCHI_L_2)*invCHI

            CovDDalpha_LL_00 = -Gam_ULL_000*dDALPHA_L_0 - Gam_ULL_100*dDALPHA_L_1 - &
                               Gam_ULL_200*dDALPHA_L_2 + dDDALPHA_LL_00
            CovDDalpha_LL_01 = -Gam_ULL_001*dDALPHA_L_0 - Gam_ULL_101*dDALPHA_L_1 - &
                               Gam_ULL_201*dDALPHA_L_2 + dDDALPHA_LL_01
            CovDDalpha_LL_02 = -Gam_ULL_002*dDALPHA_L_0 - Gam_ULL_102*dDALPHA_L_1 - &
                               Gam_ULL_202*dDALPHA_L_2 + dDDALPHA_LL_02
            CovDDalpha_LL_11 = -Gam_ULL_011*dDALPHA_L_0 - Gam_ULL_111*dDALPHA_L_1 - &
                               Gam_ULL_211*dDALPHA_L_2 + dDDALPHA_LL_11
            CovDDalpha_LL_12 = -Gam_ULL_012*dDALPHA_L_0 - Gam_ULL_112*dDALPHA_L_1 - &
                               Gam_ULL_212*dDALPHA_L_2 + dDDALPHA_LL_12
            CovDDalpha_LL_22 = -Gam_ULL_022*dDALPHA_L_0 - Gam_ULL_122*dDALPHA_L_1 - &
                               Gam_ULL_222*dDALPHA_L_2 + dDDALPHA_LL_22

            gam_LL_00 = GAMTILDE_LL_00*invCHI
            gam_LL_01 = GAMTILDE_LL_01*invCHI
            gam_LL_02 = GAMTILDE_LL_02*invCHI
            gam_LL_11 = GAMTILDE_LL_11*invCHI
            gam_LL_12 = GAMTILDE_LL_12*invCHI
            gam_LL_22 = GAMTILDE_LL_22*invCHI

            invgam_UU_00 = CHI*invGAMTILDE_UU_00
            invgam_UU_01 = CHI*invGAMTILDE_UU_01
            invgam_UU_02 = CHI*invGAMTILDE_UU_02
            invgam_UU_11 = CHI*invGAMTILDE_UU_11
            invgam_UU_12 = CHI*invGAMTILDE_UU_12
            invgam_UU_22 = CHI*invGAMTILDE_UU_22

            Atilde_UU_00 = ATILDE_LL_00*invGAMTILDE_UU_00**2 + &
                           2*ATILDE_LL_01*invGAMTILDE_UU_00*invGAMTILDE_UU_01 + &
                           2*ATILDE_LL_02*invGAMTILDE_UU_00*invGAMTILDE_UU_02 + &
                           ATILDE_LL_11*invGAMTILDE_UU_01**2 + &
                           2*ATILDE_LL_12*invGAMTILDE_UU_01*invGAMTILDE_UU_02 + &
                           ATILDE_LL_22*invGAMTILDE_UU_02**2
            Atilde_UU_01 = ATILDE_LL_00*invGAMTILDE_UU_00*invGAMTILDE_UU_01 + &
                           ATILDE_LL_01*invGAMTILDE_UU_00*invGAMTILDE_UU_11 + &
                           ATILDE_LL_01*invGAMTILDE_UU_01**2 + &
                           ATILDE_LL_02*invGAMTILDE_UU_00*invGAMTILDE_UU_12 + &
                           ATILDE_LL_02*invGAMTILDE_UU_01*invGAMTILDE_UU_02 + &
                           ATILDE_LL_11*invGAMTILDE_UU_01*invGAMTILDE_UU_11 + &
                           ATILDE_LL_12*invGAMTILDE_UU_01*invGAMTILDE_UU_12 + &
                           ATILDE_LL_12*invGAMTILDE_UU_02*invGAMTILDE_UU_11 + &
                           ATILDE_LL_22*invGAMTILDE_UU_02*invGAMTILDE_UU_12
            Atilde_UU_02 = ATILDE_LL_00*invGAMTILDE_UU_00*invGAMTILDE_UU_02 + &
                           ATILDE_LL_01*invGAMTILDE_UU_00*invGAMTILDE_UU_12 + &
                           ATILDE_LL_01*invGAMTILDE_UU_01*invGAMTILDE_UU_02 + &
                           ATILDE_LL_02*invGAMTILDE_UU_00*invGAMTILDE_UU_22 + &
                           ATILDE_LL_02*invGAMTILDE_UU_02**2 + &
                           ATILDE_LL_11*invGAMTILDE_UU_01*invGAMTILDE_UU_12 + &
                           ATILDE_LL_12*invGAMTILDE_UU_01*invGAMTILDE_UU_22 + &
                           ATILDE_LL_12*invGAMTILDE_UU_02*invGAMTILDE_UU_12 + &
                           ATILDE_LL_22*invGAMTILDE_UU_02*invGAMTILDE_UU_22
            Atilde_UU_11 = ATILDE_LL_00*invGAMTILDE_UU_01**2 + &
                           2*ATILDE_LL_01*invGAMTILDE_UU_01*invGAMTILDE_UU_11 + &
                           2*ATILDE_LL_02*invGAMTILDE_UU_01*invGAMTILDE_UU_12 + &
                           ATILDE_LL_11*invGAMTILDE_UU_11**2 + &
                           2*ATILDE_LL_12*invGAMTILDE_UU_11*invGAMTILDE_UU_12 + &
                           ATILDE_LL_22*invGAMTILDE_UU_12**2
            Atilde_UU_12 = ATILDE_LL_00*invGAMTILDE_UU_01*invGAMTILDE_UU_02 + &
                           ATILDE_LL_01*invGAMTILDE_UU_01*invGAMTILDE_UU_12 + &
                           ATILDE_LL_01*invGAMTILDE_UU_02*invGAMTILDE_UU_11 + &
                           ATILDE_LL_02*invGAMTILDE_UU_01*invGAMTILDE_UU_22 + &
                           ATILDE_LL_02*invGAMTILDE_UU_02*invGAMTILDE_UU_12 + &
                           ATILDE_LL_11*invGAMTILDE_UU_11*invGAMTILDE_UU_12 + &
                           ATILDE_LL_12*invGAMTILDE_UU_11*invGAMTILDE_UU_22 + &
                           ATILDE_LL_12*invGAMTILDE_UU_12**2 + &
                           ATILDE_LL_22*invGAMTILDE_UU_12*invGAMTILDE_UU_22
            Atilde_UU_22 = ATILDE_LL_00*invGAMTILDE_UU_02**2 + &
                           2*ATILDE_LL_01*invGAMTILDE_UU_02*invGAMTILDE_UU_12 + &
                           2*ATILDE_LL_02*invGAMTILDE_UU_02*invGAMTILDE_UU_22 + &
                           ATILDE_LL_11*invGAMTILDE_UU_12**2 + &
                           2*ATILDE_LL_12*invGAMTILDE_UU_12*invGAMTILDE_UU_22 + &
                           ATILDE_LL_22*invGAMTILDE_UU_22**2

            Atilde_UL_00 = ATILDE_LL_00*invGAMTILDE_UU_00 + ATILDE_LL_01*invGAMTILDE_UU_01 + &
                           ATILDE_LL_02*invGAMTILDE_UU_02
            Atilde_UL_01 = ATILDE_LL_01*invGAMTILDE_UU_00 + ATILDE_LL_11*invGAMTILDE_UU_01 + &
                           ATILDE_LL_12*invGAMTILDE_UU_02
            Atilde_UL_02 = ATILDE_LL_02*invGAMTILDE_UU_00 + ATILDE_LL_12*invGAMTILDE_UU_01 + &
                           ATILDE_LL_22*invGAMTILDE_UU_02
            Atilde_UL_10 = ATILDE_LL_00*invGAMTILDE_UU_01 + ATILDE_LL_01*invGAMTILDE_UU_11 + &
                           ATILDE_LL_02*invGAMTILDE_UU_12
            Atilde_UL_11 = ATILDE_LL_01*invGAMTILDE_UU_01 + ATILDE_LL_11*invGAMTILDE_UU_11 + &
                           ATILDE_LL_12*invGAMTILDE_UU_12
            Atilde_UL_12 = ATILDE_LL_02*invGAMTILDE_UU_01 + ATILDE_LL_12*invGAMTILDE_UU_11 + &
                           ATILDE_LL_22*invGAMTILDE_UU_12
            Atilde_UL_20 = ATILDE_LL_00*invGAMTILDE_UU_02 + ATILDE_LL_01*invGAMTILDE_UU_12 + &
                           ATILDE_LL_02*invGAMTILDE_UU_22
            Atilde_UL_21 = ATILDE_LL_01*invGAMTILDE_UU_02 + ATILDE_LL_11*invGAMTILDE_UU_12 + &
                           ATILDE_LL_12*invGAMTILDE_UU_22
            Atilde_UL_22 = ATILDE_LL_02*invGAMTILDE_UU_02 + ATILDE_LL_12*invGAMTILDE_UU_12 + &
                           ATILDE_LL_22*invGAMTILDE_UU_22

            thirdCovDDalpha = (CovDDalpha_LL_00*invgam_UU_00 + 2d0*CovDDalpha_LL_01*invgam_UU_01 + &
                               2d0*CovDDalpha_LL_02*invgam_UU_02 + CovDDalpha_LL_11*invgam_UU_11 + &
                               2d0*CovDDalpha_LL_12*invgam_UU_12 + CovDDalpha_LL_22*invgam_UU_22)/3d0

            CovDDalphaTF_LL_00 = CovDDalpha_LL_00 - gam_LL_00*thirdCovDDalpha
            CovDDalphaTF_LL_01 = CovDDalpha_LL_01 - gam_LL_01*thirdCovDDalpha
            CovDDalphaTF_LL_02 = CovDDalpha_LL_02 - gam_LL_02*thirdCovDDalpha
            CovDDalphaTF_LL_11 = CovDDalpha_LL_11 - gam_LL_11*thirdCovDDalpha
            CovDDalphaTF_LL_12 = CovDDalpha_LL_12 - gam_LL_12*thirdCovDDalpha
            CovDDalphaTF_LL_22 = CovDDalpha_LL_22 - gam_LL_22*thirdCovDDalpha

            GamtildeD_U_0 = Gamtilde_ULL_000*invGAMTILDE_UU_00 + 2*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                            2*Gamtilde_ULL_002*invGAMTILDE_UU_02 + Gamtilde_ULL_011*invGAMTILDE_UU_11 + &
                            2*Gamtilde_ULL_012*invGAMTILDE_UU_12 + Gamtilde_ULL_022*invGAMTILDE_UU_22
            GamtildeD_U_1 = Gamtilde_ULL_100*invGAMTILDE_UU_00 + 2*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                            2*Gamtilde_ULL_102*invGAMTILDE_UU_02 + Gamtilde_ULL_111*invGAMTILDE_UU_11 + &
                            2*Gamtilde_ULL_112*invGAMTILDE_UU_12 + Gamtilde_ULL_122*invGAMTILDE_UU_22
            GamtildeD_U_2 = Gamtilde_ULL_200*invGAMTILDE_UU_00 + 2*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                            2*Gamtilde_ULL_202*invGAMTILDE_UU_02 + Gamtilde_ULL_211*invGAMTILDE_UU_11 + &
                            2*Gamtilde_ULL_212*invGAMTILDE_UU_12 + Gamtilde_ULL_222*invGAMTILDE_UU_22

            Rtilde_LL_00 = GAMTILDE_LL_00*dDGAMTILDE_UL_00 + GAMTILDE_LL_01*dDGAMTILDE_UL_10 + &
                           GAMTILDE_LL_02*dDGAMTILDE_UL_20 + GamtildeD_U_0*Gamtilde_LLL_000 + &
                           GamtildeD_U_1*Gamtilde_LLL_001 + GamtildeD_U_2*Gamtilde_LLL_002 + &
                           3*Gamtilde_LLL_000*Gamtilde_ULL_000*invGAMTILDE_UU_00 + &
                           3*Gamtilde_LLL_000*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_000*Gamtilde_ULL_002*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_001*Gamtilde_ULL_000*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_001*Gamtilde_ULL_001*invGAMTILDE_UU_11 + &
                           3*Gamtilde_LLL_001*Gamtilde_ULL_002*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_001*Gamtilde_ULL_100*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_001*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_001*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_002*Gamtilde_ULL_000*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_002*Gamtilde_ULL_001*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_002*Gamtilde_ULL_002*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_002*Gamtilde_ULL_200*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_002*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_002*Gamtilde_ULL_202*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_011*Gamtilde_ULL_100*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_011*Gamtilde_ULL_101*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_011*Gamtilde_ULL_102*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_100*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_101*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_102*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_200*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_201*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_012*Gamtilde_ULL_202*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_022*Gamtilde_ULL_200*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_022*Gamtilde_ULL_201*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_022*Gamtilde_ULL_202*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_100*Gamtilde_ULL_100*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_100*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_100*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_100*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_101*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_102*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_100*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_101*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_102*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_200*Gamtilde_ULL_200*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_200*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_200*Gamtilde_ULL_202*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_200*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_201*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_202*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_200*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_201*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_202*invGAMTILDE_UU_22 - &
                           dDDGAMTILDE_LLLL_0001*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_0002*invGAMTILDE_UU_02 - &
                           dDDGAMTILDE_LLLL_0012*invGAMTILDE_UU_12 - &
                           0.5d0*dDDGAMTILDE_LLLL_0000*invGAMTILDE_UU_00 - &
                           0.5d0*dDDGAMTILDE_LLLL_0011*invGAMTILDE_UU_11 - &
                           0.5d0*dDDGAMTILDE_LLLL_0022*invGAMTILDE_UU_22

            Rtilde_LL_01 = 0.5d0*(GAMTILDE_LL_00*dDGAMTILDE_UL_01 + GAMTILDE_LL_01*dDGAMTILDE_UL_00 + &
                                  GAMTILDE_LL_01*dDGAMTILDE_UL_11 + GAMTILDE_LL_02*dDGAMTILDE_UL_21 + &
                                  GAMTILDE_LL_11*dDGAMTILDE_UL_10 + GAMTILDE_LL_12*dDGAMTILDE_UL_20 + &
                                  GamtildeD_U_0*Gamtilde_LLL_001 + GamtildeD_U_0*Gamtilde_LLL_100 + &
                                  GamtildeD_U_1*Gamtilde_LLL_011 + GamtildeD_U_1*Gamtilde_LLL_101 + &
                                  GamtildeD_U_2*Gamtilde_LLL_012 + GamtildeD_U_2*Gamtilde_LLL_102) + &
                           Gamtilde_LLL_001*Gamtilde_ULL_000*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_001*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_001*Gamtilde_ULL_002*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_000*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_001*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_002*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_000*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_001*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_002*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_100*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_101*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_111*Gamtilde_ULL_100*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_111*Gamtilde_ULL_101*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_111*Gamtilde_ULL_102*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_100*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_101*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_102*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_200*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_202*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_200*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_201*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_202*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_200*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_201*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_202*invGAMTILDE_UU_22 - 1.0d0/2.0d0 &
                           *dDDGAMTILDE_LLLL_0100*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_0101* &
                           invGAMTILDE_UU_01 - dDDGAMTILDE_LLLL_0102*invGAMTILDE_UU_02 - &
                           1.0d0/2.0d0*dDDGAMTILDE_LLLL_0111*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_0112*invGAMTILDE_UU_12 - 1.0d0/2.0d0* &
                           dDDGAMTILDE_LLLL_0122*invGAMTILDE_UU_22 + invGAMTILDE_UU_00*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_001 + Gamtilde_LLL_100* &
                           Gamtilde_ULL_000) + invGAMTILDE_UU_00*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_101 + Gamtilde_LLL_101*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_002*Gamtilde_ULL_201 + &
                                              Gamtilde_LLL_102*Gamtilde_ULL_200) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_011 + Gamtilde_LLL_100* &
                           Gamtilde_ULL_001) + invGAMTILDE_UU_01*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_001 + Gamtilde_LLL_101*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_111 + &
                                              Gamtilde_LLL_101*Gamtilde_ULL_101) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_211 + Gamtilde_LLL_102* &
                           Gamtilde_ULL_201) + invGAMTILDE_UU_01*(Gamtilde_LLL_011* &
                                                                  Gamtilde_ULL_101 + Gamtilde_LLL_111*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_012*Gamtilde_ULL_201 + &
                                              Gamtilde_LLL_112*Gamtilde_ULL_200) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_012 + Gamtilde_LLL_100* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_02*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_112 + Gamtilde_LLL_101*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_001 + &
                                              Gamtilde_LLL_102*Gamtilde_ULL_000) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_212 + Gamtilde_LLL_102* &
                           Gamtilde_ULL_202) + invGAMTILDE_UU_02*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_101 + Gamtilde_LLL_112*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_022*Gamtilde_ULL_201 + &
                                              Gamtilde_LLL_122*Gamtilde_ULL_200) + invGAMTILDE_UU_11*( &
                           Gamtilde_LLL_001*Gamtilde_ULL_011 + Gamtilde_LLL_101* &
                           Gamtilde_ULL_001) + invGAMTILDE_UU_11*(Gamtilde_LLL_011* &
                                                                  Gamtilde_ULL_111 + Gamtilde_LLL_111*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_012*Gamtilde_ULL_211 + &
                                              Gamtilde_LLL_112*Gamtilde_ULL_201) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_001*Gamtilde_ULL_012 + Gamtilde_LLL_101* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_12*(Gamtilde_LLL_002* &
                                                                  Gamtilde_ULL_011 + Gamtilde_LLL_102*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_011*Gamtilde_ULL_112 + &
                                              Gamtilde_LLL_111*Gamtilde_ULL_102) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_012*Gamtilde_ULL_111 + Gamtilde_LLL_112* &
                           Gamtilde_ULL_101) + invGAMTILDE_UU_12*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_212 + Gamtilde_LLL_112*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_022*Gamtilde_ULL_211 + &
                                              Gamtilde_LLL_122*Gamtilde_ULL_201) + invGAMTILDE_UU_22*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_012 + Gamtilde_LLL_102* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_22*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_112 + Gamtilde_LLL_112*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_022*Gamtilde_ULL_212 + &
                                              Gamtilde_LLL_122*Gamtilde_ULL_202)
            Rtilde_LL_02 = 0.5d0*GAMTILDE_LL_00*dDGAMTILDE_UL_02 + 0.5d0* &
                           GAMTILDE_LL_01*dDGAMTILDE_UL_12 + 0.5d0*GAMTILDE_LL_02* &
                           dDGAMTILDE_UL_00 + 0.5d0*GAMTILDE_LL_02*dDGAMTILDE_UL_22 &
                           + 0.5d0*GAMTILDE_LL_12*dDGAMTILDE_UL_10 + 0.5d0* &
                           GAMTILDE_LL_22*dDGAMTILDE_UL_20 + 0.5d0*GamtildeD_U_0* &
                           Gamtilde_LLL_002 + 0.5d0*GamtildeD_U_0*Gamtilde_LLL_200 + &
                           0.5d0*GamtildeD_U_1*Gamtilde_LLL_012 + 0.5d0* &
                           GamtildeD_U_1*Gamtilde_LLL_201 + 0.5d0*GamtildeD_U_2* &
                           Gamtilde_LLL_022 + 0.5d0*GamtildeD_U_2*Gamtilde_LLL_202 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_000*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_002*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_000*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_001*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_002*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_000*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_001*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_002*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_100*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_100*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_101*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_102*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_100*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_101*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_102*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_200*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_202*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_200*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_201*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_202*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_200*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_201*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_202*invGAMTILDE_UU_22 - 1.0d0/2.0d0 &
                           *dDDGAMTILDE_LLLL_0200*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_0201* &
                           invGAMTILDE_UU_01 - dDDGAMTILDE_LLLL_0202*invGAMTILDE_UU_02 - &
                           1.0d0/2.0d0*dDDGAMTILDE_LLLL_0211*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_0212*invGAMTILDE_UU_12 - 1.0d0/2.0d0* &
                           dDDGAMTILDE_LLLL_0222*invGAMTILDE_UU_22 + invGAMTILDE_UU_00*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_002 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_000) + invGAMTILDE_UU_00*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_201*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_002*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_202*Gamtilde_ULL_200) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_012 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_001) + invGAMTILDE_UU_01*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_002 + Gamtilde_LLL_201*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_112 + &
                                              Gamtilde_LLL_201*Gamtilde_ULL_101) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_212 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_201) + invGAMTILDE_UU_01*(Gamtilde_LLL_011* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_211*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_012*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_212*Gamtilde_ULL_200) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_000*Gamtilde_ULL_022 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_02*(Gamtilde_LLL_001* &
                                                                  Gamtilde_ULL_122 + Gamtilde_LLL_201*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_002 + &
                                              Gamtilde_LLL_202*Gamtilde_ULL_000) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_222 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_202) + invGAMTILDE_UU_02*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_212*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_022*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_200) + invGAMTILDE_UU_11*( &
                           Gamtilde_LLL_001*Gamtilde_ULL_012 + Gamtilde_LLL_201* &
                           Gamtilde_ULL_001) + invGAMTILDE_UU_11*(Gamtilde_LLL_011* &
                                                                  Gamtilde_ULL_112 + Gamtilde_LLL_211*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_012*Gamtilde_ULL_212 + &
                                              Gamtilde_LLL_212*Gamtilde_ULL_201) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_001*Gamtilde_ULL_022 + Gamtilde_LLL_201* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_12*(Gamtilde_LLL_002* &
                                                                  Gamtilde_ULL_012 + Gamtilde_LLL_202*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_011*Gamtilde_ULL_122 + &
                                              Gamtilde_LLL_211*Gamtilde_ULL_102) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_012*Gamtilde_ULL_112 + Gamtilde_LLL_212* &
                           Gamtilde_ULL_101) + invGAMTILDE_UU_12*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_222 + Gamtilde_LLL_212*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_022*Gamtilde_ULL_212 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_201) + invGAMTILDE_UU_22*( &
                           Gamtilde_LLL_002*Gamtilde_ULL_022 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_002) + invGAMTILDE_UU_22*(Gamtilde_LLL_012* &
                                                                  Gamtilde_ULL_122 + Gamtilde_LLL_212*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_022*Gamtilde_ULL_222 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_202)
            Rtilde_LL_11 = GAMTILDE_LL_01*dDGAMTILDE_UL_01 + GAMTILDE_LL_11*dDGAMTILDE_UL_11 &
                           + GAMTILDE_LL_12*dDGAMTILDE_UL_21 + GamtildeD_U_0* &
                           Gamtilde_LLL_101 + GamtildeD_U_1*Gamtilde_LLL_111 + GamtildeD_U_2 &
                           *Gamtilde_LLL_112 + Gamtilde_LLL_001*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_00 + Gamtilde_LLL_001*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_001*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_011*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_011*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_11 + Gamtilde_LLL_011*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_012*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_012*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_012*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_22 + 2*Gamtilde_LLL_100*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_00 + 2*Gamtilde_LLL_100*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_100*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_101*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_101*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_11 + 2*Gamtilde_LLL_101*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_12 + 3*Gamtilde_LLL_101*Gamtilde_ULL_101* &
                           invGAMTILDE_UU_00 + 3*Gamtilde_LLL_101*Gamtilde_ULL_111* &
                           invGAMTILDE_UU_01 + 3*Gamtilde_LLL_101*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_102*Gamtilde_ULL_001* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_102*Gamtilde_ULL_011* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_102*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_22 + 2*Gamtilde_LLL_102*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_00 + 2*Gamtilde_LLL_102*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_102*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_02 + 3*Gamtilde_LLL_111*Gamtilde_ULL_101* &
                           invGAMTILDE_UU_01 + 3*Gamtilde_LLL_111*Gamtilde_ULL_111* &
                           invGAMTILDE_UU_11 + 3*Gamtilde_LLL_111*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_12 + 3*Gamtilde_LLL_112*Gamtilde_ULL_101* &
                           invGAMTILDE_UU_02 + 3*Gamtilde_LLL_112*Gamtilde_ULL_111* &
                           invGAMTILDE_UU_12 + 3*Gamtilde_LLL_112*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_22 + 2*Gamtilde_LLL_112*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_112*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_11 + 2*Gamtilde_LLL_112*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_122*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_122*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_122*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_22 + Gamtilde_LLL_201*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_00 + Gamtilde_LLL_201*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_201*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_211*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_211*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_11 + Gamtilde_LLL_211*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_212*Gamtilde_ULL_201* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_212*Gamtilde_ULL_211* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_212*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_22 - 1.0d0/2.0d0*dDDGAMTILDE_LLLL_1100* &
                           invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_1101*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_1102*invGAMTILDE_UU_02 - 1.0d0/2.0d0* &
                           dDDGAMTILDE_LLLL_1111*invGAMTILDE_UU_11 - dDDGAMTILDE_LLLL_1112* &
                           invGAMTILDE_UU_12 - 1.0d0/2.0d0*dDDGAMTILDE_LLLL_1122* &
                           invGAMTILDE_UU_22
            Rtilde_LL_12 = 0.5d0*GAMTILDE_LL_01*dDGAMTILDE_UL_02 + 0.5d0* &
                           GAMTILDE_LL_02*dDGAMTILDE_UL_01 + 0.5d0*GAMTILDE_LL_11* &
                           dDGAMTILDE_UL_12 + 0.5d0*GAMTILDE_LL_12*dDGAMTILDE_UL_11 &
                           + 0.5d0*GAMTILDE_LL_12*dDGAMTILDE_UL_22 + 0.5d0* &
                           GAMTILDE_LL_22*dDGAMTILDE_UL_21 + 0.5d0*GamtildeD_U_0* &
                           Gamtilde_LLL_102 + 0.5d0*GamtildeD_U_0*Gamtilde_LLL_201 + &
                           0.5d0*GamtildeD_U_1*Gamtilde_LLL_112 + 0.5d0* &
                           GamtildeD_U_1*Gamtilde_LLL_211 + 0.5d0*GamtildeD_U_2* &
                           Gamtilde_LLL_122 + 0.5d0*GamtildeD_U_2*Gamtilde_LLL_212 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_001*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_011*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_012*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_011*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_012*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_001*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_011*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_012*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_101*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_111*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_112*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_111*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_112*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_101*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_111*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_112*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_201*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_211*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_202*Gamtilde_ULL_212*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_211*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_212*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_201*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_211*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_222*Gamtilde_ULL_212*invGAMTILDE_UU_22 - 1.0d0/2.0d0 &
                           *dDDGAMTILDE_LLLL_1200*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_1201* &
                           invGAMTILDE_UU_01 - dDDGAMTILDE_LLLL_1202*invGAMTILDE_UU_02 - &
                           1.0d0/2.0d0*dDDGAMTILDE_LLLL_1211*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_1212*invGAMTILDE_UU_12 - 1.0d0/2.0d0* &
                           dDDGAMTILDE_LLLL_1222*invGAMTILDE_UU_22 + invGAMTILDE_UU_00*( &
                           Gamtilde_LLL_100*Gamtilde_ULL_002 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_001) + invGAMTILDE_UU_00*(Gamtilde_LLL_101* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_201*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_102*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_202*Gamtilde_ULL_201) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_100*Gamtilde_ULL_012 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_011) + invGAMTILDE_UU_01*(Gamtilde_LLL_101* &
                                                                  Gamtilde_ULL_002 + Gamtilde_LLL_201*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_101*Gamtilde_ULL_112 + &
                                              Gamtilde_LLL_201*Gamtilde_ULL_111) + invGAMTILDE_UU_01*( &
                           Gamtilde_LLL_102*Gamtilde_ULL_212 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_211) + invGAMTILDE_UU_01*(Gamtilde_LLL_111* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_211*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_112*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_212*Gamtilde_ULL_201) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_100*Gamtilde_ULL_022 + Gamtilde_LLL_200* &
                           Gamtilde_ULL_012) + invGAMTILDE_UU_02*(Gamtilde_LLL_101* &
                                                                  Gamtilde_ULL_122 + Gamtilde_LLL_201*Gamtilde_ULL_112) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_102*Gamtilde_ULL_002 + &
                                              Gamtilde_LLL_202*Gamtilde_ULL_001) + invGAMTILDE_UU_02*( &
                           Gamtilde_LLL_102*Gamtilde_ULL_222 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_212) + invGAMTILDE_UU_02*(Gamtilde_LLL_112* &
                                                                  Gamtilde_ULL_102 + Gamtilde_LLL_212*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_122*Gamtilde_ULL_202 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_201) + invGAMTILDE_UU_11*( &
                           Gamtilde_LLL_101*Gamtilde_ULL_012 + Gamtilde_LLL_201* &
                           Gamtilde_ULL_011) + invGAMTILDE_UU_11*(Gamtilde_LLL_111* &
                                                                  Gamtilde_ULL_112 + Gamtilde_LLL_211*Gamtilde_ULL_111) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_112*Gamtilde_ULL_212 + &
                                              Gamtilde_LLL_212*Gamtilde_ULL_211) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_101*Gamtilde_ULL_022 + Gamtilde_LLL_201* &
                           Gamtilde_ULL_012) + invGAMTILDE_UU_12*(Gamtilde_LLL_102* &
                                                                  Gamtilde_ULL_012 + Gamtilde_LLL_202*Gamtilde_ULL_011) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_111*Gamtilde_ULL_122 + &
                                              Gamtilde_LLL_211*Gamtilde_ULL_112) + invGAMTILDE_UU_12*( &
                           Gamtilde_LLL_112*Gamtilde_ULL_112 + Gamtilde_LLL_212* &
                           Gamtilde_ULL_111) + invGAMTILDE_UU_12*(Gamtilde_LLL_112* &
                                                                  Gamtilde_ULL_222 + Gamtilde_LLL_212*Gamtilde_ULL_212) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_122*Gamtilde_ULL_212 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_211) + invGAMTILDE_UU_22*( &
                           Gamtilde_LLL_102*Gamtilde_ULL_022 + Gamtilde_LLL_202* &
                           Gamtilde_ULL_012) + invGAMTILDE_UU_22*(Gamtilde_LLL_112* &
                                                                  Gamtilde_ULL_122 + Gamtilde_LLL_212*Gamtilde_ULL_112) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_122*Gamtilde_ULL_222 + &
                                              Gamtilde_LLL_222*Gamtilde_ULL_212)
            Rtilde_LL_22 = GAMTILDE_LL_02*dDGAMTILDE_UL_02 + GAMTILDE_LL_12*dDGAMTILDE_UL_12 &
                           + GAMTILDE_LL_22*dDGAMTILDE_UL_22 + GamtildeD_U_0* &
                           Gamtilde_LLL_202 + GamtildeD_U_1*Gamtilde_LLL_212 + GamtildeD_U_2 &
                           *Gamtilde_LLL_222 + Gamtilde_LLL_002*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_00 + Gamtilde_LLL_002*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_002*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_012*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_012*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_11 + Gamtilde_LLL_012*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_022*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_022*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_022*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_22 + Gamtilde_LLL_102*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_00 + Gamtilde_LLL_102*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_102*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_112*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_01 + Gamtilde_LLL_112*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_11 + Gamtilde_LLL_112*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_122*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_02 + Gamtilde_LLL_122*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_12 + Gamtilde_LLL_122*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_22 + 2*Gamtilde_LLL_200*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_00 + 2*Gamtilde_LLL_200*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_200*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_201*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_201*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_11 + 2*Gamtilde_LLL_201*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_201*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_00 + 2*Gamtilde_LLL_201*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_201*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_202*Gamtilde_ULL_002* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_202*Gamtilde_ULL_012* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_202*Gamtilde_ULL_022* &
                           invGAMTILDE_UU_22 + 3*Gamtilde_LLL_202*Gamtilde_ULL_202* &
                           invGAMTILDE_UU_00 + 3*Gamtilde_LLL_202*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_01 + 3*Gamtilde_LLL_202*Gamtilde_ULL_222* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_211*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_01 + 2*Gamtilde_LLL_211*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_11 + 2*Gamtilde_LLL_211*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_212*Gamtilde_ULL_102* &
                           invGAMTILDE_UU_02 + 2*Gamtilde_LLL_212*Gamtilde_ULL_112* &
                           invGAMTILDE_UU_12 + 2*Gamtilde_LLL_212*Gamtilde_ULL_122* &
                           invGAMTILDE_UU_22 + 3*Gamtilde_LLL_212*Gamtilde_ULL_202* &
                           invGAMTILDE_UU_01 + 3*Gamtilde_LLL_212*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_11 + 3*Gamtilde_LLL_212*Gamtilde_ULL_222* &
                           invGAMTILDE_UU_12 + 3*Gamtilde_LLL_222*Gamtilde_ULL_202* &
                           invGAMTILDE_UU_02 + 3*Gamtilde_LLL_222*Gamtilde_ULL_212* &
                           invGAMTILDE_UU_12 + 3*Gamtilde_LLL_222*Gamtilde_ULL_222* &
                           invGAMTILDE_UU_22 - 1.0d0/2.0d0*dDDGAMTILDE_LLLL_2200* &
                           invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_2201*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_2202*invGAMTILDE_UU_02 - 1.0d0/2.0d0* &
                           dDDGAMTILDE_LLLL_2211*invGAMTILDE_UU_11 - dDDGAMTILDE_LLLL_2212* &
                           invGAMTILDE_UU_12 - 1.0d0/2.0d0*dDDGAMTILDE_LLLL_2222* &
                           invGAMTILDE_UU_22

            CovDtildechi_L_0 = dDCHI_L_0
            CovDtildechi_L_1 = dDCHI_L_1
            CovDtildechi_L_2 = dDCHI_L_2

            CovDtildeDtildechi_LL_00 = -Gamtilde_ULL_000*dDCHI_L_0 - Gamtilde_ULL_100*dDCHI_L_1 - &
                                       Gamtilde_ULL_200*dDCHI_L_2 + dDDCHI_LL_00
            CovDtildeDtildechi_LL_01 = -Gamtilde_ULL_001*dDCHI_L_0 - Gamtilde_ULL_101*dDCHI_L_1 - &
                                       Gamtilde_ULL_201*dDCHI_L_2 + dDDCHI_LL_01
            CovDtildeDtildechi_LL_02 = -Gamtilde_ULL_002*dDCHI_L_0 - Gamtilde_ULL_102*dDCHI_L_1 - &
                                       Gamtilde_ULL_202*dDCHI_L_2 + dDDCHI_LL_02
            CovDtildeDtildechi_LL_11 = -Gamtilde_ULL_011*dDCHI_L_0 - Gamtilde_ULL_111*dDCHI_L_1 - &
                                       Gamtilde_ULL_211*dDCHI_L_2 + dDDCHI_LL_11
            CovDtildeDtildechi_LL_12 = -Gamtilde_ULL_012*dDCHI_L_0 - Gamtilde_ULL_112*dDCHI_L_1 - &
                                       Gamtilde_ULL_212*dDCHI_L_2 + dDDCHI_LL_12
            CovDtildeDtildechi_LL_22 = -Gamtilde_ULL_022*dDCHI_L_0 - Gamtilde_ULL_122*dDCHI_L_1 - &
                                       Gamtilde_ULL_222*dDCHI_L_2 + dDDCHI_LL_22

            Rchi_LL_00 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_00* &
                         invGAMTILDE_UU_00/CHI + 0.5d0*CovDtildeDtildechi_LL_00/ &
                         CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_00*invGAMTILDE_UU_01/ &
                         CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_00*invGAMTILDE_UU_02/ &
                         CHI + 0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_00* &
                         invGAMTILDE_UU_11/CHI + CovDtildeDtildechi_LL_12*GAMTILDE_LL_00* &
                         invGAMTILDE_UU_12/CHI + 0.5d0*CovDtildeDtildechi_LL_22* &
                         GAMTILDE_LL_00*invGAMTILDE_UU_22/CHI - 3.0d0/4.0d0* &
                         CovDtildechi_L_0**2*GAMTILDE_LL_00*invGAMTILDE_UU_00/CHI**2 - &
                         1.0d0/4.0d0*CovDtildechi_L_0**2/CHI**2 - 3.0d0/2.0d0* &
                         CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_00* &
                         invGAMTILDE_UU_01/CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0* &
                         CovDtildechi_L_2*GAMTILDE_LL_00*invGAMTILDE_UU_02/CHI**2 - 3.0d0/ &
                         4.0d0*CovDtildechi_L_1**2*GAMTILDE_LL_00*invGAMTILDE_UU_11/CHI**2 &
                         - 3.0d0/2.0d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_00* &
                         invGAMTILDE_UU_12/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_2**2* &
                         GAMTILDE_LL_00*invGAMTILDE_UU_22/CHI**2
            Rchi_LL_01 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_00/CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_01/CHI + 0.5d0*CovDtildeDtildechi_LL_01/ &
                         CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_01*invGAMTILDE_UU_02/ &
                         CHI + 0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_11/CHI + CovDtildeDtildechi_LL_12*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_12/CHI + 0.5d0*CovDtildeDtildechi_LL_22* &
                         GAMTILDE_LL_01*invGAMTILDE_UU_22/CHI - 3.0d0/4.0d0* &
                         CovDtildechi_L_0**2*GAMTILDE_LL_01*invGAMTILDE_UU_00/CHI**2 - &
                         3.0d0/2.0d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_01/CHI**2 - 1.0d0/4.0d0*CovDtildechi_L_0* &
                         CovDtildechi_L_1/CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0* &
                         CovDtildechi_L_2*GAMTILDE_LL_01*invGAMTILDE_UU_02/CHI**2 - 3.0d0/ &
                         4.0d0*CovDtildechi_L_1**2*GAMTILDE_LL_01*invGAMTILDE_UU_11/CHI**2 &
                         - 3.0d0/2.0d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_01* &
                         invGAMTILDE_UU_12/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_2**2* &
                         GAMTILDE_LL_01*invGAMTILDE_UU_22/CHI**2
            Rchi_LL_02 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_00/CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_01/CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_02/CHI + 0.5d0*CovDtildeDtildechi_LL_02/ &
                         CHI + 0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_11/CHI + CovDtildeDtildechi_LL_12*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_12/CHI + 0.5d0*CovDtildeDtildechi_LL_22* &
                         GAMTILDE_LL_02*invGAMTILDE_UU_22/CHI - 3.0d0/4.0d0* &
                         CovDtildechi_L_0**2*GAMTILDE_LL_02*invGAMTILDE_UU_00/CHI**2 - &
                         3.0d0/2.0d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_01/CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0* &
                         CovDtildechi_L_2*GAMTILDE_LL_02*invGAMTILDE_UU_02/CHI**2 - 1.0d0/ &
                         4.0d0*CovDtildechi_L_0*CovDtildechi_L_2/CHI**2 - 3.0d0/4.0d0* &
                         CovDtildechi_L_1**2*GAMTILDE_LL_02*invGAMTILDE_UU_11/CHI**2 - &
                         3.0d0/2.0d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_02* &
                         invGAMTILDE_UU_12/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_2**2* &
                         GAMTILDE_LL_02*invGAMTILDE_UU_22/CHI**2
            Rchi_LL_11 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_11* &
                         invGAMTILDE_UU_00/CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_11* &
                         invGAMTILDE_UU_01/CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_11* &
                         invGAMTILDE_UU_02/CHI + 0.5d0*CovDtildeDtildechi_LL_11* &
                         GAMTILDE_LL_11*invGAMTILDE_UU_11/CHI + 0.5d0* &
                         CovDtildeDtildechi_LL_11/CHI + CovDtildeDtildechi_LL_12* &
                         GAMTILDE_LL_11*invGAMTILDE_UU_12/CHI + 0.5d0* &
                         CovDtildeDtildechi_LL_22*GAMTILDE_LL_11*invGAMTILDE_UU_22/CHI - &
                         3.0d0/4.0d0*CovDtildechi_L_0**2*GAMTILDE_LL_11*invGAMTILDE_UU_00/ &
                         CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0*CovDtildechi_L_1* &
                         GAMTILDE_LL_11*invGAMTILDE_UU_01/CHI**2 - 3.0d0/2.0d0* &
                         CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_11* &
                         invGAMTILDE_UU_02/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_1**2* &
                         GAMTILDE_LL_11*invGAMTILDE_UU_11/CHI**2 - 1.0d0/4.0d0* &
                         CovDtildechi_L_1**2/CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_1* &
                         CovDtildechi_L_2*GAMTILDE_LL_11*invGAMTILDE_UU_12/CHI**2 - 3.0d0/ &
                         4.0d0*CovDtildechi_L_2**2*GAMTILDE_LL_11*invGAMTILDE_UU_22/CHI**2
            Rchi_LL_12 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_12* &
                         invGAMTILDE_UU_00/CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_12* &
                         invGAMTILDE_UU_01/CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_12* &
                         invGAMTILDE_UU_02/CHI + 0.5d0*CovDtildeDtildechi_LL_11* &
                         GAMTILDE_LL_12*invGAMTILDE_UU_11/CHI + CovDtildeDtildechi_LL_12* &
                         GAMTILDE_LL_12*invGAMTILDE_UU_12/CHI + 0.5d0* &
                         CovDtildeDtildechi_LL_12/CHI + 0.5d0* &
                         CovDtildeDtildechi_LL_22*GAMTILDE_LL_12*invGAMTILDE_UU_22/CHI - &
                         3.0d0/4.0d0*CovDtildechi_L_0**2*GAMTILDE_LL_12*invGAMTILDE_UU_00/ &
                         CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0*CovDtildechi_L_1* &
                         GAMTILDE_LL_12*invGAMTILDE_UU_01/CHI**2 - 3.0d0/2.0d0* &
                         CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_12* &
                         invGAMTILDE_UU_02/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_1**2* &
                         GAMTILDE_LL_12*invGAMTILDE_UU_11/CHI**2 - 3.0d0/2.0d0* &
                         CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_12* &
                         invGAMTILDE_UU_12/CHI**2 - 1.0d0/4.0d0*CovDtildechi_L_1* &
                         CovDtildechi_L_2/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_2**2* &
                         GAMTILDE_LL_12*invGAMTILDE_UU_22/CHI**2
            Rchi_LL_22 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_22* &
                         invGAMTILDE_UU_00/CHI + CovDtildeDtildechi_LL_01*GAMTILDE_LL_22* &
                         invGAMTILDE_UU_01/CHI + CovDtildeDtildechi_LL_02*GAMTILDE_LL_22* &
                         invGAMTILDE_UU_02/CHI + 0.5d0*CovDtildeDtildechi_LL_11* &
                         GAMTILDE_LL_22*invGAMTILDE_UU_11/CHI + CovDtildeDtildechi_LL_12* &
                         GAMTILDE_LL_22*invGAMTILDE_UU_12/CHI + 0.5d0* &
                         CovDtildeDtildechi_LL_22*GAMTILDE_LL_22*invGAMTILDE_UU_22/CHI + ( &
                         1.0d0/2.0d0)*CovDtildeDtildechi_LL_22/CHI - 3.0d0/4.0d0* &
                         CovDtildechi_L_0**2*GAMTILDE_LL_22*invGAMTILDE_UU_00/CHI**2 - &
                         3.0d0/2.0d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_22* &
                         invGAMTILDE_UU_01/CHI**2 - 3.0d0/2.0d0*CovDtildechi_L_0* &
                         CovDtildechi_L_2*GAMTILDE_LL_22*invGAMTILDE_UU_02/CHI**2 - 3.0d0/ &
                         4.0d0*CovDtildechi_L_1**2*GAMTILDE_LL_22*invGAMTILDE_UU_11/CHI**2 &
                         - 3.0d0/2.0d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_22* &
                         invGAMTILDE_UU_12/CHI**2 - 3.0d0/4.0d0*CovDtildechi_L_2**2* &
                         GAMTILDE_LL_22*invGAMTILDE_UU_22/CHI**2 - 1.0d0/4.0d0* &
                         CovDtildechi_L_2**2/CHI**2

            R_LL_00 = Rchi_LL_00 + Rtilde_LL_00
            R_LL_01 = Rchi_LL_01 + Rtilde_LL_01
            R_LL_02 = Rchi_LL_02 + Rtilde_LL_02
            R_LL_11 = Rchi_LL_11 + Rtilde_LL_11
            R_LL_12 = Rchi_LL_12 + Rtilde_LL_12
            R_LL_22 = Rchi_LL_22 + Rtilde_LL_22

            Rsclr = R_LL_00*invgam_UU_00 + 2*R_LL_01*invgam_UU_01 + 2*R_LL_02* &
                    invgam_UU_02 + R_LL_11*invgam_UU_11 + 2*R_LL_12*invgam_UU_12 + &
                    R_LL_22*invgam_UU_22

            RTF_LL_00 = R_LL_00 - 1.0d0/3.0d0*Rsclr*gam_LL_00
            RTF_LL_01 = R_LL_01 - 1.0d0/3.0d0*Rsclr*gam_LL_01
            RTF_LL_02 = R_LL_02 - 1.0d0/3.0d0*Rsclr*gam_LL_02
            RTF_LL_11 = R_LL_11 - 1.0d0/3.0d0*Rsclr*gam_LL_11
            RTF_LL_12 = R_LL_12 - 1.0d0/3.0d0*Rsclr*gam_LL_12
            RTF_LL_22 = R_LL_22 - 1.0d0/3.0d0*Rsclr*gam_LL_22

            CovDtildebeta_UL_00 = BETA_U_0*Gamtilde_ULL_000 + BETA_U_1*Gamtilde_ULL_001 + BETA_U_2* &
                                  Gamtilde_ULL_002 + dDBETA_UL_00
            CovDtildebeta_UL_01 = BETA_U_0*Gamtilde_ULL_001 + BETA_U_1*Gamtilde_ULL_011 + BETA_U_2* &
                                  Gamtilde_ULL_012 + dDBETA_UL_01
            CovDtildebeta_UL_02 = BETA_U_0*Gamtilde_ULL_002 + BETA_U_1*Gamtilde_ULL_012 + BETA_U_2* &
                                  Gamtilde_ULL_022 + dDBETA_UL_02
            CovDtildebeta_UL_10 = BETA_U_0*Gamtilde_ULL_100 + BETA_U_1*Gamtilde_ULL_101 + BETA_U_2* &
                                  Gamtilde_ULL_102 + dDBETA_UL_10
            CovDtildebeta_UL_11 = BETA_U_0*Gamtilde_ULL_101 + BETA_U_1*Gamtilde_ULL_111 + BETA_U_2* &
                                  Gamtilde_ULL_112 + dDBETA_UL_11
            CovDtildebeta_UL_12 = BETA_U_0*Gamtilde_ULL_102 + BETA_U_1*Gamtilde_ULL_112 + BETA_U_2* &
                                  Gamtilde_ULL_122 + dDBETA_UL_12
            CovDtildebeta_UL_20 = BETA_U_0*Gamtilde_ULL_200 + BETA_U_1*Gamtilde_ULL_201 + BETA_U_2* &
                                  Gamtilde_ULL_202 + dDBETA_UL_20
            CovDtildebeta_UL_21 = BETA_U_0*Gamtilde_ULL_201 + BETA_U_1*Gamtilde_ULL_211 + BETA_U_2* &
                                  Gamtilde_ULL_212 + dDBETA_UL_21
            CovDtildebeta_UL_22 = BETA_U_0*Gamtilde_ULL_202 + BETA_U_1*Gamtilde_ULL_212 + BETA_U_2* &
                                  Gamtilde_ULL_222 + dDBETA_UL_22

            divTildebeta = CovDtildebeta_UL_00 + CovDtildebeta_UL_11 + CovDtildebeta_UL_22

            KOSigma = 0.100000000000000d0

            rhs(CHI_RHS, i, j, k) = AdvDBETACHI + (2.0d0/3.0d0)*CHI*(ALPHA*(KHAT + 2*THETAFUNC) - &
                                                                     divTildebeta) + KOSigma*dKODCHI

            rhs(GAMTILDE_LL_00_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_00 + AdvDBETAGAMTILDE_LL_00 + (4.0d0/3.0d0)* &
                                               GAMTILDE_LL_00*dDBETA_UL_00 - 2.0d0/3.0d0*GAMTILDE_LL_00* &
                                               dDBETA_UL_11 - 2.0d0/3.0d0*GAMTILDE_LL_00*dDBETA_UL_22 + 2* &
                                               GAMTILDE_LL_01*dDBETA_UL_10 + 2*GAMTILDE_LL_02*dDBETA_UL_20 + &
                                               KOSigma*dKODGAMTILDE_LL_00

            rhs(GAMTILDE_LL_01_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_01 + AdvDBETAGAMTILDE_LL_01 + GAMTILDE_LL_00* &
                                               dDBETA_UL_01 + (1.0d0/3.0d0)*GAMTILDE_LL_01*dDBETA_UL_00 + (1.0d0 &
                                                                 /3.0d0)*GAMTILDE_LL_01*dDBETA_UL_11 - 2.0d0/3.0d0*GAMTILDE_LL_01* &
                                               dDBETA_UL_22 + GAMTILDE_LL_02*dDBETA_UL_21 + GAMTILDE_LL_11* &
                                               dDBETA_UL_10 + GAMTILDE_LL_12*dDBETA_UL_20 + KOSigma* &
                                               dKODGAMTILDE_LL_01
            rhs(GAMTILDE_LL_02_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_02 + AdvDBETAGAMTILDE_LL_02 + GAMTILDE_LL_00* &
                                               dDBETA_UL_02 + GAMTILDE_LL_01*dDBETA_UL_12 + (1.0d0/3.0d0)* &
                                               GAMTILDE_LL_02*dDBETA_UL_00 - 2.0d0/3.0d0*GAMTILDE_LL_02* &
                                               dDBETA_UL_11 + (1.0d0/3.0d0)*GAMTILDE_LL_02*dDBETA_UL_22 + &
                                               GAMTILDE_LL_12*dDBETA_UL_10 + GAMTILDE_LL_22*dDBETA_UL_20 + &
                                               KOSigma*dKODGAMTILDE_LL_02
            rhs(GAMTILDE_LL_11_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_11 + AdvDBETAGAMTILDE_LL_11 + 2*GAMTILDE_LL_01* &
                                               dDBETA_UL_01 - 2.0d0/3.0d0*GAMTILDE_LL_11*dDBETA_UL_00 + (4.0d0/ &
                                                                  3.0d0)*GAMTILDE_LL_11*dDBETA_UL_11 - 2.0d0/3.0d0*GAMTILDE_LL_11* &
                                               dDBETA_UL_22 + 2*GAMTILDE_LL_12*dDBETA_UL_21 + KOSigma* &
                                               dKODGAMTILDE_LL_11
            rhs(GAMTILDE_LL_12_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_12 + AdvDBETAGAMTILDE_LL_12 + GAMTILDE_LL_01* &
                                               dDBETA_UL_02 + GAMTILDE_LL_02*dDBETA_UL_01 + GAMTILDE_LL_11* &
                                               dDBETA_UL_12 - 2.0d0/3.0d0*GAMTILDE_LL_12*dDBETA_UL_00 + (1.0d0/ &
                                                                 3.0d0)*GAMTILDE_LL_12*dDBETA_UL_11 + (1.0d0/3.0d0)*GAMTILDE_LL_12 &
                                               *dDBETA_UL_22 + GAMTILDE_LL_22*dDBETA_UL_21 + KOSigma* &
                                               dKODGAMTILDE_LL_12
            rhs(GAMTILDE_LL_22_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_22 + AdvDBETAGAMTILDE_LL_22 + 2*GAMTILDE_LL_02* &
                                               dDBETA_UL_02 + 2*GAMTILDE_LL_12*dDBETA_UL_12 - 2.0d0/3.0d0* &
                                               GAMTILDE_LL_22*dDBETA_UL_00 - 2.0d0/3.0d0*GAMTILDE_LL_22* &
                                               dDBETA_UL_11 + (4.0d0/3.0d0)*GAMTILDE_LL_22*dDBETA_UL_22 + &
                                               KOSigma*dKODGAMTILDE_LL_22

            kappa1 = 1.0d0

            kappa2 = 0

            rhs(KHAT_RHS, i, j, k) = ALPHA*ATILDE_LL_00*Atilde_UU_00 + 2*ALPHA*ATILDE_LL_01* &
                                     Atilde_UU_01 + 2*ALPHA*ATILDE_LL_02*Atilde_UU_02 + ALPHA* &
                                     ATILDE_LL_11*Atilde_UU_11 + 2*ALPHA*ATILDE_LL_12*Atilde_UU_12 + &
                                     ALPHA*ATILDE_LL_22*Atilde_UU_22 + ALPHA*THETAFUNC*kappa1*(1 - &
                                                                           kappa2) + (1.0d0/3.0d0)*ALPHA*(KHAT + 2*THETAFUNC)**2 + &
                                     AdvDBETAKHAT - CovDDalpha_LL_00*invgam_UU_00 - 2*CovDDalpha_LL_01 &
                                     *invgam_UU_01 - 2*CovDDalpha_LL_02*invgam_UU_02 - &
                                     CovDDalpha_LL_11*invgam_UU_11 - 2*CovDDalpha_LL_12*invgam_UU_12 - &
                                     CovDDalpha_LL_22*invgam_UU_22 + KOSigma*dKODKHAT

            rhs(ATILDE_LL_00_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_00*Atilde_UL_00 + ALPHA*ATILDE_LL_00*(KHAT + 2* &
                                                                         THETAFUNC) - 2*ALPHA*ATILDE_LL_01*Atilde_UL_10 - 2*ALPHA* &
                                             ATILDE_LL_02*Atilde_UL_20 + (4.0d0/3.0d0)*ATILDE_LL_00* &
                                             dDBETA_UL_00 - 2.0d0/3.0d0*ATILDE_LL_00*dDBETA_UL_11 - 2.0d0/ &
                                             3.0d0*ATILDE_LL_00*dDBETA_UL_22 + 2*ATILDE_LL_01*dDBETA_UL_10 + 2 &
                                             *ATILDE_LL_02*dDBETA_UL_20 + AdvDBETAATILDE_LL_00 + CHI*(ALPHA* &
                                                                          RTF_LL_00 - CovDDalphaTF_LL_00) + KOSigma*dKODATILDE_LL_00
            rhs(ATILDE_LL_01_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_00*Atilde_UL_01 - 2*ALPHA*ATILDE_LL_01* &
                                             Atilde_UL_11 + ALPHA*ATILDE_LL_01*(KHAT + 2*THETAFUNC) - 2*ALPHA* &
                                             ATILDE_LL_02*Atilde_UL_21 + ATILDE_LL_00*dDBETA_UL_01 + (1.0d0/ &
                                                                    3.0d0)*ATILDE_LL_01*dDBETA_UL_00 + (1.0d0/3.0d0)*ATILDE_LL_01* &
                                             dDBETA_UL_11 - 2.0d0/3.0d0*ATILDE_LL_01*dDBETA_UL_22 + &
                                             ATILDE_LL_02*dDBETA_UL_21 + ATILDE_LL_11*dDBETA_UL_10 + &
                                             ATILDE_LL_12*dDBETA_UL_20 + AdvDBETAATILDE_LL_01 + CHI*(ALPHA* &
                                                                          RTF_LL_01 - CovDDalphaTF_LL_01) + KOSigma*dKODATILDE_LL_01
            rhs(ATILDE_LL_02_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_00*Atilde_UL_02 - 2*ALPHA*ATILDE_LL_01* &
                                             Atilde_UL_12 - 2*ALPHA*ATILDE_LL_02*Atilde_UL_22 + ALPHA* &
                                             ATILDE_LL_02*(KHAT + 2*THETAFUNC) + ATILDE_LL_00*dDBETA_UL_02 + &
                                             ATILDE_LL_01*dDBETA_UL_12 + (1.0d0/3.0d0)*ATILDE_LL_02* &
                                             dDBETA_UL_00 - 2.0d0/3.0d0*ATILDE_LL_02*dDBETA_UL_11 + (1.0d0/ &
                                                                    3.0d0)*ATILDE_LL_02*dDBETA_UL_22 + ATILDE_LL_12*dDBETA_UL_10 + &
                                             ATILDE_LL_22*dDBETA_UL_20 + AdvDBETAATILDE_LL_02 + CHI*(ALPHA* &
                                                                          RTF_LL_02 - CovDDalphaTF_LL_02) + KOSigma*dKODATILDE_LL_02
            rhs(ATILDE_LL_11_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_01*Atilde_UL_01 - 2*ALPHA*ATILDE_LL_11* &
                                             Atilde_UL_11 + ALPHA*ATILDE_LL_11*(KHAT + 2*THETAFUNC) - 2*ALPHA* &
                                             ATILDE_LL_12*Atilde_UL_21 + 2*ATILDE_LL_01*dDBETA_UL_01 - 2.0d0/ &
                                             3.0d0*ATILDE_LL_11*dDBETA_UL_00 + (4.0d0/3.0d0)*ATILDE_LL_11* &
                                             dDBETA_UL_11 - 2.0d0/3.0d0*ATILDE_LL_11*dDBETA_UL_22 + 2* &
                                             ATILDE_LL_12*dDBETA_UL_21 + AdvDBETAATILDE_LL_11 + CHI*(ALPHA* &
                                                                          RTF_LL_11 - CovDDalphaTF_LL_11) + KOSigma*dKODATILDE_LL_11
            rhs(ATILDE_LL_12_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_01*Atilde_UL_02 - 2*ALPHA*ATILDE_LL_11* &
                                             Atilde_UL_12 - 2*ALPHA*ATILDE_LL_12*Atilde_UL_22 + ALPHA* &
                                             ATILDE_LL_12*(KHAT + 2*THETAFUNC) + ATILDE_LL_01*dDBETA_UL_02 + &
                                             ATILDE_LL_02*dDBETA_UL_01 + ATILDE_LL_11*dDBETA_UL_12 - 2.0d0/ &
                                             3.0d0*ATILDE_LL_12*dDBETA_UL_00 + (1.0d0/3.0d0)*ATILDE_LL_12* &
                                             dDBETA_UL_11 + (1.0d0/3.0d0)*ATILDE_LL_12*dDBETA_UL_22 + &
                                             ATILDE_LL_22*dDBETA_UL_21 + AdvDBETAATILDE_LL_12 + CHI*(ALPHA* &
                                                                          RTF_LL_12 - CovDDalphaTF_LL_12) + KOSigma*dKODATILDE_LL_12
            rhs(ATILDE_LL_22_RHS, i, j, k) = -2*ALPHA*ATILDE_LL_02*Atilde_UL_02 - 2*ALPHA*ATILDE_LL_12* &
                                             Atilde_UL_12 - 2*ALPHA*ATILDE_LL_22*Atilde_UL_22 + ALPHA* &
                                             ATILDE_LL_22*(KHAT + 2*THETAFUNC) + 2*ATILDE_LL_02*dDBETA_UL_02 + &
                                             2*ATILDE_LL_12*dDBETA_UL_12 - 2.0d0/3.0d0*ATILDE_LL_22* &
                                             dDBETA_UL_00 - 2.0d0/3.0d0*ATILDE_LL_22*dDBETA_UL_11 + (4.0d0/ &
                                                                   3.0d0)*ATILDE_LL_22*dDBETA_UL_22 + AdvDBETAATILDE_LL_22 + CHI*( &
                                             ALPHA*RTF_LL_22 - CovDDalphaTF_LL_22) + KOSigma*dKODATILDE_LL_22

            rhs(THETAFUNC_RHS, i, j, k) = -1.0d0/2.0d0*ALPHA*ATILDE_LL_00*Atilde_UU_00 - ALPHA*ATILDE_LL_01* &
                                          Atilde_UU_01 - ALPHA*ATILDE_LL_02*Atilde_UU_02 - 1.0d0/2.0d0* &
                                          ALPHA*ATILDE_LL_11*Atilde_UU_11 - ALPHA*ATILDE_LL_12*Atilde_UU_12 &
                                          - 1.0d0/2.0d0*ALPHA*ATILDE_LL_22*Atilde_UU_22 + 0.5d0* &
                                          ALPHA*(Rsclr - 2*THETAFUNC*kappa1*(kappa2 + 2) + (2.0d0/3.0d0)*( &
                                                 KHAT + 2*THETAFUNC)**2) + AdvDBETATHETAFUNC + KOSigma* &
                                          dKODTHETAFUNC

            rhs(GAMTILDE_U_0_RHS, i, j, k) = 2*ALPHA*Atilde_UU_00*Gamtilde_ULL_000 + 4*ALPHA*Atilde_UU_01* &
                                             Gamtilde_ULL_001 + 4*ALPHA*Atilde_UU_02*Gamtilde_ULL_002 + 2* &
                                             ALPHA*Atilde_UU_11*Gamtilde_ULL_011 + 4*ALPHA*Atilde_UU_12* &
                                             Gamtilde_ULL_012 + 2*ALPHA*Atilde_UU_22*Gamtilde_ULL_022 - 2* &
                                             ALPHA*kappa1*(GAMTILDE_U_0 - GamtildeD_U_0) + 2*ALPHA*(-3.0d0/ &
                                                                 2.0d0*Atilde_UU_00*dDCHI_L_0/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_00* &
                                                                        (2*dDKHAT_L_0 + dDTHETAFUNC_L_0)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_01*dDCHI_L_1/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_01*(2* &
                                                                           dDKHAT_L_1 + dDTHETAFUNC_L_1)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_02*dDCHI_L_2/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_02*(2* &
                                                                        dDKHAT_L_2 + dDTHETAFUNC_L_2)) + AdvDBETAGAMTILDE_U_0 - 2* &
                                             Atilde_UU_00*dDALPHA_L_0 - 2*Atilde_UU_01*dDALPHA_L_1 - 2* &
                                             Atilde_UU_02*dDALPHA_L_2 - 1.0d0/3.0d0*GamtildeD_U_0*dDBETA_UL_00 &
                                             + (2.0d0/3.0d0)*GamtildeD_U_0*dDBETA_UL_11 + (2.0d0/3.0d0)* &
                                             GamtildeD_U_0*dDBETA_UL_22 - GamtildeD_U_1*dDBETA_UL_01 - &
                                             GamtildeD_U_2*dDBETA_UL_02 + KOSigma*dKODGAMTILDE_U_0 + (4.0d0/ &
                                                                         3.0d0)*dDDBETA_ULL_000*invGAMTILDE_UU_00 + (7.0d0/3.0d0)* &
                                             dDDBETA_ULL_001*invGAMTILDE_UU_01 + (7.0d0/3.0d0)*dDDBETA_ULL_002 &
                                             *invGAMTILDE_UU_02 + dDDBETA_ULL_011*invGAMTILDE_UU_11 + 2* &
                                             dDDBETA_ULL_012*invGAMTILDE_UU_12 + dDDBETA_ULL_022* &
                                             invGAMTILDE_UU_22 + (1.0d0/3.0d0)*dDDBETA_ULL_101* &
                                             invGAMTILDE_UU_00 + (1.0d0/3.0d0)*dDDBETA_ULL_111* &
                                             invGAMTILDE_UU_01 + (1.0d0/3.0d0)*dDDBETA_ULL_112* &
                                             invGAMTILDE_UU_02 + (1.0d0/3.0d0)*dDDBETA_ULL_202* &
                                             invGAMTILDE_UU_00 + (1.0d0/3.0d0)*dDDBETA_ULL_212* &
                                             invGAMTILDE_UU_01 + (1.0d0/3.0d0)*dDDBETA_ULL_222* &
                                             invGAMTILDE_UU_02
            rhs(GAMTILDE_U_1_RHS, i, j, k) = 2*ALPHA*Atilde_UU_00*Gamtilde_ULL_100 + 4*ALPHA*Atilde_UU_01* &
                                             Gamtilde_ULL_101 + 4*ALPHA*Atilde_UU_02*Gamtilde_ULL_102 + 2* &
                                             ALPHA*Atilde_UU_11*Gamtilde_ULL_111 + 4*ALPHA*Atilde_UU_12* &
                                             Gamtilde_ULL_112 + 2*ALPHA*Atilde_UU_22*Gamtilde_ULL_122 - 2* &
                                             ALPHA*kappa1*(GAMTILDE_U_1 - GamtildeD_U_1) + 2*ALPHA*(-3.0d0/ &
                                                                 2.0d0*Atilde_UU_01*dDCHI_L_0/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_01* &
                                                                        (2*dDKHAT_L_0 + dDTHETAFUNC_L_0)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_11*dDCHI_L_1/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_11*(2* &
                                                                           dDKHAT_L_1 + dDTHETAFUNC_L_1)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_12*dDCHI_L_2/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_12*(2* &
                                                                        dDKHAT_L_2 + dDTHETAFUNC_L_2)) + AdvDBETAGAMTILDE_U_1 - 2* &
                                             Atilde_UU_01*dDALPHA_L_0 - 2*Atilde_UU_11*dDALPHA_L_1 - 2* &
                                             Atilde_UU_12*dDALPHA_L_2 - GamtildeD_U_0*dDBETA_UL_10 + (2.0d0/ &
                                                                    3.0d0)*GamtildeD_U_1*dDBETA_UL_00 - 1.0d0/3.0d0*GamtildeD_U_1* &
                                             dDBETA_UL_11 + (2.0d0/3.0d0)*GamtildeD_U_1*dDBETA_UL_22 - &
                                             GamtildeD_U_2*dDBETA_UL_12 + KOSigma*dKODGAMTILDE_U_1 + (1.0d0/ &
                                                                         3.0d0)*dDDBETA_ULL_000*invGAMTILDE_UU_01 + (1.0d0/3.0d0)* &
                                             dDDBETA_ULL_001*invGAMTILDE_UU_11 + (1.0d0/3.0d0)*dDDBETA_ULL_002 &
                                             *invGAMTILDE_UU_12 + dDDBETA_ULL_100*invGAMTILDE_UU_00 + (7.0d0/ &
                                                                     3.0d0)*dDDBETA_ULL_101*invGAMTILDE_UU_01 + 2*dDDBETA_ULL_102* &
                                             invGAMTILDE_UU_02 + (4.0d0/3.0d0)*dDDBETA_ULL_111* &
                                             invGAMTILDE_UU_11 + (7.0d0/3.0d0)*dDDBETA_ULL_112* &
                                             invGAMTILDE_UU_12 + dDDBETA_ULL_122*invGAMTILDE_UU_22 + (1.0d0/ &
                                                                         3.0d0)*dDDBETA_ULL_202*invGAMTILDE_UU_01 + (1.0d0/3.0d0)* &
                                             dDDBETA_ULL_212*invGAMTILDE_UU_11 + (1.0d0/3.0d0)*dDDBETA_ULL_222 &
                                             *invGAMTILDE_UU_12
            rhs(GAMTILDE_U_2_RHS, i, j, k) = 2*ALPHA*Atilde_UU_00*Gamtilde_ULL_200 + 4*ALPHA*Atilde_UU_01* &
                                             Gamtilde_ULL_201 + 4*ALPHA*Atilde_UU_02*Gamtilde_ULL_202 + 2* &
                                             ALPHA*Atilde_UU_11*Gamtilde_ULL_211 + 4*ALPHA*Atilde_UU_12* &
                                             Gamtilde_ULL_212 + 2*ALPHA*Atilde_UU_22*Gamtilde_ULL_222 - 2* &
                                             ALPHA*kappa1*(GAMTILDE_U_2 - GamtildeD_U_2) + 2*ALPHA*(-3.0d0/ &
                                                                 2.0d0*Atilde_UU_02*dDCHI_L_0/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_02* &
                                                                        (2*dDKHAT_L_0 + dDTHETAFUNC_L_0)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_12*dDCHI_L_1/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_12*(2* &
                                                                           dDKHAT_L_1 + dDTHETAFUNC_L_1)) + 2*ALPHA*(-3.0d0/2.0d0* &
                                                                    Atilde_UU_22*dDCHI_L_2/CHI - 1.0d0/3.0d0*invGAMTILDE_UU_22*(2* &
                                                                        dDKHAT_L_2 + dDTHETAFUNC_L_2)) + AdvDBETAGAMTILDE_U_2 - 2* &
                                             Atilde_UU_02*dDALPHA_L_0 - 2*Atilde_UU_12*dDALPHA_L_1 - 2* &
                                             Atilde_UU_22*dDALPHA_L_2 - GamtildeD_U_0*dDBETA_UL_20 - &
                                             GamtildeD_U_1*dDBETA_UL_21 + (2.0d0/3.0d0)*GamtildeD_U_2* &
                                             dDBETA_UL_00 + (2.0d0/3.0d0)*GamtildeD_U_2*dDBETA_UL_11 - 1.0d0/ &
                                             3.0d0*GamtildeD_U_2*dDBETA_UL_22 + KOSigma*dKODGAMTILDE_U_2 + ( &
                                             1.0d0/3.0d0)*dDDBETA_ULL_000*invGAMTILDE_UU_02 + (1.0d0/3.0d0)* &
                                             dDDBETA_ULL_001*invGAMTILDE_UU_12 + (1.0d0/3.0d0)*dDDBETA_ULL_002 &
                                             *invGAMTILDE_UU_22 + (1.0d0/3.0d0)*dDDBETA_ULL_101* &
                                             invGAMTILDE_UU_02 + (1.0d0/3.0d0)*dDDBETA_ULL_111* &
                                             invGAMTILDE_UU_12 + (1.0d0/3.0d0)*dDDBETA_ULL_112* &
                                             invGAMTILDE_UU_22 + dDDBETA_ULL_200*invGAMTILDE_UU_00 + 2* &
                                             dDDBETA_ULL_201*invGAMTILDE_UU_01 + (7.0d0/3.0d0)*dDDBETA_ULL_202 &
                                             *invGAMTILDE_UU_02 + dDDBETA_ULL_211*invGAMTILDE_UU_11 + (7.0d0/ &
                                                                         3.0d0)*dDDBETA_ULL_212*invGAMTILDE_UU_12 + (4.0d0/3.0d0)* &
                                             dDDBETA_ULL_222*invGAMTILDE_UU_22

            mul = 1.0d0 !2.0d0/ALPHA !

            mus = ALPHA**(-2)

            rhs(ALPHA_RHS, i, j, k) = -ALPHA**2*KHAT*mul + AdvDBETAALPHA + KOSigma*dKODALPHA

            eta = 2.0d0

            rhs(BETA_U_0_RHS, i, j, k) = 0
            rhs(BETA_U_1_RHS, i, j, k) = 0
            rhs(BETA_U_2_RHS, i, j, k) = 0

            !rhs(BETA_U_0_RHS, i, j, k) =       ALPHA**2*GAMTILDE_U_0*mus + AdvDBETABETA_U_0 - BETA_U_0*eta + &
            !KOSigma*dKODBETA_U_0
            !rhs(BETA_U_1_RHS, i, j, k) =       ALPHA**2*GAMTILDE_U_1*mus + AdvDBETABETA_U_1 - BETA_U_1*eta + &
            !KOSigma*dKODBETA_U_1
            !rhs(BETA_U_2_RHS, i, j, k) =       ALPHA**2*GAMTILDE_U_2*mus + AdvDBETABETA_U_2 - BETA_U_2*eta + &
            !KOSigma*dKODBETA_U_2

         end do ! i
      end do ! j
   end do ! k

   call MoL_releaseDataPtr(tileDesc, rhs, activeRHS)
   call MoL_releaseDataPtr(tileDesc, vars, MOL_EVOLVED)

end subroutine Spacetime_molFastRHS_tile
