subroutine z4c_calculateConstraintViolation(vars, lim, del)

#include "Z4c.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   integer, dimension(LOW:HIGH, MDIM) :: lim
   real, dimension(MDIM) :: del

   integer :: i, j, k
   real ::  dx0, dx1, dx2, idx0, idx1, idx2

   real :: CHI, KHAT, THETAFUNC
   real :: GAMTILDE_LL_00, GAMTILDE_LL_01, GAMTILDE_LL_02, GAMTILDE_LL_11, GAMTILDE_LL_12, GAMTILDE_LL_22
   real :: ATILDE_LL_00, ATILDE_LL_01, ATILDE_LL_02, ATILDE_LL_11, ATILDE_LL_12, ATILDE_LL_22
   real :: GAMTILDE_U_0, GAMTILDE_U_1, GAMTILDE_U_2

   real :: dDCHI_L_0, dDCHI_L_1, dDCHI_L_2
   real :: dDDCHI_LL_00, dDDCHI_LL_01, dDDCHI_LL_02, &
           dDDCHI_LL_11, dDDCHI_LL_12, dDDCHI_LL_22

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

   real :: dDKHAT_L_0, dDKHAT_L_1, dDKHAT_L_2

   real :: dDATILDE_LLL_000, dDATILDE_LLL_001, dDATILDE_LLL_002, &
           dDATILDE_LLL_010, dDATILDE_LLL_011, dDATILDE_LLL_012, &
           dDATILDE_LLL_020, dDATILDE_LLL_021, dDATILDE_LLL_022, &
           dDATILDE_LLL_110, dDATILDE_LLL_111, dDATILDE_LLL_112, &
           dDATILDE_LLL_120, dDATILDE_LLL_121, dDATILDE_LLL_122, &
           dDATILDE_LLL_220, dDATILDE_LLL_221, dDATILDE_LLL_222

   real :: dDTHETAFUNC_L_0, dDTHETAFUNC_L_1, dDTHETAFUNC_L_2

   real :: dDGAMTILDE_UL_00, dDGAMTILDE_UL_01, dDGAMTILDE_UL_02, &
           dDGAMTILDE_UL_10, dDGAMTILDE_UL_11, dDGAMTILDE_UL_12, &
           dDGAMTILDE_UL_20, dDGAMTILDE_UL_21, dDGAMTILDE_UL_22

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

   real :: Gam_ULL_000, Gam_ULL_001, Gam_ULL_002, &
           Gam_ULL_011, Gam_ULL_012, Gam_ULL_022, &
           Gam_ULL_100, Gam_ULL_101, Gam_ULL_102, &
           Gam_ULL_111, Gam_ULL_112, Gam_ULL_122, &
           Gam_ULL_200, Gam_ULL_201, Gam_ULL_202, &
           Gam_ULL_211, Gam_ULL_212, Gam_ULL_222

   real :: gam_LL_00, gam_LL_01, gam_LL_02, &
           gam_LL_11, gam_LL_12, gam_LL_22

   real :: invgam_UU_00, invgam_UU_01, invgam_UU_02, &
           invgam_UU_11, invgam_UU_12, invgam_UU_22

   real :: GamtildeD_U_0, GamtildeD_U_1, GamtildeD_U_2

   real :: Rtilde_LL_00, Rtilde_LL_01, Rtilde_LL_02, &
           Rtilde_LL_11, Rtilde_LL_12, Rtilde_LL_22

   real :: CovDtildechi_L_0, CovDtildechi_L_1, CovDtildechi_L_2

   real :: CovDtildeDtildechi_LL_00, CovDtildeDtildechi_LL_01, CovDtildeDtildechi_LL_02, &
           CovDtildeDtildechi_LL_11, CovDtildeDtildechi_LL_12, CovDtildeDtildechi_LL_22

   real :: Rchi_LL_00, Rchi_LL_01, Rchi_LL_02, Rchi_LL_11, Rchi_LL_12, Rchi_LL_22
   real :: R_LL_00, R_LL_01, R_LL_02, R_LL_11, R_LL_12, R_LL_22

   real :: K_LL_00, K_LL_01, K_LL_02, K_LL_11, K_LL_12, K_LL_22
   ! real :: K_UU_00, K_UU_01, K_UU_02, K_UU_11, K_UU_12, K_UU_22

   ! real :: K_LU_00, K_LU_01, K_LU_02, &
   !         K_LU_10, K_LU_11, K_LU_12, &
   !         K_LU_20, K_LU_21, K_LU_22

   real :: CovDK_LLL_000, CovDK_LLL_001, CovDK_LLL_002, &
           CovDK_LLL_010, CovDK_LLL_011, CovDK_LLL_012, &
           CovDK_LLL_020, CovDK_LLL_021, CovDK_LLL_022, &
           CovDK_LLL_110, CovDK_LLL_111, CovDK_LLL_112, &
           CovDK_LLL_120, CovDK_LLL_121, CovDK_LLL_122, &
           CovDK_LLL_220, CovDK_LLL_221, CovDK_LLL_222

   real :: CovDK_L_0, CovDK_L_1, CovDK_L_2

   real :: M_L_0, M_L_1, M_L_2
   real :: Z_L_0, Z_L_1, Z_L_2

   real :: Ksclr
   real :: dDKsclr_L_0, dDKsclr_L_1, dDKsclr_L_2
   real :: invCHI, Rsclr, K2
   real :: H, M2, Z2, C2

   ! For finite differences (stencils given as dimensions)
   real, dimension(-2:2), parameter :: dDi = [1d0/12d0, -2d0/3d0, 0d0, 2d0/3d0, -1d0/12d0]
   real, dimension(-2:2), parameter :: dDDii = [-1d0/12d0, 4d0/3d0, -5d0/2d0, 4d0/3d0, -1d0/12d0]
   real, dimension(-2:2, -2:2), parameter :: dDDij = reshape([1d0/144d0, -1d0/18d0, 0d0, 1d0/18d0, -1d0/144d0, &
                                                              -1d0/18d0, 4d0/9d0, 0d0, -4d0/9d0, 1d0/18d0, &
                                                              0d0, 0d0, 0d0, 0d0, 0d0, &
                                                              1d0/18d0, -4d0/9d0, 0d0, 4d0/9d0, -1d0/18d0, &
                                                              -1d0/144d0, 1d0/18d0, 0d0, -1d0/18d0, 1d0/144d0], [5, 5])

   dx0 = del(IAXIS)
   dx1 = del(JAXIS)
   dx2 = del(KAXIS)

   idx0 = 1d0/dx0
   idx1 = 1d0/dx1
   idx2 = 1d0/dx2

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

            dDCHI_L_0 = sum(dDi*vars(CHI_VAR, i - 2:i + 2, j, k))*idx0
            dDCHI_L_1 = sum(dDi*vars(CHI_VAR, i, j - 2:j + 2, k))*idx1
            dDCHI_L_2 = sum(dDi*vars(CHI_VAR, i, j, k - 2:k + 2))*idx2

            dDDCHI_LL_00 = sum(dDDii*vars(CHI_VAR, i - 2:i + 2, j, k))*idx0**2
            dDDCHI_LL_01 = sum(dDDij*vars(CHI_VAR, i - 2:i + 2, j - 2:j + 2, k))*idx0*idx1
            dDDCHI_LL_02 = sum(dDDij*vars(CHI_VAR, i - 2:i + 2, j, k - 2:k + 2))*idx0*idx2
            dDDCHI_LL_11 = sum(dDDii*vars(CHI_VAR, i, j - 2:j + 2, k))*idx1**2
            dDDCHI_LL_12 = sum(dDDij*vars(CHI_VAR, i, j - 2:j + 2, k - 2:k + 2))*idx1*idx2
            dDDCHI_LL_22 = sum(dDDii*vars(CHI_VAR, i, j, k - 2:k + 2))*idx2**2

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

            dDGAMTILDE_UL_00 = sum(dDi*vars(GAMTILDE_U_0_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_UL_10 = sum(dDi*vars(GAMTILDE_U_1_VAR, i - 2:i + 2, j, k))*idx0
            dDGAMTILDE_UL_20 = sum(dDi*vars(GAMTILDE_U_2_VAR, i - 2:i + 2, j, k))*idx0

            dDGAMTILDE_UL_01 = sum(dDi*vars(GAMTILDE_U_0_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_UL_11 = sum(dDi*vars(GAMTILDE_U_1_VAR, i, j - 2:j + 2, k))*idx1
            dDGAMTILDE_UL_21 = sum(dDi*vars(GAMTILDE_U_2_VAR, i, j - 2:j + 2, k))*idx1

            dDGAMTILDE_UL_02 = sum(dDi*vars(GAMTILDE_U_0_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_UL_12 = sum(dDi*vars(GAMTILDE_U_1_VAR, i, j, k - 2:k + 2))*idx2
            dDGAMTILDE_UL_22 = sum(dDi*vars(GAMTILDE_U_2_VAR, i, j, k - 2:k + 2))*idx2

            dDKHAT_L_0 = sum(dDi*vars(KHAT_VAR, i - 2:i + 2, j, k))*idx0
            dDKHAT_L_1 = sum(dDi*vars(KHAT_VAR, i, j - 2:j + 2, k))*idx1
            dDKHAT_L_2 = sum(dDi*vars(KHAT_VAR, i, j, k - 2:k + 2))*idx2

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

            dDTHETAFUNC_L_0 = sum(dDi*vars(THETAFUNC_VAR, i - 2:i + 2, j, k))*idx0
            dDTHETAFUNC_L_1 = sum(dDi*vars(THETAFUNC_VAR, i, j - 2:j + 2, k))*idx1
            dDTHETAFUNC_L_2 = sum(dDi*vars(THETAFUNC_VAR, i, j, k - 2:k + 2))*idx2

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
                           0.5d0*(dDDGAMTILDE_LLLL_0000*invGAMTILDE_UU_00 + &
                                  dDDGAMTILDE_LLLL_0011*invGAMTILDE_UU_11 + &
                                  dDDGAMTILDE_LLLL_0022*invGAMTILDE_UU_22)

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
                           Gamtilde_LLL_212*Gamtilde_ULL_202*invGAMTILDE_UU_22 - &
                           0.5d0*dDDGAMTILDE_LLLL_0100*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_0101*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_0102*invGAMTILDE_UU_02 - 0.5d0*dDDGAMTILDE_LLLL_0111*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_0112*invGAMTILDE_UU_12 - 0.5d0*dDDGAMTILDE_LLLL_0122*invGAMTILDE_UU_22 + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_000*Gamtilde_ULL_001 + Gamtilde_LLL_100*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_001*Gamtilde_ULL_101 + Gamtilde_LLL_101*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_002*Gamtilde_ULL_201 + Gamtilde_LLL_102*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_000*Gamtilde_ULL_011 + Gamtilde_LLL_100*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_001 + Gamtilde_LLL_101*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_111 + Gamtilde_LLL_101*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_002*Gamtilde_ULL_211 + Gamtilde_LLL_102*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_011*Gamtilde_ULL_101 + Gamtilde_LLL_111*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_012*Gamtilde_ULL_201 + Gamtilde_LLL_112*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_000*Gamtilde_ULL_012 + Gamtilde_LLL_100*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_001*Gamtilde_ULL_112 + Gamtilde_LLL_101*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_001 + Gamtilde_LLL_102*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_212 + Gamtilde_LLL_102*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_012*Gamtilde_ULL_101 + Gamtilde_LLL_112*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_022*Gamtilde_ULL_201 + Gamtilde_LLL_122*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_001*Gamtilde_ULL_011 + Gamtilde_LLL_101*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_011*Gamtilde_ULL_111 + Gamtilde_LLL_111*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_012*Gamtilde_ULL_211 + Gamtilde_LLL_112*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_001*Gamtilde_ULL_012 + Gamtilde_LLL_101*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_002*Gamtilde_ULL_011 + Gamtilde_LLL_102*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_011*Gamtilde_ULL_112 + Gamtilde_LLL_111*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_012*Gamtilde_ULL_111 + Gamtilde_LLL_112*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_012*Gamtilde_ULL_212 + Gamtilde_LLL_112*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_022*Gamtilde_ULL_211 + Gamtilde_LLL_122*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_002*Gamtilde_ULL_012 + Gamtilde_LLL_102*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_012*Gamtilde_ULL_112 + Gamtilde_LLL_112*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_022*Gamtilde_ULL_212 + Gamtilde_LLL_122*Gamtilde_ULL_202)
            Rtilde_LL_02 = 0.5d0*(GAMTILDE_LL_00*dDGAMTILDE_UL_02 + GAMTILDE_LL_01*dDGAMTILDE_UL_12 + &
                                  GAMTILDE_LL_02*dDGAMTILDE_UL_00 + GAMTILDE_LL_02*dDGAMTILDE_UL_22 + &
                                  GAMTILDE_LL_12*dDGAMTILDE_UL_10 + GAMTILDE_LL_22*dDGAMTILDE_UL_20 + &
                                  GamtildeD_U_0*Gamtilde_LLL_002 + GamtildeD_U_0*Gamtilde_LLL_200 + &
                                  GamtildeD_U_1*Gamtilde_LLL_012 + GamtildeD_U_1*Gamtilde_LLL_201 + &
                                  GamtildeD_U_2*Gamtilde_LLL_022 + GamtildeD_U_2*Gamtilde_LLL_202) + &
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
                           Gamtilde_LLL_222*Gamtilde_ULL_202*invGAMTILDE_UU_22 - &
                           0.5d0*dDDGAMTILDE_LLLL_0200*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_0201*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_0202*invGAMTILDE_UU_02 - 0.5d0*dDDGAMTILDE_LLLL_0211*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_0212*invGAMTILDE_UU_12 - 0.5d0*dDDGAMTILDE_LLLL_0222*invGAMTILDE_UU_22 + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_000*Gamtilde_ULL_002 + Gamtilde_LLL_200*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_001*Gamtilde_ULL_102 + Gamtilde_LLL_201*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_002*Gamtilde_ULL_202 + Gamtilde_LLL_202*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_000*Gamtilde_ULL_012 + Gamtilde_LLL_200*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_002 + Gamtilde_LLL_201*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_001*Gamtilde_ULL_112 + Gamtilde_LLL_201*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_002*Gamtilde_ULL_212 + Gamtilde_LLL_202*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_011*Gamtilde_ULL_102 + Gamtilde_LLL_211*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_012*Gamtilde_ULL_202 + Gamtilde_LLL_212*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_000*Gamtilde_ULL_022 + Gamtilde_LLL_200*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_001*Gamtilde_ULL_122 + Gamtilde_LLL_201*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_002 + Gamtilde_LLL_202*Gamtilde_ULL_000) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_002*Gamtilde_ULL_222 + Gamtilde_LLL_202*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_012*Gamtilde_ULL_102 + Gamtilde_LLL_212*Gamtilde_ULL_100) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_022*Gamtilde_ULL_202 + Gamtilde_LLL_222*Gamtilde_ULL_200) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_001*Gamtilde_ULL_012 + Gamtilde_LLL_201*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_011*Gamtilde_ULL_112 + Gamtilde_LLL_211*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_012*Gamtilde_ULL_212 + Gamtilde_LLL_212*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_001*Gamtilde_ULL_022 + Gamtilde_LLL_201*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_002*Gamtilde_ULL_012 + Gamtilde_LLL_202*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_011*Gamtilde_ULL_122 + Gamtilde_LLL_211*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_012*Gamtilde_ULL_112 + Gamtilde_LLL_212*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_012*Gamtilde_ULL_222 + Gamtilde_LLL_212*Gamtilde_ULL_202) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_022*Gamtilde_ULL_212 + Gamtilde_LLL_222*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_002*Gamtilde_ULL_022 + Gamtilde_LLL_202*Gamtilde_ULL_002) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_012*Gamtilde_ULL_122 + Gamtilde_LLL_212*Gamtilde_ULL_102) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_022*Gamtilde_ULL_222 + Gamtilde_LLL_222*Gamtilde_ULL_202)
            Rtilde_LL_11 = GAMTILDE_LL_01*dDGAMTILDE_UL_01 + GAMTILDE_LL_11*dDGAMTILDE_UL_11 + &
                           GAMTILDE_LL_12*dDGAMTILDE_UL_21 + GamtildeD_U_0*Gamtilde_LLL_101 + &
                           GamtildeD_U_1*Gamtilde_LLL_111 + GamtildeD_U_2*Gamtilde_LLL_112 + &
                           Gamtilde_LLL_001*Gamtilde_ULL_001*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_001*Gamtilde_ULL_011*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_001*Gamtilde_ULL_012*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_011*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_011*Gamtilde_ULL_012*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_001*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_011*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_012*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_100*Gamtilde_ULL_001*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_100*Gamtilde_ULL_011*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_100*Gamtilde_ULL_012*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_101*Gamtilde_ULL_001*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_101*Gamtilde_ULL_011*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_101*Gamtilde_ULL_012*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_101*Gamtilde_ULL_101*invGAMTILDE_UU_00 + &
                           3*Gamtilde_LLL_101*Gamtilde_ULL_111*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_101*Gamtilde_ULL_112*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_001*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_011*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_012*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_201*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_211*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_102*Gamtilde_ULL_212*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_111*Gamtilde_ULL_101*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_111*Gamtilde_ULL_111*invGAMTILDE_UU_11 + &
                           3*Gamtilde_LLL_111*Gamtilde_ULL_112*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_112*Gamtilde_ULL_101*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_112*Gamtilde_ULL_111*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_112*Gamtilde_ULL_112*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_112*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_112*Gamtilde_ULL_211*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_112*Gamtilde_ULL_212*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_122*Gamtilde_ULL_201*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_122*Gamtilde_ULL_211*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_122*Gamtilde_ULL_212*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_201*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_211*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_201*Gamtilde_ULL_212*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_201*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_211*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_211*Gamtilde_ULL_212*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_201*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_211*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_212*Gamtilde_ULL_212*invGAMTILDE_UU_22 - &
                           0.5d0*dDDGAMTILDE_LLLL_1100*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_1101*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_1102*invGAMTILDE_UU_02 - 0.5d0*dDDGAMTILDE_LLLL_1111*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_1112*invGAMTILDE_UU_12 - 0.5d0*dDDGAMTILDE_LLLL_1122*invGAMTILDE_UU_22
            Rtilde_LL_12 = 0.5d0*(GAMTILDE_LL_01*dDGAMTILDE_UL_02 + GAMTILDE_LL_02*dDGAMTILDE_UL_01 + &
                                  GAMTILDE_LL_11*dDGAMTILDE_UL_12 + GAMTILDE_LL_12*dDGAMTILDE_UL_11 + &
                                  GAMTILDE_LL_12*dDGAMTILDE_UL_22 + GAMTILDE_LL_22*dDGAMTILDE_UL_21 + &
                                  GamtildeD_U_0*Gamtilde_LLL_102 + GamtildeD_U_0*Gamtilde_LLL_201 + &
                                  GamtildeD_U_1*Gamtilde_LLL_112 + GamtildeD_U_1*Gamtilde_LLL_211 + &
                                  GamtildeD_U_2*Gamtilde_LLL_122 + GamtildeD_U_2*Gamtilde_LLL_212) + &
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
                           Gamtilde_LLL_222*Gamtilde_ULL_212*invGAMTILDE_UU_22 - &
                           0.5d0*dDDGAMTILDE_LLLL_1200*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_1201*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_1202*invGAMTILDE_UU_02 - 0.5d0*dDDGAMTILDE_LLLL_1211*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_1212*invGAMTILDE_UU_12 - 0.5d0*dDDGAMTILDE_LLLL_1222*invGAMTILDE_UU_22 + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_100*Gamtilde_ULL_002 + Gamtilde_LLL_200*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_101*Gamtilde_ULL_102 + Gamtilde_LLL_201*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_00*(Gamtilde_LLL_102*Gamtilde_ULL_202 + Gamtilde_LLL_202*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_100*Gamtilde_ULL_012 + Gamtilde_LLL_200*Gamtilde_ULL_011) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_101*Gamtilde_ULL_002 + Gamtilde_LLL_201*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_101*Gamtilde_ULL_112 + Gamtilde_LLL_201*Gamtilde_ULL_111) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_102*Gamtilde_ULL_212 + Gamtilde_LLL_202*Gamtilde_ULL_211) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_111*Gamtilde_ULL_102 + Gamtilde_LLL_211*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_01*(Gamtilde_LLL_112*Gamtilde_ULL_202 + Gamtilde_LLL_212*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_100*Gamtilde_ULL_022 + Gamtilde_LLL_200*Gamtilde_ULL_012) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_101*Gamtilde_ULL_122 + Gamtilde_LLL_201*Gamtilde_ULL_112) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_102*Gamtilde_ULL_002 + Gamtilde_LLL_202*Gamtilde_ULL_001) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_102*Gamtilde_ULL_222 + Gamtilde_LLL_202*Gamtilde_ULL_212) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_112*Gamtilde_ULL_102 + Gamtilde_LLL_212*Gamtilde_ULL_101) + &
                           invGAMTILDE_UU_02*(Gamtilde_LLL_122*Gamtilde_ULL_202 + Gamtilde_LLL_222*Gamtilde_ULL_201) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_101*Gamtilde_ULL_012 + Gamtilde_LLL_201*Gamtilde_ULL_011) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_111*Gamtilde_ULL_112 + Gamtilde_LLL_211*Gamtilde_ULL_111) + &
                           invGAMTILDE_UU_11*(Gamtilde_LLL_112*Gamtilde_ULL_212 + Gamtilde_LLL_212*Gamtilde_ULL_211) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_101*Gamtilde_ULL_022 + Gamtilde_LLL_201*Gamtilde_ULL_012) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_102*Gamtilde_ULL_012 + Gamtilde_LLL_202*Gamtilde_ULL_011) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_111*Gamtilde_ULL_122 + Gamtilde_LLL_211*Gamtilde_ULL_112) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_112*Gamtilde_ULL_112 + Gamtilde_LLL_212*Gamtilde_ULL_111) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_112*Gamtilde_ULL_222 + Gamtilde_LLL_212*Gamtilde_ULL_212) + &
                           invGAMTILDE_UU_12*(Gamtilde_LLL_122*Gamtilde_ULL_212 + Gamtilde_LLL_222*Gamtilde_ULL_211) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_102*Gamtilde_ULL_022 + Gamtilde_LLL_202*Gamtilde_ULL_012) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_112*Gamtilde_ULL_122 + Gamtilde_LLL_212*Gamtilde_ULL_112) + &
                           invGAMTILDE_UU_22*(Gamtilde_LLL_122*Gamtilde_ULL_222 + Gamtilde_LLL_222*Gamtilde_ULL_212)
            Rtilde_LL_22 = GAMTILDE_LL_02*dDGAMTILDE_UL_02 + GAMTILDE_LL_12*dDGAMTILDE_UL_12 + &
                           GAMTILDE_LL_22*dDGAMTILDE_UL_22 + GamtildeD_U_0*Gamtilde_LLL_202 + &
                           GamtildeD_U_1*Gamtilde_LLL_212 + GamtildeD_U_2*Gamtilde_LLL_222 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_002*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_012*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_002*Gamtilde_ULL_022*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_002*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_012*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_012*Gamtilde_ULL_022*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_002*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_012*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_022*Gamtilde_ULL_022*invGAMTILDE_UU_22 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_102*invGAMTILDE_UU_00 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_112*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_102*Gamtilde_ULL_122*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_102*invGAMTILDE_UU_01 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_112*invGAMTILDE_UU_11 + &
                           Gamtilde_LLL_112*Gamtilde_ULL_122*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_112*invGAMTILDE_UU_12 + &
                           Gamtilde_LLL_122*Gamtilde_ULL_122*invGAMTILDE_UU_22 + &
                           2*Gamtilde_LLL_200*Gamtilde_ULL_002*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_200*Gamtilde_ULL_012*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_200*Gamtilde_ULL_022*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_002*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_012*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_022*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_102*invGAMTILDE_UU_00 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_112*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_201*Gamtilde_ULL_122*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_202*Gamtilde_ULL_002*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_202*Gamtilde_ULL_012*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_202*Gamtilde_ULL_022*invGAMTILDE_UU_22 + &
                           3*Gamtilde_LLL_202*Gamtilde_ULL_202*invGAMTILDE_UU_00 + &
                           3*Gamtilde_LLL_202*Gamtilde_ULL_212*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_202*Gamtilde_ULL_222*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_211*Gamtilde_ULL_102*invGAMTILDE_UU_01 + &
                           2*Gamtilde_LLL_211*Gamtilde_ULL_112*invGAMTILDE_UU_11 + &
                           2*Gamtilde_LLL_211*Gamtilde_ULL_122*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_212*Gamtilde_ULL_102*invGAMTILDE_UU_02 + &
                           2*Gamtilde_LLL_212*Gamtilde_ULL_112*invGAMTILDE_UU_12 + &
                           2*Gamtilde_LLL_212*Gamtilde_ULL_122*invGAMTILDE_UU_22 + &
                           3*Gamtilde_LLL_212*Gamtilde_ULL_202*invGAMTILDE_UU_01 + &
                           3*Gamtilde_LLL_212*Gamtilde_ULL_212*invGAMTILDE_UU_11 + &
                           3*Gamtilde_LLL_212*Gamtilde_ULL_222*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_222*Gamtilde_ULL_202*invGAMTILDE_UU_02 + &
                           3*Gamtilde_LLL_222*Gamtilde_ULL_212*invGAMTILDE_UU_12 + &
                           3*Gamtilde_LLL_222*Gamtilde_ULL_222*invGAMTILDE_UU_22 - &
                           0.5d0*dDDGAMTILDE_LLLL_2200*invGAMTILDE_UU_00 - dDDGAMTILDE_LLLL_2201*invGAMTILDE_UU_01 - &
                           dDDGAMTILDE_LLLL_2202*invGAMTILDE_UU_02 - 0.5d0*dDDGAMTILDE_LLLL_2211*invGAMTILDE_UU_11 - &
                           dDDGAMTILDE_LLLL_2212*invGAMTILDE_UU_12 - 0.5d0*dDDGAMTILDE_LLLL_2222*invGAMTILDE_UU_22

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

            Rchi_LL_00 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_00*invGAMTILDE_UU_00*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_00*invGAMTILDE_UU_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_00*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_00*invGAMTILDE_UU_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_00*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_00*invGAMTILDE_UU_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_00*invGAMTILDE_UU_00*invCHI**2 - &
                         0.25d0*CovDtildechi_L_0**2*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_00*invGAMTILDE_UU_01*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_00*invGAMTILDE_UU_02*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_00*invGAMTILDE_UU_11*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_00*invGAMTILDE_UU_12*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_00*invGAMTILDE_UU_22*invCHI**2
            Rchi_LL_01 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_01*invGAMTILDE_UU_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_01*invGAMTILDE_UU_01*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_01*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_01*invGAMTILDE_UU_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_01*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_01*invGAMTILDE_UU_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_01*invGAMTILDE_UU_00*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_01*invGAMTILDE_UU_01*invCHI**2 - &
                         0.25d0*CovDtildechi_L_0*CovDtildechi_L_1*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_01*invGAMTILDE_UU_02*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_01*invGAMTILDE_UU_11*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_01*invGAMTILDE_UU_12*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_01*invGAMTILDE_UU_22*invCHI**2
            Rchi_LL_02 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_02*invGAMTILDE_UU_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_02*invGAMTILDE_UU_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_02*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_02*invGAMTILDE_UU_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_02*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_02*invGAMTILDE_UU_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_02*invGAMTILDE_UU_00*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_02*invGAMTILDE_UU_01*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_02*invGAMTILDE_UU_02*invCHI**2 - &
                         0.25d0*CovDtildechi_L_0*CovDtildechi_L_2*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_02*invGAMTILDE_UU_11*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_02*invGAMTILDE_UU_12*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_02*invGAMTILDE_UU_22*invCHI**2
            Rchi_LL_11 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_11*invGAMTILDE_UU_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_11*invGAMTILDE_UU_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_11*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_11*invGAMTILDE_UU_11*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_11*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_11*invGAMTILDE_UU_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_11*invGAMTILDE_UU_00*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_11*invGAMTILDE_UU_01*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_11*invGAMTILDE_UU_02*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_11*invGAMTILDE_UU_11*invCHI**2 - &
                         0.25d0*CovDtildechi_L_1**2*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_11*invGAMTILDE_UU_12*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_11*invGAMTILDE_UU_22*invCHI**2
            Rchi_LL_12 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_12*invGAMTILDE_UU_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_12*invGAMTILDE_UU_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_12*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_12*invGAMTILDE_UU_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_12*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_12*invGAMTILDE_UU_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_12*invGAMTILDE_UU_00*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_12*invGAMTILDE_UU_01*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_12*invGAMTILDE_UU_02*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_12*invGAMTILDE_UU_11*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_12*invGAMTILDE_UU_12*invCHI**2 - &
                         0.25d0*CovDtildechi_L_1*CovDtildechi_L_2*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_12*invGAMTILDE_UU_22*invCHI**2
            Rchi_LL_22 = 0.5d0*CovDtildeDtildechi_LL_00*GAMTILDE_LL_22*invGAMTILDE_UU_00*invCHI + &
                         CovDtildeDtildechi_LL_01*GAMTILDE_LL_22*invGAMTILDE_UU_01*invCHI + &
                         CovDtildeDtildechi_LL_02*GAMTILDE_LL_22*invGAMTILDE_UU_02*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_11*GAMTILDE_LL_22*invGAMTILDE_UU_11*invCHI + &
                         CovDtildeDtildechi_LL_12*GAMTILDE_LL_22*invGAMTILDE_UU_12*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*GAMTILDE_LL_22*invGAMTILDE_UU_22*invCHI + &
                         0.5d0*CovDtildeDtildechi_LL_22*invCHI - &
                         0.75d0*CovDtildechi_L_0**2*GAMTILDE_LL_22*invGAMTILDE_UU_00*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_1*GAMTILDE_LL_22*invGAMTILDE_UU_01*invCHI**2 - &
                         1.5d0*CovDtildechi_L_0*CovDtildechi_L_2*GAMTILDE_LL_22*invGAMTILDE_UU_02*invCHI**2 - &
                         0.75d0*CovDtildechi_L_1**2*GAMTILDE_LL_22*invGAMTILDE_UU_11*invCHI**2 - &
                         1.5d0*CovDtildechi_L_1*CovDtildechi_L_2*GAMTILDE_LL_22*invGAMTILDE_UU_12*invCHI**2 - &
                         0.75d0*CovDtildechi_L_2**2*GAMTILDE_LL_22*invGAMTILDE_UU_22*invCHI**2 - &
                         0.25d0*CovDtildechi_L_2**2*invCHI**2

            R_LL_00 = Rchi_LL_00 + Rtilde_LL_00
            R_LL_01 = Rchi_LL_01 + Rtilde_LL_01
            R_LL_02 = Rchi_LL_02 + Rtilde_LL_02
            R_LL_11 = Rchi_LL_11 + Rtilde_LL_11
            R_LL_12 = Rchi_LL_12 + Rtilde_LL_12
            R_LL_22 = Rchi_LL_22 + Rtilde_LL_22

            Rsclr = R_LL_00*invgam_UU_00 + 2*R_LL_01*invgam_UU_01 + 2*R_LL_02*invgam_UU_02 + &
                    R_LL_11*invgam_UU_11 + 2*R_LL_12*invgam_UU_12 + R_LL_22*invgam_UU_22

            Ksclr = KHAT + 2d0*THETAFUNC

            K_LL_00 = ATILDE_LL_00*invCHI + gam_LL_00*Ksclr/3d0
            K_LL_01 = ATILDE_LL_01*invCHI + gam_LL_01*Ksclr/3d0
            K_LL_02 = ATILDE_LL_02*invCHI + gam_LL_02*Ksclr/3d0
            K_LL_11 = ATILDE_LL_11*invCHI + gam_LL_11*Ksclr/3d0
            K_LL_12 = ATILDE_LL_12*invCHI + gam_LL_12*Ksclr/3d0
            K_LL_22 = ATILDE_LL_22*invCHI + gam_LL_22*Ksclr/3d0

            ! K_LU_00 = K_LL_00*invgam_UU_00 + K_LL_01*invgam_UU_01 + K_LL_02*invgam_UU_02
            ! K_LU_01 = K_LL_00*invgam_UU_01 + K_LL_01*invgam_UU_11 + K_LL_02*invgam_UU_12
            ! K_LU_02 = K_LL_00*invgam_UU_02 + K_LL_01*invgam_UU_12 + K_LL_02*invgam_UU_22

            ! K_LU_10 = K_LL_01*invgam_UU_00 + K_LL_11*invgam_UU_01 + K_LL_12*invgam_UU_02
            ! K_LU_11 = K_LL_01*invgam_UU_01 + K_LL_11*invgam_UU_11 + K_LL_12*invgam_UU_12
            ! K_LU_12 = K_LL_01*invgam_UU_02 + K_LL_11*invgam_UU_12 + K_LL_12*invgam_UU_22

            ! K_LU_20 = K_LL_02*invgam_UU_00 + K_LL_12*invgam_UU_01 + K_LL_22*invgam_UU_02
            ! K_LU_21 = K_LL_02*invgam_UU_01 + K_LL_12*invgam_UU_11 + K_LL_22*invgam_UU_12
            ! K_LU_22 = K_LL_02*invgam_UU_02 + K_LL_12*invgam_UU_12 + K_LL_22*invgam_UU_22

            ! K_UU_00 = invgam_UU_00*K_LU_00 + invgam_UU_01*K_LU_10 + invgam_UU_02*K_LU_20
            ! K_UU_01 = invgam_UU_00*K_LU_01 + invgam_UU_01*K_LU_11 + invgam_UU_02*K_LU_21
            ! K_UU_02 = invgam_UU_00*K_LU_02 + invgam_UU_01*K_LU_12 + invgam_UU_02*K_LU_22
            ! K_UU_11 = invgam_UU_01*K_LU_01 + invgam_UU_11*K_LU_11 + invgam_UU_12*K_LU_21
            ! K_UU_12 = invgam_UU_01*K_LU_02 + invgam_UU_11*K_LU_12 + invgam_UU_12*K_LU_22
            ! K_UU_22 = invgam_UU_02*K_LU_02 + invgam_UU_12*K_LU_12 + invgam_UU_22*K_LU_22

            ! K2 = K_UU_00*K_LL_00 + K_UU_11*K_LL_11 + K_UU_22 + K_LL_22 + &
            !      2d0*(K_UU_01*K_LL_01 + K_UU_02*K_LL_02 + K_UU_12 + K_LL_12)

            K2 = invgam_UU_00**2*K_LL_00**2 + 2*invgam_UU_02**2*K_LL_02**2 + &
                 2*invgam_UU_00*(2*invgam_UU_01*K_LL_00*K_LL_01 + invgam_UU_11*K_LL_01**2 + &
                                 K_LL_02*(2*invgam_UU_02*K_LL_00 + 2*invgam_UU_12*K_LL_01 + invgam_UU_22*K_LL_02)) + &
                 invgam_UU_11**2*K_LL_11**2 + 2*invgam_UU_01**2*(K_LL_01**2 + K_LL_00*K_LL_11) + &
                 4*(invgam_UU_02*invgam_UU_11*K_LL_01 + invgam_UU_02*invgam_UU_12*K_LL_02 + &
                    invgam_UU_11*invgam_UU_12*K_LL_11)*K_LL_12 + &
                 2*(invgam_UU_12**2 + invgam_UU_11*invgam_UU_22)*K_LL_12**2 + &
                 4*invgam_UU_01*(invgam_UU_02*K_LL_01*K_LL_02 + invgam_UU_11*K_LL_01*K_LL_11 + &
                                 invgam_UU_12*K_LL_02*K_LL_11 + invgam_UU_02*K_LL_00*K_LL_12 + &
                                 invgam_UU_12*K_LL_01*K_LL_12 + invgam_UU_22*K_LL_02*K_LL_12) + &
                 2*(invgam_UU_02**2*K_LL_00 + 2*invgam_UU_02*(invgam_UU_12*K_LL_01 + invgam_UU_22*K_LL_02) + &
                    invgam_UU_12*(invgam_UU_12*K_LL_11 + 2*invgam_UU_22*K_LL_12))*K_LL_22 + &
                 invgam_UU_22**2*K_LL_22**2

            dDKsclr_L_0 = dDKHAT_L_0 + 2d0*dDTHETAFUNC_L_0
            dDKsclr_L_1 = dDKHAT_L_1 + 2d0*dDTHETAFUNC_L_1
            dDKsclr_L_2 = dDKHAT_L_2 + 2d0*dDTHETAFUNC_L_2

            CovDK_LLL_000 = (dDATILDE_LLL_000 + (dDGAMTILDE_LLL_000*Ksclr + GAMTILDE_LL_00*dDKsclr_L_0)/3d0)*invCHI - &
                            2d0*K_LL_00*Gam_ULL_000 - &
                            2d0*K_LL_01*Gam_ULL_100 - &
                            2d0*K_LL_02*Gam_ULL_200 - &
                            K_LL_00*dDCHI_L_0*invCHI

            CovDK_LLL_001 = (dDATILDE_LLL_001 + (dDGAMTILDE_LLL_001*Ksclr + GAMTILDE_LL_00*dDKsclr_L_1)/3d0)*invCHI - &
                            2d0*K_LL_00*Gam_ULL_001 - &
                            2d0*K_LL_01*Gam_ULL_101 - &
                            2d0*K_LL_02*Gam_ULL_201 - &
                            K_LL_00*dDCHI_L_1*invCHI

            CovDK_LLL_002 = (dDATILDE_LLL_002 + (dDGAMTILDE_LLL_002*Ksclr + GAMTILDE_LL_00*dDKsclr_L_2)/3d0)*invCHI - &
                            2d0*K_LL_00*Gam_ULL_002 - &
                            2d0*K_LL_01*Gam_ULL_102 - &
                            2d0*K_LL_02*Gam_ULL_202 - &
                            K_LL_00*dDCHI_L_2*invCHI

            CovDK_LLL_010 = (dDATILDE_LLL_010 + (dDGAMTILDE_LLL_010*Ksclr + GAMTILDE_LL_01*dDKsclr_L_0)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_001 - &
                            K_LL_01*(Gam_ULL_101 + Gam_ULL_000) - &
                            K_LL_02*Gam_ULL_201 - &
                            K_LL_11*Gam_ULL_100 - &
                            K_LL_12*Gam_ULL_200 - &
                            K_LL_01*dDCHI_L_0*invCHI

            CovDK_LLL_011 = (dDATILDE_LLL_011 + (dDGAMTILDE_LLL_011*Ksclr + GAMTILDE_LL_01*dDKsclr_L_1)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_011 - &
                            K_LL_01*(Gam_ULL_111 + Gam_ULL_001) - &
                            K_LL_02*Gam_ULL_211 - &
                            K_LL_11*Gam_ULL_101 - &
                            K_LL_12*Gam_ULL_201 - &
                            K_LL_01*dDCHI_L_1*invCHI

            CovDK_LLL_012 = (dDATILDE_LLL_012 + (dDGAMTILDE_LLL_012*Ksclr + GAMTILDE_LL_01*dDKsclr_L_2)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_012 - &
                            K_LL_01*(Gam_ULL_112 + Gam_ULL_002) - &
                            K_LL_02*Gam_ULL_212 - &
                            K_LL_11*Gam_ULL_102 - &
                            K_LL_12*Gam_ULL_202 - &
                            K_LL_01*dDCHI_L_2*invCHI

            CovDK_LLL_020 = (dDATILDE_LLL_020 + (dDGAMTILDE_LLL_020*Ksclr + GAMTILDE_LL_02*dDKsclr_L_0)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_002 - &
                            K_LL_01*Gam_ULL_102 - &
                            K_LL_02*(Gam_ULL_202 + Gam_ULL_000) - &
                            K_LL_12*Gam_ULL_100 - &
                            K_LL_22*Gam_ULL_200 - &
                            K_LL_02*dDCHI_L_0*invCHI

            CovDK_LLL_021 = (dDATILDE_LLL_021 + (dDGAMTILDE_LLL_021*Ksclr + GAMTILDE_LL_02*dDKsclr_L_1)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_012 - &
                            K_LL_01*Gam_ULL_112 - &
                            K_LL_02*(Gam_ULL_212 + Gam_ULL_001) - &
                            K_LL_12*Gam_ULL_101 - &
                            K_LL_22*Gam_ULL_201 - &
                            K_LL_02*dDCHI_L_1*invCHI

            CovDK_LLL_022 = (dDATILDE_LLL_022 + (dDGAMTILDE_LLL_022*Ksclr + GAMTILDE_LL_02*dDKsclr_L_2)/3d0)*invCHI - &
                            K_LL_00*Gam_ULL_022 - &
                            K_LL_01*Gam_ULL_122 - &
                            K_LL_02*(Gam_ULL_222 + Gam_ULL_002) - &
                            K_LL_12*Gam_ULL_102 - &
                            K_LL_22*Gam_ULL_202 - &
                            K_LL_02*dDCHI_L_2*invCHI

            CovDK_LLL_110 = (dDATILDE_LLL_110 + (dDGAMTILDE_LLL_110*Ksclr + GAMTILDE_LL_11*dDKsclr_L_0)/3d0)*invCHI - &
                            2d0*K_LL_01*Gam_ULL_001 - &
                            2d0*K_LL_11*Gam_ULL_101 - &
                            2d0*K_LL_12*Gam_ULL_201 - &
                            K_LL_11*dDCHI_L_0*invCHI

            CovDK_LLL_111 = (dDATILDE_LLL_111 + (dDGAMTILDE_LLL_111*Ksclr + GAMTILDE_LL_11*dDKsclr_L_1)/3d0)*invCHI - &
                            2d0*K_LL_01*Gam_ULL_011 - &
                            2d0*K_LL_11*Gam_ULL_111 - &
                            2d0*K_LL_12*Gam_ULL_211 - &
                            K_LL_11*dDCHI_L_1*invCHI

            CovDK_LLL_112 = (dDATILDE_LLL_112 + (dDGAMTILDE_LLL_112*Ksclr + GAMTILDE_LL_11*dDKsclr_L_2)/3d0)*invCHI - &
                            2d0*K_LL_01*Gam_ULL_012 - &
                            2d0*K_LL_11*Gam_ULL_112 - &
                            2d0*K_LL_12*Gam_ULL_212 - &
                            K_LL_11*dDCHI_L_2*invCHI

            CovDK_LLL_120 = (dDATILDE_LLL_120 + (dDGAMTILDE_LLL_120*Ksclr + GAMTILDE_LL_12*dDKsclr_L_0)/3d0)*invCHI - &
                            K_LL_01*Gam_ULL_002 - &
                            K_LL_02*Gam_ULL_001 - &
                            K_LL_11*Gam_ULL_102 - &
                            K_LL_12*(Gam_ULL_202 + Gam_ULL_101) - &
                            K_LL_22*Gam_ULL_201 - &
                            K_LL_12*dDCHI_L_0*invCHI

            CovDK_LLL_121 = (dDATILDE_LLL_121 + (dDGAMTILDE_LLL_121*Ksclr + GAMTILDE_LL_12*dDKsclr_L_1)/3d0)*invCHI - &
                            K_LL_01*Gam_ULL_012 - &
                            K_LL_02*Gam_ULL_011 - &
                            K_LL_11*Gam_ULL_112 - &
                            K_LL_12*(Gam_ULL_212 + Gam_ULL_111) - &
                            K_LL_22*Gam_ULL_211 - &
                            K_LL_12*dDCHI_L_1*invCHI

            CovDK_LLL_122 = (dDATILDE_LLL_122 + (dDGAMTILDE_LLL_122*Ksclr + GAMTILDE_LL_12*dDKsclr_L_2)/3d0)*invCHI - &
                            K_LL_01*Gam_ULL_022 - &
                            K_LL_02*Gam_ULL_012 - &
                            K_LL_11*Gam_ULL_122 - &
                            K_LL_12*(Gam_ULL_222 + Gam_ULL_112) - &
                            K_LL_22*Gam_ULL_212 - &
                            K_LL_12*dDCHI_L_2*invCHI

            CovDK_LLL_220 = (dDATILDE_LLL_220 + (dDGAMTILDE_LLL_220*Ksclr + GAMTILDE_LL_22*dDKsclr_L_0)/3d0)*invCHI - &
                            2d0*K_LL_02*Gam_ULL_002 - &
                            2d0*K_LL_12*Gam_ULL_102 - &
                            2d0*K_LL_22*Gam_ULL_202 - &
                            K_LL_22*dDCHI_L_0*invCHI

            CovDK_LLL_221 = (dDATILDE_LLL_221 + (dDGAMTILDE_LLL_221*Ksclr + GAMTILDE_LL_22*dDKsclr_L_1)/3d0)*invCHI - &
                            2d0*K_LL_02*Gam_ULL_012 - &
                            2d0*K_LL_12*Gam_ULL_112 - &
                            2d0*K_LL_22*Gam_ULL_212 - &
                            K_LL_22*dDCHI_L_1*invCHI

            CovDK_LLL_222 = (dDATILDE_LLL_222 + (dDGAMTILDE_LLL_222*Ksclr + GAMTILDE_LL_22*dDKsclr_L_2)/3d0)*invCHI - &
                            2d0*K_LL_02*Gam_ULL_022 - &
                            2d0*K_LL_12*Gam_ULL_122 - &
                            2d0*K_LL_22*Gam_ULL_222 - &
                            K_LL_22*dDCHI_L_2*invCHI

            CovDK_L_0 = CovDK_LLL_000*invgam_UU_00 + CovDK_LLL_011*invgam_UU_11 + CovDK_LLL_022*invgam_UU_22 + &
                        (CovDK_LLL_001 + CovDK_LLL_010)*invgam_UU_01 + &
                        (CovDK_LLL_002 + CovDK_LLL_020)*invgam_UU_02 + &
                        (CovDK_LLL_012 + CovDK_LLL_021)*invgam_UU_12

            CovDK_L_1 = CovDK_LLL_010*invgam_UU_00 + CovDK_LLL_111*invgam_UU_11 + CovDK_LLL_122*invgam_UU_22 + &
                        (CovDK_LLL_011 + CovDK_LLL_110)*invgam_UU_01 + &
                        (CovDK_LLL_012 + CovDK_LLL_120)*invgam_UU_02 + &
                        (CovDK_LLL_112 + CovDK_LLL_121)*invgam_UU_12

            CovDK_L_2 = CovDK_LLL_020*invgam_UU_00 + CovDK_LLL_121*invgam_UU_11 + CovDK_LLL_222*invgam_UU_22 + &
                        (CovDK_LLL_021 + CovDK_LLL_120)*invgam_UU_01 + &
                        (CovDK_LLL_022 + CovDK_LLL_220)*invgam_UU_02 + &
                        (CovDK_LLL_122 + CovDK_LLL_221)*invgam_UU_12

            H = Rsclr - K2 + Ksclr**2 ! - 16d0*PI*vars(TE_VAR, i, j, k)

            M_L_0 = CovDK_L_0 - dDKsclr_L_0 ! - 8d0*PI*vars(TSX_VAR, i, j, k)
            M_L_1 = CovDK_L_1 - dDKsclr_L_1 ! - 8d0*PI*vars(TSY_VAR, i, j, k)
            M_L_2 = CovDK_L_2 - dDKsclr_L_2 ! - 8d0*PI*vars(TSZ_VAR, i, j, k)

            M2 = invgam_UU_00*M_L_0**2 + invgam_UU_11*M_L_1**2 + invgam_UU_22*M_L_2**2 + &
                 2d0*(invgam_UU_01*M_L_0*M_L_1 + invgam_UU_02*M_L_0*M_L_2 + invgam_UU_12*M_L_1*M_L_2)

            Z_L_0 = 0.5d0*(GAMTILDE_LL_00*GAMTILDE_U_0 + GAMTILDE_LL_01*GAMTILDE_U_1 + GAMTILDE_LL_02*GAMTILDE_U_2) - &
                    0.5d0*(dDGAMTILDE_LLL_000*invGAMTILDE_UU_00 + dDGAMTILDE_LLL_001*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_002*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_010*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_011*invGAMTILDE_UU_11 + dDGAMTILDE_LLL_012*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_020*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_021*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_022*invGAMTILDE_UU_22)

            Z_L_1 = 0.5d0*(GAMTILDE_LL_01*GAMTILDE_U_0 + GAMTILDE_LL_11*GAMTILDE_U_1 + GAMTILDE_LL_12*GAMTILDE_U_2) - &
                    0.5d0*(dDGAMTILDE_LLL_010*invGAMTILDE_UU_00 + dDGAMTILDE_LLL_011*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_012*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_110*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_111*invGAMTILDE_UU_11 + dDGAMTILDE_LLL_112*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_120*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_121*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_122*invGAMTILDE_UU_22)

            Z_L_2 = 0.5d0*(GAMTILDE_LL_02*GAMTILDE_U_0 + GAMTILDE_LL_12*GAMTILDE_U_1 + GAMTILDE_LL_22*GAMTILDE_U_2) - &
                    0.5d0*(dDGAMTILDE_LLL_020*invGAMTILDE_UU_00 + dDGAMTILDE_LLL_021*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_022*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_120*invGAMTILDE_UU_01 + &
                           dDGAMTILDE_LLL_121*invGAMTILDE_UU_11 + dDGAMTILDE_LLL_122*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_220*invGAMTILDE_UU_02 + dDGAMTILDE_LLL_221*invGAMTILDE_UU_12 + &
                           dDGAMTILDE_LLL_222*invGAMTILDE_UU_22)

            Z2 = invgam_UU_00*Z_L_0**2 + invgam_UU_11*Z_L_1**2 + invgam_UU_22*Z_L_2**2 + &
                 2d0*(invgam_UU_01*Z_L_0*Z_L_1 + invgam_UU_02*Z_L_0*Z_L_2 + invgam_UU_12*Z_L_1*Z_L_2)

            C2 = H**2 + M2 + 4d0*Z2 + THETAFUNC**2

            vars(Z4C_H_VAR, i, j, k) = H
            vars(Z4C_M_VAR, i, j, k) = sqrt(M2)
            vars(Z4C_Z_VAR, i, j, k) = sqrt(Z2)
            vars(Z4C_C_VAR, i, j, k) = sqrt(C2)
         end do ! 1
      end do ! j
   end do ! k
end subroutine z4c_calculateConstraintViolation
