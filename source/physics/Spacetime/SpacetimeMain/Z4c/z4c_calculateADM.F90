subroutine z4c_calculateADM(vars, lim)

#include "Z4c.h"
#include "constants.h"

   implicit none

   real, dimension(:, :, :, :), pointer :: vars
   integer, dimension(LOW:HIGH, MDIM) :: lim

   integer :: i, j, k

   real :: invCHI, trK

   do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
      do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
         do i = lim(LOW, IAXIS), lim(HIGH, IAXIS)
            invCHI = 1d0/vars(CHI_VAR, i, j, k)
            trK = vars(KHAT_VAR, i, j, k) + 2d0*vars(THETAFUNC_VAR, i, j, k)

            vars(GXX_VAR, i, j, k) = vars(GAMTILDE_LL_00_VAR, i, j, k)*invCHI
            vars(GXY_VAR, i, j, k) = vars(GAMTILDE_LL_01_VAR, i, j, k)*invCHI
            vars(GXZ_VAR, i, j, k) = vars(GAMTILDE_LL_02_VAR, i, j, k)*invCHI
            vars(GYY_VAR, i, j, k) = vars(GAMTILDE_LL_11_VAR, i, j, k)*invCHI
            vars(GYZ_VAR, i, j, k) = vars(GAMTILDE_LL_12_VAR, i, j, k)*invCHI
            vars(GZZ_VAR, i, j, k) = vars(GAMTILDE_LL_22_VAR, i, j, k)*invCHI

            vars(KXX_VAR, i, j, k) = vars(ATILDE_LL_00_VAR, i, j, k)*invCHI + trK*vars(GXX_VAR, i, j, k)
            vars(KXY_VAR, i, j, k) = vars(ATILDE_LL_01_VAR, i, j, k)*invCHI + trK*vars(GXY_VAR, i, j, k)
            vars(KXZ_VAR, i, j, k) = vars(ATILDE_LL_02_VAR, i, j, k)*invCHI + trK*vars(GXZ_VAR, i, j, k)
            vars(KYY_VAR, i, j, k) = vars(ATILDE_LL_11_VAR, i, j, k)*invCHI + trK*vars(GYY_VAR, i, j, k)
            vars(KYZ_VAR, i, j, k) = vars(ATILDE_LL_12_VAR, i, j, k)*invCHI + trK*vars(GYZ_VAR, i, j, k)
            vars(KZZ_VAR, i, j, k) = vars(ATILDE_LL_22_VAR, i, j, k)*invCHI + trK*vars(GZZ_VAR, i, j, k)

            ! TODO: re-use lapse/shift variables for Z4c evolved
            vars(ALP_VAR, i, j, k) = vars(ALPHA_VAR, i, j, k)

            vars(BETAX_VAR, i, j, k) = vars(BETA_U_0_VAR, i, j, k)
            vars(BETAY_VAR, i, j, k) = vars(BETA_U_1_VAR, i, j, k)
            vars(BETAZ_VAR, i, j, k) = vars(BETA_U_2_VAR, i, j, k)
         end do ! j
      end do ! j
   end do ! k
end subroutine z4c_calculateADM
