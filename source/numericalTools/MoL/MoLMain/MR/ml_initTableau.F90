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
!! @brief ml_initTableau implementation for MR

!> @ingroup MoLMR
!!
!! @brief Implements ml_initTableau for MR
!!
!! @stubref{ml_initTableau}
subroutine ml_initTableau()
   use ml_mrData, only: ml_AF, ml_bF, ml_cF, ml_gamK, ml_wK, &
                        ml_gamBar, ml_wBar, ml_cS, ml_kmax, &
                        ml_nstages_slow, ml_nstages_fast, &
                        ml_slowOrder, ml_fastOrder, ml_slowMethod, ml_fastMethod
   use ml_interface, only: ml_error

   use erk, only: euler_init, rk2_heun_init, rk3_ssp_init, rk4_init
   use gark, only: gark3_init, gark4_init

   implicit none

   integer :: k

   if (allocated(ml_AF)) deallocate (ml_AF)
   if (allocated(ml_bF)) deallocate (ml_bF)
   if (allocated(ml_cF)) deallocate (ml_cF)
   if (allocated(ml_gamK)) deallocate (ml_gamK)
   if (allocated(ml_wK)) deallocate (ml_wK)
   if (allocated(ml_gamBar)) deallocate (ml_gamBar)
   if (allocated(ml_wBar)) deallocate (ml_wBar)
   if (allocated(ml_cS)) deallocate (ml_cS)

   select case (ml_fastMethod)

   case ("erk-euler")
      call euler_init(ml_AF, ml_bF, ml_cF, ml_fastOrder, ml_nstages_fast)

   case ("erk-rk2-heun")
      call rk2_heun_init(ml_AF, ml_bF, ml_cF, ml_fastOrder, ml_nstages_fast)

   case ("erk-rk3-ssp")
      call rk3_ssp_init(ml_AF, ml_bF, ml_cF, ml_fastOrder, ml_nstages_fast)

   case ("erk-rk4")
      call rk4_init(ml_AF, ml_bF, ml_cF, ml_fastOrder, ml_nstages_fast)

   case default
      call ml_error("Unkown ERK method")
   end select

   select case (ml_slowMethod)
   case ("mr-gark3")
      call gark3_init(ml_kmax, ml_gamK, ml_wK, ml_cS, ml_slowOrder, ml_nstages_slow)

   case ("mr-gark4")
      call gark4_init(ml_kmax, ml_gamK, ml_wK, ml_cS, ml_slowOrder, ml_nstages_slow)
   end select

   allocate (ml_gamBar(ml_nstages_slow, ml_nstages_slow))
   allocate (ml_wBar(ml_nstages_slow, ml_nstages_slow))

   ml_gamBar = 0.0
   ml_wBar = 0.0

   do k = 1, ml_kmax
      ml_gamBar = ml_gamBar + ml_gamK(:, :, k)/k
      ml_wBar = ml_wBar + ml_wK(:, :, k)/k
   end do ! k

end subroutine ml_initTableau
