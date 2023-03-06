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
!! @brief ml_initTableau implementation for IMEX

!> @ingroup MoLIMEX
!!
!! @brief Implements ml_initTableau for IMEX
!!
!! @stubref{ml_initTableau}
subroutine ml_initTableau()
   use ml_imexData, only: ml_method, ml_order, ml_nstages, ml_AE, ml_bE, ml_cE, ml_AI, ml_bI, ml_cI

   use ml_interface, only: ml_error

   use imex_tableau, only: fbe_init, ssp2_222_init, ssp2_322_init, ssp2_332_init, ssp3_332_init, ssp3_433_init, &
                           ark_111_init, ark_121_init, ark_122_init, ark_222_init, ark_232_init, ark_233_init, &
                           ark_343_init, ark_443_init

   implicit none

   integer :: k

   if (allocated(ml_AE)) deallocate (ml_AE)
   if (allocated(ml_bE)) deallocate (ml_bE)
   if (allocated(ml_cE)) deallocate (ml_cE)
   if (allocated(ml_AI)) deallocate (ml_AI)
   if (allocated(ml_bI)) deallocate (ml_bI)
   if (allocated(ml_cI)) deallocate (ml_cI)

   select case (ml_method)

   case ("imex-fbe")
      call fbe_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ssp2-222")
      call ssp2_222_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ssp2-322")
      call ssp2_322_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ssp2-332")
      call ssp2_332_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ssp3-332")
      call ssp3_332_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ssp3-433")
      call ssp3_433_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-111")
      call ark_111_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-121")
      call ark_121_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-122")
      call ark_122_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-222")
      call ark_222_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-232")
      call ark_232_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-233")
      call ark_233_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-343")
      call ark_343_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case ("imex-ark-443")
      call ark_443_init(ml_AI, ml_bI, ml_cI, ml_AE, ml_bE, ml_cE, ml_order, ml_nstages)

   case default
      call ml_error("Unkown IMEX method")
   end select

end subroutine ml_initTableau
