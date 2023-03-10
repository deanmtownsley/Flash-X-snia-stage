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
!! @brief ml_finalize implementation for MR

!> @ingroup MoLMR
!!
!! @brief Implements ml_finalize for MR
!!
!! @stubref{ml_finalize}
subroutine ml_finalize()
   use ml_mrData, only: FE, FI, FF, ml_gamK, ml_wK, &
                        ml_gamBar, ml_wBar, ml_cS, ml_AF, ml_bF, ml_cF

   implicit none

   if (allocated(FE)) deallocate (FE)
   if (allocated(FI)) deallocate (FI)
   if (allocated(FF)) deallocate (FF)

   if (allocated(ml_gamK)) deallocate (ml_gamK)
   if (allocated(ml_wK)) deallocate (ml_wK)

   if (allocated(ml_gamBar)) deallocate (ml_gamBar)
   if (allocated(ml_wBar)) deallocate (ml_wBar)

   if (allocated(ml_cS)) deallocate (ml_cS)

   if (allocated(ml_AF)) deallocate (ml_AF)
   if (allocated(ml_bF)) deallocate (ml_bF)
   if (allocated(ml_cF)) deallocate (ml_cF)
end subroutine ml_finalize
