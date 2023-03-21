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
!! @brief ml_finalize implementation for IMEX

!> @ingroup MoLIMEX
!!
!! @brief Implements ml_finalize for IMEX
!!
!! @stubref{ml_finalize}
subroutine ml_finalize()
   use ml_imexData, only: ml_AE, ml_bE, ml_cE, ml_AI, ml_bI, ml_cI, FE, FI

   implicit none

   if (allocated(ml_AE)) deallocate (ml_AE)
   if (allocated(ml_bE)) deallocate (ml_bE)
   if (allocated(ml_cE)) deallocate (ml_cE)
   if (allocated(ml_AI)) deallocate (ml_AI)
   if (allocated(ml_bI)) deallocate (ml_bI)
   if (allocated(ml_cI)) deallocate (ml_cI)

   if (allocated(FE)) deallocate (FE)
   if (allocated(FI)) deallocate (FI)
end subroutine ml_finalize
