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
!! @brief ml_finalize implementation for ERK

!> @ingroup MoLERK
!!
!! @brief Implements ml_finalize for ERK
!!
!! @stubref{ml_finalize}
subroutine ml_finalize()
   use ml_erkData, only: ml_A, ml_b, ml_c, ml_K

   implicit none

   if (allocated(ml_A)) deallocate (ml_A)
   if (allocated(ml_b)) deallocate (ml_b)
   if (allocated(ml_c)) deallocate (ml_c)
   if (allocated(ml_K)) deallocate (ml_K)
end subroutine ml_finalize
