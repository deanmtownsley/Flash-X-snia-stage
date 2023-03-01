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
!! @brief ml_advance implementation for ERK

!> @ingroup MoLERK
!!
!! @brief Implements ml_advance for ERK
!!
!! @stubref{ml_advance}
subroutine ml_advance(t, dt)
   use ml_erkData, only: ml_stages, ml_A, ml_b, ml_c, ml_K, ml_stages
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcs(1:ml_stages + 1)
   real    :: facs(1:ml_stages + 1)

   integer :: s

   srcs(1) = MOL_INITIAL
   srcs(2:) = ml_K

   facs(1) = 1.0

   do s = 1, ml_stages
      if (ml_c(s) .gt. 0.0) then
         facs(2:) = ml_A(s, :)*dt

         call ml_memAddToVars(MOL_EVOLVED, 0.0, s, srcs(:s), facs(:s))

         call ml_postUpdate(t + ml_c(s)*dt)
         call ml_postUpdateFast(t + ml_c(s)*dt)
      end if

      call ml_calcRHS(MOL_RHS_EXPLICIT, ml_K(s), t + ml_c(s)*dt, ml_b(s)*dt)
   end do ! s

   ! Final linear combination
   facs(2:) = ml_b*dt
   call ml_memAddToVars(MOL_EVOLVED, 0.0, ml_stages + 1, srcs, facs)

   call ml_postUpdate(t + dt)
   call ml_postUpdateFast(t + dt)
end subroutine ml_advance
