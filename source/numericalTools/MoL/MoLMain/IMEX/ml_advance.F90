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
   use ml_imexData, only: ml_nstages, ml_AE, ml_bE, ml_cE, ml_AI, ml_bI, ml_cI, FE, FI
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast, ml_implicitUpdate
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcs(2*ml_nstages + 1)
   real    :: facs(2*ml_nstages + 1)

   integer :: s

   srcs(1) = MOL_INITIAL
   facs(1) = 1.0

   do s = 1, ml_nstages
      srcs(2*s) = FE(s)
      srcs(2*s + 1) = FI(s)
   end do ! s

   do s = 1, ml_nstages
      if ((any(ml_AE(s, 1:s - 1) .ne. 0.0)) .or. (any(ml_AI(s, 1:s - 1) .ne. 0.0))) then
         facs(2::2) = ml_AE(s, :)*dt
         facs(3::2) = ml_AI(s, :)*dt

         call ml_memAddToVars(MOL_EVOLVED, 0.0, 2*s - 1, srcs(1:2*s - 1), facs(1:2*s - 1))

         call ml_postUpdate(t + ml_cE(s)*dt)
         call ml_postUpdateFast(t + ml_cE(s)*dt)
      end if

      if (ml_AI(s, s) .ne. 0.0) then
         call ml_implicitUpdate(t + ml_cI(s)*dt, ml_AI(s, s)*dt)
      end if

      call ml_calcRHS(MOL_RHS_EXPLICIT, FE(s), t + ml_cE(s)*dt, ml_bE(s)*dt)
      call ml_calcRHS(MOL_RHS_IMPLICIT, FI(s), t + ml_cI(s)*dt, ml_bI(s)*dt)
   end do ! s

   ! Final linear combination
   if (any(ml_bE .ne. 0.0) .or. any(ml_bI .ne. 0.0)) then
      do s = 1, ml_nstages
         facs(2*s) = ml_bE(s)*dt
         facs(2*s + 1) = ml_bI(s)*dt
      end do ! s

      call ml_memAddToVars(MOL_EVOLVED, 0.0, 2*ml_nstages + 1, srcs, facs)

      call ml_postUpdate(t + dt)
      call ml_postUpdateFast(t + dt)
   end if
end subroutine ml_advance
