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
!! @brief ml_advance implementation for MR

!> @ingroup MoLMR
!!
!! @brief Implements ml_advance for MR
!!
!! @stubref{ml_advance}
subroutine ml_advance(t, dt)
   use ml_mrData, only: ml_nstages_slow, ml_nstages_fast, FAST_INITIAL, &
                        FF, FE, FI, ml_nsubcycle, ml_gamBar, ml_wBar, ml_cS, ml_cF, ml_bF, ml_AF, ml_gamK, ml_wK, ml_kmax
   use ml_functions, only: ml_implicitUpdate, ml_postUpdate, ml_postUpdateFast
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars, ml_memCopy
   use gark, only: gamTau, wTau

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcsS(ml_nstages_slow)
   real    :: facsS(ml_nstages_slow)

   integer :: srcsF(ml_nstages_fast + 1)
   real    :: facsF(ml_nstages_fast + 1)

   integer :: sS, sF, n, j
   real :: t_stage, dc, theta, dtheta, t_fast, t_fast_stage, tau

   srcsF(1) = FAST_INITIAL
   srcsF(2:) = FF

   facsF(1) = 1.0

   dtheta = dt/ml_nsubcycle

   do sS = 1, ml_nstages_slow
      if (mod(sS, 2) .ne. 0) then
         ! Slow stage

         ! The time that the RHS for this stage is evaluated at
         t_stage = t + ml_cS(sS)*dt

         ! If necessary, calculate the intermediate state for the RHS evaluation
         if (ml_cS(sS) .gt. 0.0) then

            ! Store source terms and scaling factors for the linear combination
            do j = 1, sS - 1, 2
               srcsS(j) = FE(j)
               srcsS(j + 1) = FI(j)

               facsS(j) = ml_wBar(sS, j)*dt
               facsS(j + 1) = ml_gamBar(sS, j)*dt
            end do

            ! U^j = U^n + dt*A^ji rhs_i
            call ml_memAddToVars(MOL_EVOLVED, 1.0, sS - 1, srcsS(:sS - 1), facsS(:sS - 1))

            ! Perform post-update work, e.g. con2prim, filling guard cells
            call ml_postUpdate(t_stage)

            ! Perform an necessary implicit update
            if (ml_gamBar(sS, sS) .gt. 0.0) then
               call ml_implicitUpdate(t_stage, ml_gamBar(sS, sS)*dt)
            end if
         end if

         ! Calculate the RHS terms for this stage
         call ml_calcRHS(MOL_RHS_EXPLICIT, FE(sS), t_stage, ml_wBar(ml_nstages_slow, sS)*dt)
         call ml_calcRHS(MOL_RHS_IMPLICIT, FI(sS), t_stage, ml_gamBar(ml_nstages_slow, sS)*dt)
      else
         ! Fast stage
         theta = 0.0

         ! The time that this stage starts at
         t_stage = t + ml_cS(sS - 1)*dt
         t_fast = t_stage

         ! Scaling factor from theta \in [0,dt] --> t \in [t^n,t^n+dt]
         dc = ml_cS(sS) - ml_cS(sS - 1)

         ! Integrate for some specified number of steps in the fast scheme
         do n = 1, ml_nsubcycle
            ! Save the initial step at the start of each fast step
            call ml_memCopy(FAST_INITIAL, MOL_EVOLVED)

            do sF = 1, ml_nstages_fast
               ! The time that this fast stage evaluation takes place at
               t_fast_stage = t_fast + dc*ml_cF(sF)*dtheta

               ! If necessary, compute the intermediate state for the RHS evaluation
               if (ml_cF(sF) .gt. 0.0) then
                  ! Scaling factors from the tableau
                  facsF(2:sF) = ml_AF(sF, :sF - 1)*dtheta

                  ! Linear combination of RHS terms for this stage
                  call ml_memAddToVars(MOL_EVOLVED, 0.0, sF, srcsF(:sF), facsF(:sF))

                  ! Guard-cell filling, etc.
                  call ml_postUpdateFast(t_fast_stage)
               end if

               if (dc .gt. 0.0) then
                  call ml_calcRHS(MOL_RHS_FAST, FF(sF), t_fast_stage, ml_bF(sF)*dtheta)
               end if

               ! Scaled time for forcing term
               tau = theta/dt

               ! Sources and scaling factors for the forcing term
               do j = 1, sS, 2
                  srcsS(j) = FE(j)
                  srcsS(j + 1) = FI(j)

                  facsS(j) = wTau(sS, j, tau, ml_wK, ml_kmax)
                  facsS(j + 1) = gamTau(sS, j, tau, ml_gamK, ml_kmax)
               end do

               call ml_memAddToVars(FF(sF), dc, sS, srcsS(:sS), facsS(:sS))
            end do ! sF

            ! Final linear combination
            facsF(2:) = ml_bF*dtheta
            call ml_memAddToVars(MOL_EVOLVED, 0.0, ml_nstages_fast + 1, srcsF, facsF)

            ! Update the time for the next fast (sub)step
            theta = theta + dtheta
            t_fast = t_stage + dc*theta

            call ml_postUpdateFast(t_fast)
         end do ! n
      end if
   end do ! sS

end subroutine ml_advance
