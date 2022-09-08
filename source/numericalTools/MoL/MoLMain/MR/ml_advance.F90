!!****if* source/numericalTools/MoL/MoLMain/MR/ml_advance
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!      ml_advance
!!
!!  SYNOPSIS
!!
!!      call ml_advance(real, intent(in) :: t
!!                      real, intent(in) :: dt)
!!
!!  DESCRIPTION
!!
!!      Take a timestep from t to t+dt
!!
!!  ARGUMENTS
!!
!!      t  : current time
!!      dt : size of timestep
!!
!!  TODO
!!
!!      When/if subcyling is available in Flash-X, this will extend
!!      to include a level-indicator as well
!!***
subroutine ml_advance(t, dt)
   use mr_data, only: mr_nstages_slow, mr_nstages_fast, FAST_INITIAL, &
                      FF, FE, FI, mr_nsubcycle, mr_gamBar, mr_wBar, mr_cS, mr_cF, mr_bF, mr_AF
   use MoL_functions
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars, ml_memCopy
   use gark, only: mr_gamTau, mr_wTau

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcsS(mr_nstages_slow)
   real    :: facsS(mr_nstages_slow)

   integer :: srcsF(mr_nstages_fast + 1)
   real    :: facsF(mr_nstages_fast + 1)

   integer :: sS, sF, n, j
   real :: t_stage, dc, theta, dtheta, t_fast, t_fast_stage, tau

   srcsF(1) = FAST_INITIAL
   srcsF(2:) = FF

   facsF(1) = 1d0

   dtheta = dt/mr_nsubcycle

   do sS = 1, mr_nstages_slow
      if (mod(sS, 2) .ne. 0) then
         ! Slow stage

         ! The time that the RHS for this stage is evaluated at
         t_stage = t + mr_cS(sS)*dt

         ! If necessary, calculate the intermediate state for the RHS evaluation
         if (mr_cS(sS) .gt. 0d0) then

            ! Store source terms and scaling factors for the linear combination
            do j = 1, sS - 1, 2
               srcsS(j) = FE(j)
               srcsS(j + 1) = FI(j)

               facsS(j) = mr_wBar(sS, j)*dt
               facsS(j + 1) = mr_gamBar(sS, j)*dt
            end do

            ! U^j = U^n + dt*A^ji rhs_i
            call ml_memAddToVars(MOL_EVOLVED, 1d0, sS - 1, srcsS(:sS - 1), facsS(:sS - 1))

            ! Perform post-update work, e.g. con2prim, filling guard cells
            call MoL_postUpdate(t_stage)

            ! Perform an necessary implicit update
            if (mr_gamBar(sS, sS) .gt. 0d0) then
               call MoL_implicitUpdate(t_stage, mr_gamBar(sS, sS)*dt)
            end if
         end if

         ! Calculate the RHS terms for this stage
         call ml_calcRHS(MOL_RHS_EXPLICIT, FE(sS), t_stage)
         call ml_calcRHS(MOL_RHS_IMPLICIT, FI(sS), t_stage)
      else
         ! Fast stage
         theta = 0d0

         ! The time that this stage starts at
         t_stage = t + mr_cS(sS - 1)*dt
         t_fast = t_stage

         ! Scaling factor from theta \in [0,dt] --> t \in [t^n,t^n+dt]
         dc = mr_cS(sS) - mr_cS(sS - 1)

         ! Integrate for some specified number of steps in the fast scheme
         do n = 1, mr_nsubcycle
            ! Save the initial step at the start of each fast step
            call ml_memCopy(FAST_INITIAL, MOL_EVOLVED)

            do sF = 1, mr_nstages_fast
               ! The time that this fast stage evaluation takes place at
               t_fast_stage = t_fast + dc*mr_cF(sF)*dtheta

               ! If necessary, compute the intermediate state for the RHS evaluation
               if (mr_cF(sF) .gt. 0d0) then
                  ! Scaling factors from the tableau
                  facsF(2:sF) = mr_AF(sF, :sF - 1)*dtheta

                  ! Linear combination of RHS terms for this stage
                  call ml_memAddToVars(MOL_EVOLVED, 0d0, sF, srcsF(:sF), facsF(:sF))

                  ! Guard-cell filling, etc.
                  call MoL_postUpdateFast(t_fast_stage)
               end if

               if (dc .gt. 0d0) then
                  call ml_calcRHS(MOL_RHS_FAST, FF(sF), t_fast_stage)
               end if

               ! Scaled time for forcing term
               tau = theta/dt

               ! Sources and scaling factors for the forcing term
               do j = 1, sS, 2
                  srcsS(j) = FE(j)
                  srcsS(j + 1) = FI(j)

                  facsS(j) = mr_wTau(sS, j, tau)
                  facsS(j + 1) = mr_gamTau(sS, j, tau)
               end do

               call ml_memAddToVars(FF(sF), dc, sS, srcsS(:sS), facsS(:sS))
            end do ! sF

            ! Final linear combination
            facsF(2:) = mr_bF*dtheta
            call ml_memAddToVars(MOL_EVOLVED, 0d0, mr_nstages_fast + 1, srcsF, facsF)

            ! Update the time for the next fast (sub)step
            theta = theta + dtheta
            t_fast = t_stage + dc*theta

            call MoL_postUpdateFast(t_fast)
         end do ! n
      end if
   end do ! sS

end subroutine ml_advance
