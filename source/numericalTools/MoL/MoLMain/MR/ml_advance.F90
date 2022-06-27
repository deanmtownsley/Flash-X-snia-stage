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
    use mr_data
    use MoL_functions
    use ml_interface, only: ml_calcRHS
    use ml_memInterface, only: ml_memAddToVars, ml_memCopy
    use gark, only: gamTau, wTau

#include "MoL.h"

    implicit none

    real, intent(in) :: t, dt

    integer :: srcsS(mr_nstages_slow)
    real    :: facsS(mr_nstages_slow)

    integer :: srcsF(mr_nstages_fast+1)
    real    :: facsF(mr_nstages_fast+1)

    integer :: sS, sF, n, j
    real :: t_stage, dc, theta, dtheta, t_fast, t_fast_stage, tau

    srcsF(1) = FAST_INITIAL
    srcsF(2:) = FF

    facsF(1) = 1d0

    dtheta = dt/mr_nsubcycle

    do sS = 1, mr_nstages_slow
        if (mod(sS,2) .ne. 0) then
            ! Slow stage
            t_stage = t + cS(sS)*dt

            if (cS(sS) .gt. 0d0) then
                do j = 1, sS-1, 2
                    srcsS(j)   = FE(j)
                    srcsS(j+1) = FI(j)

                    facsS(j)   = wBar(sS,j)*dt
                    facsS(j+1) = gamBar(sS,j)*dt
                end do
    
                call ml_memAddToVars(MOL_EVOLVED, 1d0, sS-1, srcsS(:sS-1), facsS(:sS-1))
    
                call MoL_postUpdate(t_stage)
    
                call MoL_implicitUpdate(t_stage, gamBar(sS,sS)*dt)
            end if

            call ml_calcRHS(MOL_RHS_EXPLICIT, FE(sS), t_stage)
            call ml_calcRHS(MOL_RHS_IMPLICIT, FI(sS), t_stage)
        else
            ! Fast step
            theta = 0d0

            t_stage = t + cS(sS-1)*dt
            t_fast = t_stage

            dc = cS(sS) - cS(sS-1)

            do n = 1, mr_nsubcycle
                call ml_memCopy(FAST_INITIAL, MOL_EVOLVED)

                do sF = 1, mr_nstages_fast
                    t_fast_stage = t_fast + dc*cF(sF)*dtheta
                    if (cF(sF) .gt. 0d0) then
                        facsF(2:sF) = AF(sF,:sF-1)*dtheta

                        call ml_memAddToVars(MOL_EVOLVED, 0d0, sF, srcsF(:sF), facsF(:sF))

                        call MoL_postUpdateFast(t_fast_stage)
                    end if

                    ! if (dc .gt. 0d0) call ml_calcRHS(MOL_RHS_FAST, FF(sF), t_fast_stage)
                    call ml_calcRHS(MOL_RHS_FAST, FF(sF), t_fast_stage)

                    tau = theta/dt
                    do j = 1, sS, 2
                        srcsS(j)   = FE(j)
                        srcsS(j+1) = FI(j)

                        facsS(j)   = wTau(sS,j,tau)
                        facsS(j+1) = gamTau(sS,j,tau)
                    end do

                    call ml_memAddToVars(FF(sF), dc, sS, srcsS(:sS), facsS(:sS))
                end do ! sF

                ! Final linear combination
                facsF(2:) = bF*dtheta
                call ml_memAddToVars(MOL_EVOLVED, 0d0, mr_nstages_fast+1, srcsF, facsF)

                theta = theta + dtheta
                t_fast = t_stage + dc*theta

                call MoL_postUpdateFast(t_fast)
            end do ! n
        end if
    end do ! sS

end subroutine ml_advance
