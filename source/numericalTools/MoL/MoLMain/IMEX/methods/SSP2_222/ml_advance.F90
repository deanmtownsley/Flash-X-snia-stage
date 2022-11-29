!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief ml_advance implementation for IMEX-SSP2-222

!> @ingroup MoLIMEX
!!
!! @brief Implements ml_advance for IMEX-SSP2-222
!!
!! @stubref{ml_advance}
!!
!! @details
!! Implements the IMEX-SSP2(2,2,2) scheme from
!!
!!    Lorenzo Pareschi and Giovanni Russo,
!!    Implicit-explicit Runge-Kutta schmes and applications to hyperbolic systems with relaxation,
!!    Journal on Scientific Computing 25, 129-155 (2005)
!!    https://doi.org/10.1007/BF02728986
subroutine ml_advance(t, dt)
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast, ml_implicitUpdate
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer, parameter :: Y0 = MOL_INITIAL
   integer, parameter :: K0 = MOL_RHS + 0
   integer, parameter :: K1 = MOL_RHS + 1
   integer, parameter :: K2 = MOL_RHS + 2
   integer, parameter :: K3 = MOL_RHS + 3

   real, parameter :: gam = 1d0 - 1d0/sqrt(2d0)

   !! First stage !!

   ! Calculate and store explicit RHS at t (weight is 1/2)
   call ml_calcRHS(MOL_RHS_EXPLICIT, K0, t, 0.5d0*dt)

   ! Perform an implicit update of the implicit RHS at t + gam*dt
   call ml_implicitUpdate(t + gam*dt, gam*dt)

   ! Calculate and store implicit RHS at t+gam*dt (weight is 1/2)
   call ml_calcRHS(MOL_RHS_IMPLICIT, K1, t + gam*dt, 0.5d0*dt)

   !! Second stage !!

   ! Intermediate solution for use in implicit update
   call ml_memAddToVars(MOL_EVOLVED, 0d0, 2, [Y0, K1], [1d0, dt*(1d0 - 2d0*gam)])

   call ml_postUpdate(t + (1d0 - gam)*dt)
   call ml_postUpdateFast(t + (1d0 - gam)*dt)

   ! Perform an implicit update of the implicit RHS at t + (1-gam)*dt
   call ml_implicitUpdate(t + (1d0 - gam)*dt, gam*dt)

   ! Calculate and store implicit RHS at t+(1-gam)*dt (weight is 1/2)
   call ml_calcRHS(MOL_RHS_IMPLICIT, K2, t + (1d0 - gam)*dt, 0.5d0*dt)

   ! Explicit update to t+dt
   call ml_memAddToVars(MOL_EVOLVED, 0d0, 2, [Y0, K0], [1d0, dt])

   call ml_postUpdate(t + dt)
   call ml_postUpdateFast(t + dt)

   ! Calculate and store explicit RHS at t+dt (weight is 1/2)
   call ml_calcRHS(MOL_RHS_EXPLICIT, K3, t + dt, 0.5d0*dt)

   !! Final linear combination
   call ml_memAddToVars(MOL_EVOLVED, 0d0, 5, [Y0, K0, K1, K2, K3], [1d0, 0.5d0*dt, 0.5d0*dt, 0.5d0*dt, 0.5d0*dt])
   call ml_postUpdate(t + dt)
   call ml_postUpdateFast(t + dt)
end subroutine ml_advance
