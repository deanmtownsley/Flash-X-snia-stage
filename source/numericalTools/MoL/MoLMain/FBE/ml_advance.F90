!!****if* source/numericalTools/MoL/MoLMain/FBE/ml_advance
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
   use ml_functions, only: ml_postUpdate, ml_postUpdateFast, ml_implicitUpdate
   use ml_interface, only: ml_calcRHS
   use ml_memInterface, only: ml_memAddToVars

#include "MoL.h"

   implicit none

   real, intent(in) :: t, dt

   integer :: srcs(1)
   real :: facs(1)

   srcs = MOL_RHS
   facs = dt

   call ml_calcRHS(MOL_RHS_EXPLICIT, MOL_RHS, t)

   call ml_memAddToVars(MOL_EVOLVED, 1d0, 1, srcs, facs)

   call ml_postUpdate(t + dt)
   call ml_postUpdateFast(t + dt)

   call ml_implicitUpdate(t + dt, dt)
end subroutine ml_advance
