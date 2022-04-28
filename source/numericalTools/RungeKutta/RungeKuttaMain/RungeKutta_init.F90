!!****if* source/numericalTools/RungeKutta/RungeKuttaMain/RungeKutta_init
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
!! NAME
!!  
!!  RungeKutta_init
!!
!! SYNOPSIS
!! 
!!  call RungeKutta_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the Runge Kutta unit.
!!
!! ARGUMENTS
!!
!!***

subroutine RungeKutta_init ()

  use RungeKutta_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface,            ONLY : Driver_abort
  use rk_interface,                ONLY : rk_setButcherTableauRepository

  implicit none
!
!
!     ...Get the needed external data.
!
!
  call RuntimeParameters_get ("rk_stepSizeConfinementFactor",       rk_stepSizeConfinementFactor)
  call RuntimeParameters_get ("rk_stepSizeSafetyFactor",            rk_stepSizeSafetyFactor     )
!
!
!     ...Check for bad data.
!
!
  if (rk_maxButcherTableauDimension < 1) then
      call Driver_abort ('[RungeKutta_init] ERROR: Bad dimension for Butcher tableaus!')
  end if

  if (rk_maxNumberDependentVariables < 1) then
      call Driver_abort ('[RungeKutta_init] ERROR: No place for dependent variables!')
  end if

  if (rk_stepSizeConfinementFactor < 0.5 .or. rk_stepSizeConfinementFactor > 1.0) then
      call Driver_abort ('[RungeKutta_init] ERROR: Confinement reduction factor out of range [0.5 to 1.0]!')
  end if

  if (rk_stepSizeSafetyFactor < 0.5 .or. rk_stepSizeSafetyFactor > 1.0) then
      call Driver_abort ('[RungeKutta_init] ERROR: Step size safety factor out of range [0.5 to 1.0]!')
  end if
!
!
!     ...Calculate machine dependent paramters.
!
!
  rk_cubeRootMacheps = epsilon (1.0) ** (1.0 / 3.0)
!
!
!     ...Set the Butcher tableaus repository.
!
!
  call rk_setButcherTableauRepository ()
!
!
!    ...Ready!
!
!
  return
end subroutine RungeKutta_init
