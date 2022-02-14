!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/sim_setMaxRelativeAccuracies
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
!!  sim_setMaxRelativeAccuracies
!! 
!! SYNOPSIS
!!
!!  call sim_setMaxRelativeAccuracies
!!
!! DESCRIPTION
!!
!!  Sets the maximum allowed relative root accuracy for each x3 polynomial. This is a global
!!  accuracy for each x3 polynomial, in the sense that this maximum applies for all the 3 roots
!!  and all the real and imaginary components of the roots.
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setMaxRelativeAccuracies ()

  use  Driver_interface, ONLY: Driver_abort

  use  Simulation_data,  ONLY: sim_maxRelativeAccuracy,   &
                               sim_numberOfx3Polynomials
  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx3Polynomials > 9) then
      call Driver_abort ('[sim_setMaxRelativeAccuracies] ERROR: Not enough storage for max accuracies!')
  end if
!
!
!   ...Set the maximum relative accuracy for each x3 polynomial.
!
!
  sim_maxRelativeAccuracy (1)  = 1.0e-15
  sim_maxRelativeAccuracy (2)  = 1.0e-9
  sim_maxRelativeAccuracy (3)  = 1.0e-15
  sim_maxRelativeAccuracy (4)  = 1.0e-15
  sim_maxRelativeAccuracy (5)  = 1.0e-15
  sim_maxRelativeAccuracy (6)  = 1.0e-15
  sim_maxRelativeAccuracy (7)  = 1.0e-15
  sim_maxRelativeAccuracy (8)  = 1.0e-15
  sim_maxRelativeAccuracy (9)  = 1.0e-15
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setMaxRelativeAccuracies
