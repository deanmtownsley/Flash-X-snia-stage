!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/Simulation_init
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
!!  Simulation_init
!!
!! SYNOPSIS
!!
!!  Simulation_init ()
!!
!! DESCRIPTION
!!
!!  Initializes the parameters for the Roots x2 polynomials unit test.
!!
!!***

subroutine Simulation_init ()

  use  Simulation_data
  use  RuntimeParameters_interface, ONLY: RuntimeParameters_get

  implicit none

  integer :: ut_getFreeFileUnit
!
!
!    ...Get necessary runtime parameters.
!
!
  call RuntimeParameters_get ("basenm",           sim_baseName )
!
!
!    ...Prepare the output file.
!
!
  sim_baseName  = adjustl (sim_baseName)     ! to get ready to use 'trim'
  sim_printUnit = ut_getFreeFileUnit ()
  sim_printName = trim (sim_baseName) // "Printout.dat"

  open (unit = sim_printUnit, file = sim_printName)
!
!
!    ...Set some constants we need.
!
!
  sim_LPN     = huge (1.0)
  sim_SPN     = tiny (1.0)
  sim_sqrtLPN = sqrt (sim_LPN)
  sim_sqrtSPN = sqrt (sim_SPN)
!
!
!    ...Set the upper bound of relative accuracies, the x2 polynomial coefficients and the
!       exact roots for comparison.
!
!
  call sim_setMaxRelativeAccuracies  ()
  call sim_setx2PolynomialCoeffs     ()
  call sim_setx2PolynomialExactRoots ()
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_init
