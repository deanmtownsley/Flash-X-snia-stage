!!****if* source/Simulation/SimulationMain/unitTest/Roots/x2Polynomials/sim_setx2PolynomialCoeffs
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
!!  sim_setx2PolynomialCoeffs
!! 
!! SYNOPSIS
!!
!!  call sim_setx2PolynomialCoeffs
!!
!! DESCRIPTION
!!
!!  Sets the x2 polynomial coefficients A,B for the x^1 and x^0 terms:
!!
!!                   x2 polynomial = x^2 + Ax + B
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx2PolynomialCoeffs ()

  use  Driver_interface, ONLY: Driver_abort

  use  Simulation_data,  ONLY: sim_LPN, sim_SPN,          &
                               sim_numberOfx2Polynomials, &
                               sim_x2Polynomialx0Coeff,   &
                               sim_x2Polynomialx1Coeff

  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx2Polynomials > 16) then
      call Driver_abort ('[sim_setx2PolynomialCoeffs] ERROR: Not enough storage for coefficients!')
  end if
!
!
!   ...Set the coefficients for: x^2 + ax + b.
!
!
  sim_x2Polynomialx1Coeff (1:16) = (/ + sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      - sim_SPN  /)

  sim_x2Polynomialx0Coeff (1:16) = (/ + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_LPN, &
                                      - sim_LPN, &
                                      + sim_SPN, &
                                      - sim_SPN, &
                                      + sim_SPN, &
                                      - sim_SPN  /)
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx2PolynomialCoeffs
