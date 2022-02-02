!!****if* source/Simulation/SimulationMain/unitTest/Roots/x3Polynomials/sim_setx3PolynomialCoeffs
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
!!  sim_setx3PolynomialCoeffs
!! 
!! SYNOPSIS
!!
!!  call sim_setx3PolynomialCoeffs
!!
!! DESCRIPTION
!!
!!  Sets the x3 polynomial coefficients A,B,C for the x^2, x^1 and x^0 terms:
!!
!!              x3 polynomial = x^3 + Ax^2 + Bx + C
!!
!! ARGUMENTS 
!!
!!***

subroutine sim_setx3PolynomialCoeffs ()

  use  Driver_interface, ONLY: Driver_abort

  use  Simulation_data,  ONLY: sim_numberOfx3Polynomials, &
                               sim_x3Polynomialx0Coeff,   &
                               sim_x3Polynomialx1Coeff,   &
                               sim_x3Polynomialx2Coeff

  implicit none
!
!
!   ...Check array storage size.
!
!
  if (sim_numberOfx3Polynomials > 9) then
      call Driver_abort ('[sim_setx3PolynomialCoeffs] ERROR: Not enough storage for coefficients!')
  end if
!
!
!   ...Set the coefficients.
!
!
  sim_x3Polynomialx2Coeff (1:9) = (/ -1.00000010000001e+14, &
                                     -3.000003,             &
                                     +8.999e+80,            &
                                     +1.e+24,               &
                                     -1.99999999999999e+14, &
                                     -3.e+5,                &
                                     -3.0,                  &
                                     -2.0000001e+7,         &
                                     +0.99999999999998e+14  /)

  sim_x3Polynomialx1Coeff (1:9) = (/ +1.00000010000001e+21, &
                                     +3.000006000002,       &
                                     -1.0009e+161,          &
                                     -1.0,                  &
                                     +0.99999999999998e+28, &
                                     +3.0000000001e+10,     &
                                     +1.00000000000003e+14, &
                                     +1.00000020000001e+14, &
                                     -1.99999999999998e+14  /)

  sim_x3Polynomialx0Coeff (1:9) = (/ -1.e+21,               &
                                     -1.000003000002,       &
                                     +1.e+238,              &
                                     -1.e+24,               &
                                     +1.e+28,               &
                                     -1.0000000001e+15,     &
                                     -1.00000000000001e+14, &
                                     -1.00000000000001e+14, &
                                     +2.e+14                /)
!
!
!   ...Ready!
!
!
  return
end subroutine sim_setx3PolynomialCoeffs
