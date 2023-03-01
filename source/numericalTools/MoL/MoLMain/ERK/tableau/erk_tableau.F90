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
!! @brief ERK tableau utilities

!> @ingroup MoLERK
!! @brief Utilities for setting up an ERK tableau
!!
!! @details
!!    Available methods currently include (list by runtime parameter values for method):
!!       - `erk-euler`    : Foward Euler - first order
!!       - `erk-rk2-heun` : Huen's method second-order RK
!!       - `erk-rk3-ssp`  : Strong-stablity preserving third-order RK
!!       - `erk-rk4`      : Classic fourth-order RK
!!
!!    The tableau are all given in the form:
!!
!!    @f[
!!       \begin{array}{c|ccc}
!!          c_1    & A_{11} & \dots & A_{1n}\\
!!          \vdots & \vdots & \ddots& \vdots\\
!!          c_n    & A_{n1} & \dots & A_{nn}\\
!!          \hline
!!                 & b_1    & \dots & b_n\\
!!       \end{array}
!!    @f]
!!
!!    For the explicit methods given here, the matrices A_ij are strictly lower-triangular
!!
!!    All of the coefficients used for the above methods can be found at:
!!       https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
module erk_tableau

   implicit none

contains

   !> First-order Euler method
   !!
   !! @param A      Per-stage weighting coefficients
   !! @param b      Final weighting coefficients of intermediate states
   !! @param c      Timing coefficients
   !! @param order  Order of the method
   !! @param stages Number of stages
   subroutine euler_init(A, b, c, order, stages)

      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(1, 1))
      allocate (b(1))
      allocate (c(1))

      order = 1
      stages = 1

      A(1, 1) = 0.0

      b(1) = 1.0

      c(1) = 0.0
   end subroutine euler_init

   !> Second-order Heun's method
   !! @copydetails euler_init
   subroutine rk2_heun_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(2, 2))
      allocate (b(2))
      allocate (c(2))

      order = 2
      stages = 2

      A = 0.0
      A(2, 1) = 1.0

      b(1) = 1.0/2.0
      b(2) = 1.0/2.0

      c(1) = 0.0
      c(2) = 1.0
   end subroutine rk2_heun_init

   !> Third-order RK SSP method
   !! @copydetails euler_init
   subroutine rk3_ssp_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(3, 3))
      allocate (b(3))
      allocate (c(3))

      order = 3
      stages = 3

      A = 0.0
      A(2, 1) = 1.0
      A(3, 1) = 1.0/4.0
      A(3, 2) = 1.0/4.0

      b(1) = 1.0/6.0
      b(2) = 1.0/6.0
      b(3) = 2.0/3.0

      c(1) = 0.0
      c(2) = 1.0
      c(3) = 1.0/2.0
   end subroutine rk3_ssp_init

   !> Fourth-order RK method
   !! @copydetails euler_init
   subroutine rk4_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(4, 4))
      allocate (b(4))
      allocate (c(4))

      order = 4
      stages = 4

      A = 0.0
      A(2, 1) = 1.0/2.0
      A(3, 2) = 1.0/2.0
      A(4, 3) = 1.0

      b(1) = 1.0/6.0
      b(2) = 1.0/3.0
      b(3) = 1.0/3.0
      b(4) = 1.0/6.0

      c(1) = 0.0
      c(2) = 1.0/2.0
      c(3) = 1.0/2.0
      c(4) = 1.0
   end subroutine rk4_init

end module erk_tableau
