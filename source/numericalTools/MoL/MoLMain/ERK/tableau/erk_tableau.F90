!!****if* source/numericalTools/MoL/MoLMain/ERK/erk_tableau
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
!!    erk_tableau
!!
!!  SYNOPSIS
!!    use erk_tableau
!!
!!  DESCRIPTION
!!    Utilities for setting up a specified ERK tableau.
!!
!!    Available methods currently include (list by runtime parameter values for method):
!!       erk-euler    : Foward Euler - first order
!!       erk-rk2-heun : Huen's method second-order RK
!!       erk-rk3-ssp  : Strong-stablity preserving third-order RK
!!       erk-rk4      : Classic fourth-order RK
!!
!!    The tableau are all given in the form:
!!
!!       c_1 | a_11 ... a_1n
!!        .  |  .   .    .
!!        .  |  .    .   .
!!        .  |  .     .  .
!!       c_n | a_n1 ... a_nn
!!       -------------------
!!           | b_1  ... b_n
!!
!!    For the explicit methods given here, the matrices a_ij are strictly lower-triangular
!!
!!    All of the coefficients used for the above methods can be found at:
!!       https://en.wikipedia.org/wiki/List_of_Runge-Kutta_methods
!!
!!  NOTES
!!
!!***
module erk_tableau

   implicit none

contains

   subroutine euler_init(A, b, c, order, stages)

      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(1, 1))
      allocate (b(1))
      allocate (c(1))

      order = 1
      stages = 1

      A(1, 1) = 0d0

      b(1) = 1d0

      c(1) = 0d0
   end subroutine euler_init

   subroutine rk2_heun_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(2, 2))
      allocate (b(2))
      allocate (c(2))

      order = 2
      stages = 2

      A = 0d0
      A(2, 1) = 1d0

      b(1) = 1d0/2d0
      b(2) = 1d0/2d0

      c(1) = 0d0
      c(2) = 1d0
   end subroutine rk2_heun_init

   subroutine rk3_ssp_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(3, 3))
      allocate (b(3))
      allocate (c(3))

      order = 3
      stages = 3

      A = 0d0
      A(2, 1) = 1d0
      A(3, 1) = 1d0/4d0
      A(3, 2) = 1d0/4d0

      b(1) = 1d0/6d0
      b(2) = 1d0/6d0
      b(3) = 2d0/3d0

      c(1) = 0d0
      c(2) = 1d0
      c(3) = 1d0/2d0
   end subroutine rk3_ssp_init

   subroutine rk4_init(A, b, c, order, stages)
      implicit none

      real, allocatable, intent(out) :: A(:, :), b(:), c(:)
      integer, intent(out) :: order, stages

      allocate (A(4, 4))
      allocate (b(4))
      allocate (c(4))

      order = 4
      stages = 4

      A = 0d0
      A(2, 1) = 1d0/2d0
      A(3, 2) = 1d0/2d0
      A(4, 3) = 1d0

      b(1) = 1d0/6d0
      b(2) = 1d0/3d0
      b(3) = 1d0/3d0
      b(4) = 1d0/6d0

      c(1) = 0d0
      c(2) = 1d0/2d0
      c(3) = 1d0/2d0
      c(4) = 1d0
   end subroutine rk4_init

end module erk_tableau
