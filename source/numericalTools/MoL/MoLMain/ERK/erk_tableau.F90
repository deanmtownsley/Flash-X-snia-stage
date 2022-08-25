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
!!    Available methods currently include (list by runtime parameter values for erk_method):
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

   subroutine euler_init()
      use erk_data, only: erk_A, erk_b, erk_c, erk_order, erk_stages

      implicit none

      allocate (erk_A(1, 1))
      allocate (erk_b(1))
      allocate (erk_c(1))

      erk_order = 1
      erk_stages = 1

      erk_A(1, 1) = 0d0

      erk_b(1) = 1d0

      erk_c(1) = 0d0
   end subroutine euler_init

   subroutine rk2_heun_init()
      use erk_data, only: erk_A, erk_b, erk_c, erk_order, erk_stages

      implicit none

      allocate (erk_A(2, 2))
      allocate (erk_b(2))
      allocate (erk_c(2))

      erk_order = 2
      erk_stages = 2

      erk_A = 0d0
      erk_A(2, 1) = 1d0

      erk_b(1) = 1d0/2d0
      erk_b(2) = 1d0/2d0

      erk_c(1) = 0d0
      erk_c(2) = 1d0
   end subroutine rk2_heun_init

   subroutine rk3_ssp_init()
      use erk_data, only: erk_A, erk_b, erk_c, erk_order, erk_stages

      implicit none

      allocate (erk_A(3, 3))
      allocate (erk_b(3))
      allocate (erk_c(3))

      erk_order = 3
      erk_stages = 3

      erk_A = 0d0
      erk_A(2, 1) = 1d0
      erk_A(3, 1) = 1d0/4d0
      erk_A(3, 2) = 1d0/4d0

      erk_b(1) = 1d0/6d0
      erk_b(2) = 1d0/6d0
      erk_b(3) = 2d0/3d0

      erk_c(1) = 0d0
      erk_c(2) = 1d0
      erk_c(3) = 1d0/2d0
   end subroutine rk3_ssp_init

   subroutine rk4_init()
      use erk_data, only: erk_A, erk_b, erk_c, erk_order, erk_stages

      implicit none

      allocate (erk_A(4, 4))
      allocate (erk_b(4))
      allocate (erk_c(4))

      erk_order = 4
      erk_stages = 4

      erk_A = 0d0
      erk_A(2, 1) = 1d0/2d0
      erk_A(3, 2) = 1d0/2d0
      erk_A(4, 3) = 1d0

      erk_b(1) = 1d0/6d0
      erk_b(2) = 1d0/3d0
      erk_b(3) = 1d0/3d0
      erk_b(4) = 1d0/6d0

      erk_c(1) = 0d0
      erk_c(2) = 1d0/2d0
      erk_c(3) = 1d0/2d0
      erk_c(4) = 1d0
   end subroutine rk4_init

   subroutine erk_tableau_init
      use erk_data, only: erk_A, erk_b, erk_c, erk_method
      use ml_interface, only: ml_error

      implicit none

      if (allocated(erk_A)) deallocate (erk_A)
      if (allocated(erk_b)) deallocate (erk_b)
      if (allocated(erk_c)) deallocate (erk_c)

      select case (erk_method)

      case ("erk-euler")
         call euler_init

      case ("erk-rk2-heun")
         call rk2_heun_init

      case ("erk-rk3-ssp")
         call rk3_ssp_init

      case ("erk-rk4")
         call rk4_init

      case default
         call ml_error("Unkown ERK method")
      end select

   end subroutine erk_tableau_init

end module erk_tableau
