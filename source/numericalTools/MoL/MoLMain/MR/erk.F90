!!****if* source/numericalTools/MoL/MoLMain/MR/erk
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
!!    erk
!!
!!  SYNOPSIS
!!    use erk
!!
!!  DESCRIPTION
!!    Utilities for setting up a specified ERK tableau for "fast" integration steps.
!!
!!    Available methods currently include (list by runtime parameter values for mr_fastMethod):
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
!!    The routines provided below are identical to those in MoL/MoLMain/ERK in all but the names
!!    of the MR-specific modules and variables
!!
!!***
module erk

   implicit none

contains

   subroutine euler_init()
      use mr_data, only: mr_AF, mr_bF, mr_cF, mr_nstages_fast

      implicit none

      allocate (mr_AF(1, 1))
      allocate (mr_bF(1))
      allocate (mr_cF(1))

      mr_nstages_fast = 1

      mr_AF(1, 1) = 0d0

      mr_bF(1) = 1d0

      mr_cF(1) = 0d0
   end subroutine euler_init

   subroutine rk2_heun_init()
      use mr_data, only: mr_AF, mr_bF, mr_cF, mr_nstages_fast

      implicit none

      allocate (mr_AF(2, 2))
      allocate (mr_bF(2))
      allocate (mr_cF(2))

      mr_nstages_fast = 2

      mr_AF = 0d0
      mr_AF(2, 1) = 1d0

      mr_bF(1) = 1d0/2d0
      mr_bF(2) = 1d0/2d0

      mr_cF(1) = 0d0
      mr_cF(2) = 1d0
   end subroutine rk2_heun_init

   subroutine rk3_ssp_init()
      use mr_data, only: mr_AF, mr_bF, mr_cF, mr_nstages_fast

      implicit none

      allocate (mr_AF(3, 3))
      allocate (mr_bF(3))
      allocate (mr_cF(3))

      mr_nstages_fast = 3

      mr_AF = 0d0
      mr_AF(2, 1) = 1d0
      mr_AF(3, 1) = 1d0/4d0
      mr_AF(3, 2) = 1d0/4d0

      mr_bF(1) = 1d0/6d0
      mr_bF(2) = 1d0/6d0
      mr_bF(3) = 2d0/3d0

      mr_cF(1) = 0d0
      mr_cF(2) = 1d0
      mr_cF(3) = 1d0/2d0
   end subroutine rk3_ssp_init

   subroutine rk4_init()
      use mr_data, only: mr_AF, mr_bF, mr_cF, mr_nstages_fast

      implicit none

      allocate (mr_AF(4, 4))
      allocate (mr_bF(4))
      allocate (mr_cF(4))

      mr_nstages_fast = 4

      mr_AF = 0d0
      mr_AF(2, 1) = 1d0/2d0
      mr_AF(3, 2) = 1d0/2d0
      mr_AF(4, 3) = 1d0

      mr_bF(1) = 1d0/6d0
      mr_bF(2) = 1d0/3d0
      mr_bF(3) = 1d0/3d0
      mr_bF(4) = 1d0/6d0

      mr_cF(1) = 0d0
      mr_cF(2) = 1d0/2d0
      mr_cF(3) = 1d0/2d0
      mr_cF(4) = 1d0
   end subroutine rk4_init

   subroutine erk_init
      use mr_data, only: mr_fastMethod
      use ml_interface, only: ml_error

      implicit none

      select case (mr_fastMethod)

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

   end subroutine erk_init

end module erk
