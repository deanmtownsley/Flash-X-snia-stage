!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! @brief IMEX tableau utilities

!> @ingroup MoLIMEX
!! @brief Utilities for setting up an ERK tableau
!!
!! @details
!! Implements the IMEX methods from
!!
!!    Lorenzo Pareschi and Giovanni Russo,
!!    Implicit-explicit Runge-Kutta schmes and applications to hyperbolic systems with relaxation,
!!    Journal on Scientific Computing 25, 129-155 (2005)
!!    https://doi.org/10.1007/BF02728986
module imex_tableau

   implicit none

contains

   !> First-order forward-backward Euler method
   !!
   !! @param A      Per-stage weighting coefficients
   !! @param b      Final weighting coefficients of intermediate states
   !! @param c      Timing coefficients
   !! @param order  Order of the method
   !! @param stages Number of stages
   subroutine fbe_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 1
      stages = 2 ! To make this work like higher-order methods

      AI(1, :) = [0d0, 0d0]
      AI(2, :) = [0d0, 1d0]

      bI(:) = [0d0, 0d0]

      cI(:) = [0d0, 1d0]

      AE(1, :) = [0d0, 0d0]
      AE(2, :) = [1d0, 0d0]

      bE(:) = [0d0, 0d0]

      cE(:) = [0d0, 0d0]
   end subroutine fbe_init

   !> Second-order IMEX-SSP2(2,2,2) method
   !! @copydetails fbe_init
   subroutine ssp2_222_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = 1d0 - 1d0/sqrt(2d0)

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 2
      stages = 2

      AI(1, :) = [gam, 0d0]
      AI(2, :) = [1d0 - 2d0*gam, gam]

      bI(:) = [0.5d0, 0.5d0]

      cI(:) = [gam, 1d0 - gam]

      AE(1, :) = [0d0, 0d0]
      AE(2, :) = [1d0, 0d0]

      bE(:) = [0.5d0, 0.5d0]

      cE(:) = [0d0, 1d0]
   end subroutine ssp2_222_init

   !> Second-order IMEX-SSP2(3,2,2) method
   !! @copydetails fbe_init
   subroutine ssp2_322_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [0.5d0, 0d0, 0d0]
      AI(2, :) = [-0.5d0, 0.5d0, 0d0]
      AI(3, :) = [0d0, 0.5d0, 0.5d0]

      bI(:) = [0d0, 0.5d0, 0.5d0]

      cI(:) = [0.5d0, 0d0, 1d0]

      AE(1, :) = [0d0, 0d0, 0d0]
      AE(2, :) = [0d0, 0d0, 0d0]
      AE(3, :) = [0d0, 1d0, 0d0]

      bE(:) = [0d0, 0.5d0, 0.5d0]

      cE(:) = [0d0, 0d0, 1d0]
   end subroutine ssp2_322_init

   !> Second-order IMEX-SSP2(3,3,2) method
   !! @copydetails fbe_init
   subroutine ssp2_332_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [0.25d0, 0d0, 0d0]
      AI(2, :) = [0d0, 0.25d0, 0d0]
      AI(3, :) = [1d0/3d0, 1d0/3d0, 1d0/3d0]

      bI(:) = [1d0/3d0, 1d0/3d0, 1d0/3d0]

      cI(:) = [0.25d0, 0.25d0, 1d0]

      AE(1, :) = [0d0, 0d0, 0d0]
      AE(2, :) = [0.5d0, 0d0, 0d0]
      AE(3, :) = [0.5d0, 0.5d0, 0d0]

      bE(:) = [1d0/3d0, 1d0/3d0, 1d0/3d0]

      cE(:) = [0d0, 0.5d0, 1d0]
   end subroutine ssp2_332_init

   !> Second-order IMEX-SSP3(3,3,2) method
   !! @copydetails fbe_init
   subroutine ssp3_332_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = 1d0 - 1d0/sqrt(2d0)

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [gam, 0d0, 0d0]
      AI(2, :) = [1d0 - 2d0*gam, gam, 0d0]
      AI(3, :) = [0.5d0 - gam, 0d0, gam]

      bI(:) = [1d0/6d0, 1d0/6d0, 2d0/3d0]

      cI(:) = [gam, 1d0 - gam, 0.5d0]

      AE(1, :) = [0d0, 0d0, 0d0]
      AE(2, :) = [1d0, 0d0, 0d0]
      AE(3, :) = [0.25d0, 0.25d0, 0d0]

      bE(:) = [1d0/6d0, 1d0/6d0, 2d0/3d0]

      cE(:) = [0d0, 1d0, 0.5d0]
   end subroutine ssp3_332_init

   !> Second-order IMEX-SSP3(4,3,3) method
   !! @copydetails fbe_init
   subroutine ssp3_433_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: alp = 0.24169426078821d0
      real, parameter :: beta = 0.06042356519705d0
      real, parameter :: eta = 0.12915286960590d0

      allocate (AI(4, 4), AE(4, 4))
      allocate (bI(4), bE(4))
      allocate (cI(4), cE(4))

      order = 3
      stages = 4

      AI(1, :) = [alp, 0d0, 0d0, 0d0]
      AI(2, :) = [-alp, alp, 0d0, 0d0]
      AI(3, :) = [0d0, 1d0 - alp, alp, 0d0]
      AI(4, :) = [beta, eta, 0.5d0 - beta - eta - alp, alp]

      bI(:) = [0d0, 1d0/6d0, 1d0/6d0, 2d0/3d0]

      cI(:) = [alp, 0d0, 1d0, 0.5d0]

      AE(1, :) = [0d0, 0d0, 0d0, 0d0]
      AE(2, :) = [0d0, 0d0, 0d0, 0d0]
      AE(3, :) = [0d0, 1d0, 0d0, 0d0]
      AE(4, :) = [0d0, 0.25d0, 0.25d0, 0d0]

      bE(:) = [0d0, 1d0/6d0, 1d0/6d0, 2d0/3d0]

      cE(:) = [0d0, 0d0, 1d0, 0.5d0]
   end subroutine ssp3_433_init

end module imex_tableau
