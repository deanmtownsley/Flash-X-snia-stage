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
!!
!! and
!!
!!    Uri M. Ascher, Steven J. Ruuth, and Raymond J. Spiteri
!!    Implicit-Explicit Runge-Kutta methods for time-dependent partial differential equations
!!    Applied Numerical Mathematics 25, 151-167 (1997)
!!    https://doi.org/10.1016/S0168-9274(97)00056-1
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

      AI(1, :) = [0.0, 0.0]
      AI(2, :) = [0.0, 1.0]

      bI(:) = [0.0, 1.0]

      cI(:) = [0.0, 1.0]

      AE(1, :) = [0.0, 0.0]
      AE(2, :) = [1.0, 0.0]

      bE(:) = [1.0, 0.0]

      cE(:) = [0.0, 1.0]
   end subroutine fbe_init

   !> Second-order IMEX-SSP2(2,2,2) method
   !! @copydetails fbe_init
   subroutine ssp2_222_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = 1.0 - 1.0/sqrt(2.0)

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 2
      stages = 2

      AI(1, :) = [gam, 0.0]
      AI(2, :) = [1.0 - 2.0*gam, gam]

      bI(:) = [0.5, 0.5]

      cI(:) = [gam, 1.0 - gam]

      AE(1, :) = [0.0, 0.0]
      AE(2, :) = [1.0, 0.0]

      bE(:) = [0.5, 0.5]

      cE(:) = [0.0, 1.0]
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

      AI(1, :) = [0.5, 0.0, 0.0]
      AI(2, :) = [-0.5, 0.5, 0.0]
      AI(3, :) = [0.0, 0.5, 0.5]

      bI(:) = [0.0, 0.5, 0.5]

      cI(:) = [0.5, 0.0, 1.0]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [0.0, 0.0, 0.0]
      AE(3, :) = [0.0, 1.0, 0.0]

      bE(:) = [0.0, 0.5, 0.5]

      cE(:) = [0.0, 0.0, 1.0]
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

      AI(1, :) = [0.25, 0.0, 0.0]
      AI(2, :) = [0.0, 0.25, 0.0]
      AI(3, :) = [1.0/3.0, 1.0/3.0, 1.0/3.0]

      bI(:) = [1.0/3.0, 1.0/3.0, 1.0/3.0]

      cI(:) = [0.25, 0.25, 1.0]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [0.5, 0.0, 0.0]
      AE(3, :) = [0.5, 0.5, 0.0]

      bE(:) = [1.0/3.0, 1.0/3.0, 1.0/3.0]

      cE(:) = [0.0, 0.5, 1.0]
   end subroutine ssp2_332_init

   !> Second-order IMEX-SSP3(3,3,2) method
   !! @copydetails fbe_init
   subroutine ssp3_332_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = 1.0 - 1.0/sqrt(2.0)

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [gam, 0.0, 0.0]
      AI(2, :) = [1.0 - 2.0*gam, gam, 0.0]
      AI(3, :) = [0.5 - gam, 0.0, gam]

      bI(:) = [1.0/6.0, 1.0/6.0, 2.0/3.0]

      cI(:) = [gam, 1.0 - gam, 0.5]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [1.0, 0.0, 0.0]
      AE(3, :) = [0.25, 0.25, 0.0]

      bE(:) = [1.0/6.0, 1.0/6.0, 2.0/3.0]

      cE(:) = [0.0, 1.0, 0.5]
   end subroutine ssp3_332_init

   !> Second-order IMEX-SSP3(4,3,3) method
   !! @copydetails fbe_init
   subroutine ssp3_433_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: alp = 0.24169426078821
      real, parameter :: beta = 0.06042356519705
      real, parameter :: eta = 0.12915286960590

      allocate (AI(4, 4), AE(4, 4))
      allocate (bI(4), bE(4))
      allocate (cI(4), cE(4))

      order = 3
      stages = 4

      AI(1, :) = [alp, 0.0, 0.0, 0.0]
      AI(2, :) = [-alp, alp, 0.0, 0.0]
      AI(3, :) = [0.0, 1.0 - alp, alp, 0.0]
      AI(4, :) = [beta, eta, 0.5 - beta - eta - alp, alp]

      bI(:) = [0.0, 1.0/6.0, 1.0/6.0, 2.0/3.0]

      cI(:) = [alp, 0.0, 1.0, 0.5]

      AE(1, :) = [0.0, 0.0, 0.0, 0.0]
      AE(2, :) = [0.0, 0.0, 0.0, 0.0]
      AE(3, :) = [0.0, 1.0, 0.0, 0.0]
      AE(4, :) = [0.0, 0.25, 0.25, 0.0]

      bE(:) = [0.0, 1.0/6.0, 1.0/6.0, 2.0/3.0]

      cE(:) = [0.0, 0.0, 1.0, 0.5]
   end subroutine ssp3_433_init

   !> First-order forward-backward Euler method (1,1,1)
   !! @copydetails fbe_init
   subroutine ark_111_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 1
      stages = 2

      AI(1, :) = [0.0, 0.0]
      AI(2, :) = [0.0, 1.0]

      bI(:) = [0.0, 1.0]

      cI(:) = [0.0, 1.0]

      AE(1, :) = [0.0, 0.0]
      AE(2, :) = [1.0, 0.0]

      bE(:) = [1.0, 0.0]

      cE(:) = [0.0, 1.0]
   end subroutine ark_111_init

   !> First-order forward-backward Euler method (1,1,1)
   !! @copydetails fbe_init
   subroutine ark_121_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 1
      stages = 2

      AI(1, :) = [0.0, 0.0]
      AI(2, :) = [0.0, 1.0]

      bI(:) = [0.0, 1.0]

      cI(:) = [0.0, 1.0]

      AE(1, :) = [0.0, 0.0]
      AE(2, :) = [1.0, 0.0]

      bE(:) = [0.0, 1.0]

      cE(:) = [0.0, 1.0]
   end subroutine ark_121_init

   !> Second-order mid-point method (1,2,2)
   !! @copydetails fbe_init
   subroutine ark_122_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      allocate (AI(2, 2), AE(2, 2))
      allocate (bI(2), bE(2))
      allocate (cI(2), cE(2))

      order = 2
      stages = 2

      AI(1, :) = [0.0, 0.0]
      AI(2, :) = [0.0, 0.5]

      bI(:) = [0.0, 1.0]

      cI(:) = [0.0, 0.5]

      AE(1, :) = [0.0, 0.0]
      AE(2, :) = [0.5, 0.0]

      bE(:) = [0.0, 1.0]

      cE(:) = [0.0, 0.5]
   end subroutine ark_122_init

   !> Second-order method (2,2,2)
   !! @copydetails fbe_init
   subroutine ark_222_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = (2.0 + sqrt(2.0))/2.0
      real, parameter :: d = 1.0 - 1.0/(2.0*gam)

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [0.0, 0.0, 0.0]
      AI(2, :) = [0.0, gam, 0.0]
      AI(3, :) = [0.0, 1.0 - gam, gam]

      bI(:) = [0.0, 1.0 - gam, gam]

      cI(:) = [0.0, gam, 1.0]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [gam, 0.0, 0.0]
      AE(3, :) = [d, 1.0 - d, 0.0]

      bE(:) = [d, 1.0 - d, 0.0]

      cE(:) = [0.0, gam, 1.0]
   end subroutine ark_222_init

   !> Second-order method (2,3,2)
   !! @copydetails fbe_init
   subroutine ark_232_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = (2.0 + sqrt(2.0))/2.0
      real, parameter :: d = -2.0*sqrt(2.0)/3.0

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 2
      stages = 3

      AI(1, :) = [0.0, 0.0, 0.0]
      AI(2, :) = [0.0, gam, 0.0]
      AI(3, :) = [0.0, 1.0 - gam, gam]

      bI(:) = [0.0, 1.0 - gam, gam]

      cI(:) = [0.0, gam, 1.0]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [gam, 0.0, 0.0]
      AE(3, :) = [d, 1.0 - d, 0.0]

      bE(:) = [0.0, 1.0 - gam, gam]

      cE(:) = [0.0, gam, 1.0]
   end subroutine ark_232_init

   !> Third-order method (2,3,3)
   !! @copydetails fbe_init
   subroutine ark_233_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: gam = (3.0 + sqrt(3.0))/6.0

      allocate (AI(3, 3), AE(3, 3))
      allocate (bI(3), bE(3))
      allocate (cI(3), cE(3))

      order = 3
      stages = 3

      AI(1, :) = [0.0, 0.0, 0.0]
      AI(2, :) = [0.0, gam, 0.0]
      AI(3, :) = [0.0, 1.0 - 2.0*gam, gam]

      bI(:) = [0.0, 0.5, 0.5]

      cI(:) = [0.0, gam, 1.0 - gam]

      AE(1, :) = [0.0, 0.0, 0.0]
      AE(2, :) = [gam, 0.0, 0.0]
      AE(3, :) = [gam - 1.0, 2.0*(1.0 - gam), 0.0]

      bE(:) = [0.0, 0.5, 0.5]

      cE(:) = [0.0, gam, 1.0 - gam]
   end subroutine ark_233_init

   !> Third-order method (3,4,3)
   !! @copydetails fbe_init
   subroutine ark_343_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: eta = 0.4358665215084589994160194511935568425293
      real, parameter :: alp = 0.5529291480359398193611887297385924764949
      real, parameter :: a31 = alp*(15.0 - 60.0*eta + 21.0*eta**2)/4.0 - 3.5 + 13.0*eta - 4.5*eta**2
      real, parameter :: a32 = alp*(-15.0 + 60.0*eta - 21.0*eta**2)/4.0 + 4.0 - 12.5*eta + 4.5*eta**2
      real, parameter :: b2 = -1.5*eta**2 + 4.0*eta - 0.25
      real, parameter :: b3 = 1.5*eta**2 - 5.0*eta + 1.25

      allocate (AI(4, 4), AE(4, 4))
      allocate (bI(4), bE(4))
      allocate (cI(4), cE(4))

      order = 3
      stages = 4

      AI(1, :) = [0.0, 0.0, 0.0, 0.0]
      AI(2, :) = [0.0, eta, 0.0, 0.0]
      AI(3, :) = [0.0, 0.5*(1.0 - eta), eta, 0.0]
      AI(4, :) = [0.0, b2, b3, eta]

      bI(:) = [0.0, b2, b3, eta]

      cI(:) = [0.0, eta, 0.5*(1.0 + eta), 1.0]

      AE(1, :) = [0.0, 0.0, 0.0, 0.0]
      AE(2, :) = [eta, 0.0, 0.0, 0.0]
      AE(3, :) = [a31, a32, 0.0, 0.0]
      AE(4, :) = [1.0 - 2.0*alp, alp, alp, 0.0]

      bE(:) = [0.0, b2, b3, eta]

      cE(:) = [0.0, eta, 0.5*(1.0 + eta), 1.0]
   end subroutine ark_343_init

   !> Third-order method (4,4,3)
   !! @copydetails fbe_init
   subroutine ark_443_init(AI, bI, cI, AE, bE, cE, order, stages)

      implicit none

      real, allocatable, intent(out) :: AI(:, :), bI(:), cI(:)
      real, allocatable, intent(out) :: AE(:, :), bE(:), cE(:)
      integer, intent(out) :: order, stages

      real, parameter :: eta = 0.4358665215084589994160194511935568425293
      real, parameter :: alp = 0.5529291480359398193611887297385924764949
      real, parameter :: a31 = alp*(15.0 - 60.0*eta + 21.0*eta**2)/4.0 - 3.5 + 13.0*eta - 4.5*eta**2
      real, parameter :: a32 = alp*(-15.0 + 60.0*eta - 21.0*eta**2)/4.0 + 4.0 - 12.5*eta + 4.5*eta**2
      real, parameter :: b2 = -1.5*eta**2 + 4.0*eta - 0.25
      real, parameter :: b3 = 1.5*eta**2 - 5.0*eta + 1.25

      allocate (AI(5, 5), AE(5, 5))
      allocate (bI(5), bE(5))
      allocate (cI(5), cE(5))

      order = 3
      stages = 5

      AI(1, :) = [0.0, 0.0, 0.0, 0.0, 0.0]
      AI(2, :) = [0.0, 0.5, 0.0, 0.0, 0.0]
      AI(3, :) = [0.0, 1.0/6.0, 0.5, 0.0, 0.0]
      AI(4, :) = [0.0, -0.5, 0.5, 0.5, 0.0]
      AI(5, :) = [0.0, 1.5, -1.5, 0.5, 0.5]

      bI(:) = [0.0, 1.5, -1.5, 0.5, 0.5]

      cI(:) = [0.0, 0.5, 2.0/3.0, 0.5, 1.0]

      AE(1, :) = [0.0, 0.0, 0.0, 0.0, 0.0]
      AE(2, :) = [0.5, 0.0, 0.0, 0.0, 0.0]
      AE(3, :) = [11.0/18.0, 1.0/18.0, 0.0, 0.0, 0.0]
      AE(4, :) = [5.0/6.0, -5.0/6.0, 0.5, 0.0, 0.0]
      AE(5, :) = [0.25, 1.75, 0.75, -1.75, 0.0]

      bE(:) = [0.25, 1.75, 0.75, -1.75, 0.0]

      cE(:) = [0.0, 0.5, 2.0/3.0, 0.5, 1.0]
   end subroutine ark_443_init

end module imex_tableau
