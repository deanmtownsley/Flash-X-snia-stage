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
!! @brief `Spacetime_metric_t` derived-type implementation

!> @ingroup SpacetimeMain
!! @brief Contains a derived-type stub for managing metric-related quantities
module Spacetime_metric

   implicit none

   !> @brief Encapsulate commonly used metric quantities and operations
   !!
   !! @details
   !! Metric quantities:
   !! - Lapse function (scalar)
   !! - Shift vector (spatial, contravariant components)
   !! - Spatial metric (spatial, covariant, symmetric components)
   !! - Inverse spatial metric (spatial, contravariant, symmetric components)
   !! - Extrinsic curvature (spatial, covariant, symmetric components)
   !! - Square-root of the spatial metric determinant (scalar)
   !!
   !! Metric operations:
   !! - Raising/lowering rank-1 tensor components
   !! - Rank-1 tensor contractions
   type :: Spacetime_metric_t

      real :: alp
      real :: betax, betay, betaz
      real :: gxx, gxy, gxz, gyy, gyz, gzz
      real :: guxx, guxy, guxz, guyy, guyz, guzz
      real :: kxx, kxy, kxz, kyy, kyz, kzz
      real :: sdetg

   contains

      procedure :: raise => Spacetime_metric__raise
      procedure :: lower => Spacetime_metric__lower

      procedure :: contractUU => Spacetime_metric__contractUU
      procedure :: contractLL => Spacetime_metric__contractLL

   end type Spacetime_metric_t

contains

   elemental subroutine Spacetime_metric__raise(this, vux, vuy, vuz, vlx, vly, vlz)
      implicit none

      class(Spacetime_metric_t), intent(in) :: this
      real, intent(out) :: vux, vuy, vuz
      real, intent(in) :: vlx, vly, vlz

      vux = this%guxx*vlx + this%guxy*vly + this%guxz*vlz
      vuy = this%guxy*vlx + this%guyy*vly + this%guyz*vlz
      vuz = this%guxz*vlx + this%guyz*vly + this%guzz*vlz
   end subroutine Spacetime_metric__raise

   elemental subroutine Spacetime_metric__lower(this, vlx, vly, vlz, vux, vuy, vuz)
      implicit none

      class(Spacetime_metric_t), intent(in) :: this
      real, intent(out) :: vlx, vly, vlz
      real, intent(in) :: vux, vuy, vuz

      vlx = this%gxx*vux + this%gxy*vuy + this%gxz*vuz
      vly = this%gxy*vux + this%gyy*vuy + this%gyz*vuz
      vlz = this%gxz*vux + this%gyz*vuy + this%gzz*vuz
   end subroutine Spacetime_metric__lower

   elemental subroutine Spacetime_metric__contractUU(this, vw, vux, vuy, vuz, &
                                                     wux, wuy, wuz)
      implicit none

      class(Spacetime_metric_t), intent(in) :: this
      real, intent(out) :: vw
      real, intent(in) :: vux, vuy, vuz
      real, intent(in) :: wux, wuy, wuz

      vw = this%gxx*vux*wux + this%gyy*vuy*wuy + this%gzz*vuz*wuz &
           + this%gxy*(vux*wuy + vuy*wux) &
           + this%gxz*(vux*wuz + vuz*wux) &
           + this%gyz*(vuy*wuz + vuz*wuy)
   end subroutine Spacetime_metric__contractUU

   elemental subroutine Spacetime_metric__contractLL(this, vw, vlx, vly, &
                                                     vlz, wlx, wly, wlz)
      implicit none

      class(Spacetime_metric_t), intent(in) :: this
      real, intent(out) :: vw
      real, intent(in) :: vlx, vly, vlz
      real, intent(in) :: wlx, wly, wlz

      vw = this%guxx*vlx*wlx + this%guyy*vly*wly + this%guzz*vlz*wlz &
           + this%guxy*(vlx*wly + vly*wlx) &
           + this%guxz*(vlx*wlz + vlz*wlx) &
           + this%guyz*(vly*wlz + vlz*wly)
   end subroutine Spacetime_metric__contractLL

end module Spacetime_metric
