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
!! @brief Interfaces for MoL-registerable procedures

!> @ingroup MoL
!! Provides interfaces for MoL-registerable procedures
module MoL_functionTypes

   implicit none

   abstract interface
      !> @brief A RHS procedure
      !! @param t The time to calculate the RHS at
      subroutine MoL_rhs_t(t, activeRHS, dtWeight)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
         real, intent(in) :: dtWeight
      end subroutine MoL_rhs_t
   end interface

   abstract interface
      !> @brief An update procedure
      !! @param t  The time at the start of the update
      !! @param dt The size of the timestep to take
      subroutine MoL_update_t(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine MoL_update_t
   end interface

   abstract interface
      !> @brief A post-udate procedure
      !! @param t The time to perform the post-update at
      subroutine MoL_postUpdate_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine MoL_postUpdate_t
   end interface
end module MoL_functionTypes
