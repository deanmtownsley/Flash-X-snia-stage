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
!! @brief Procedure pointers used internally by MoL

!> @ingroup MoLMain
!! Manages registere procedures for calculating RHS terms, implicit updates
!! and post-update work
module ml_functions

   use MoL_functionTypes, only: MoL_rhs_t, MoL_update_t, MoL_postUpdate_t

   implicit none

   !> @name RHS Procedures
   !! @{

   !> Explciit RHS procedure pointer
   procedure(MoL_rhs_t), pointer :: ml_rhsE => ml_rhsE_default

   !> Implicit RHS procedure pointer
   procedure(MoL_rhs_t), pointer :: ml_rhsI => ml_rhsI_default

   !> Fast RHS procedure pointer
   procedure(MoL_rhs_t), pointer :: ml_rhsF => ml_rhsF_default
   !> @}

   !> @name Update Procedures
   !! @{

   !> Implicit update procedure pointer
   procedure(MoL_update_t), pointer :: ml_implicitUpdate => ml_implicitUpdate_default

   !> @}

   !> @name Post-Update Procedures
   !! @{

   !> Post-update procedure pointer
   procedure(MoL_postUpdate_t), pointer :: ml_postUpdate => ml_postUpdate_default

   !> Post-update (fast) procedure pointer
   procedure(MoL_postUpdate_t), pointer :: ml_postUpdateFast => ml_postUpdateFast_default
   !> @}

contains

   !> Explicit RHS default procedure
   !! @see mol_functiontypes::mol_rhs_t
   !! @param t The time to evaluate the RHS at
   subroutine ml_rhsE_default(t, activeRHS, dtWeight)
      implicit none

      real, intent(in) :: t
      integer, intent(in) :: activeRHS
      real, intent(in) :: dtWeight

      return
   end subroutine ml_rhsE_default

   !> Implicit RHS default procedure
   !! @see mol_functiontypes::mol_rhs_t
   !! @param t The time to evaluate the RHS at
   subroutine ml_rhsI_default(t, activeRHS, dtWeight)
      implicit none

      real, intent(in) :: t
      integer, intent(in) :: activeRHS
      real, intent(in) :: dtWeight

      return
   end subroutine ml_rhsI_default

   !> Fast RHS default procedure
   !! @see mol_functiontypes::mol_rhs_t
   !! @param t The time to evaluate the RHS at
   subroutine ml_rhsF_default(t, activeRHS, dtWeight)
      implicit none

      real, intent(in) :: t
      integer, intent(in) :: activeRHS
      real, intent(in) :: dtWeight

      return
   end subroutine ml_rhsF_default

   !> Implicit update default procedure
   !! @see mol_functiontypes::mol_update_t
   !! @param t  The time at the start of the implicit update
   !! @param dt The size of the timestep to take
   subroutine ml_implicitUpdate_default(t, dt)
      implicit none

      real, intent(in) :: t, dt

      return
   end subroutine ml_implicitUpdate_default

   !> Post-update default procedure
   !! @see mol_functiontypes::mol_postupdate_t
   !! @param t The time to perform post-update work at
   subroutine ml_postUpdate_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_postUpdate_default

   !> Post-update (fast) default procedure
   !! @see mol_functiontypes::mol_postupdate_t
   !! @param t The time to perform post-update work at
   subroutine ml_postUpdateFast_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_postUpdateFast_default

end module ml_functions
