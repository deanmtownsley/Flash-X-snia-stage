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
!! @brief MoL_registerVariable implementation

!> @ingroup MoLMain
!!
!! @brief Implements MoL_registerVariable
!!
!! @stubref{MoL_registerVariable}
subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
   use ml_variables, only: ml_variable_t, ml_vars, ml_nvars, &
                           ml_unk_mask, ml_scratch_mask, ml_unk_to_scratch
   use ml_interface, only: ml_warn

   implicit none

   character(len=*), intent(in) :: name
   integer, intent(in) :: evolIndex
   integer, intent(out) :: rhsIndex

   type(ml_variable_t) :: var

   integer :: ivar

   ! If this variable is registered already, output a warning (if verbose) and exit
   varLoop: do ivar = 1, ml_nvars
      if (ml_vars(ivar)%evolIndex .eq. evolIndex) then
         call ml_warn("Duplicate variable registration")
         return
      end if
   end do varLoop

   ml_nvars = ml_nvars + 1

   var%name = trim(name)
   var%evolIndex = evolIndex
   var%rhsIndex = ml_nvars

   rhsIndex = ml_nvars

   ml_unk_to_scratch(evolIndex) = rhsIndex

   if (.not. allocated(ml_vars)) then
      allocate (ml_vars(1))
      ml_vars = var
   else
      ml_vars = (/ml_vars, var/)
   end if

   if (.not. allocated(ml_unk_mask)) then
      allocate (ml_unk_mask(1))
      ml_unk_mask = evolIndex
   else
      ml_unk_mask = (/ml_unk_mask, evolIndex/)
   end if

   if (.not. allocated(ml_scratch_mask)) then
      allocate (ml_scratch_mask(1))
      ml_scratch_mask = ml_nvars
   else
      ml_scratch_mask = (/ml_scratch_mask, ml_nvars/)
   end if
end subroutine MoL_registerVariable
