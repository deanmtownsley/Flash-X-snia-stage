!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

! FORMATTING RULES
!
! The file documentation header should be at the top of the file and should need
! no modification.  C-preprocessor macros can go here if so desired.  The next
! documentation block should be placed as close to the module as possible.
! 
! @ingroup should refer to the group name associated with the current folder.
! This is declared as the first item in the @defgroup statement in the .dox file
! in the same folder as this file.

!> @ingroup MyUnitMain
!!
!! @brief Flash-X standard data module for MyUnit
module MyUnit_data

    !> My variable a is
    !! great because of
    integer, save :: a

    !> @name Group of Variables
    !! @{
    !! These variables all have the same brief description
    real, save :: b
    real, save :: c
    !> @}

    !> @name Some other variable group
    !! @{

    !> Brief description of d
    logical, save :: d

    !> Brief description of e
    type(I_thought_this_was_c_t) :: e

    !> @}

    ...
end module MyUnit_data
