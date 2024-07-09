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
! @ingroup should refer to the group name for this unit's stub group.  This is
! declared as the first item in the @defgroup statement in MyUnit.dox.

!> @ingroup MyUnit
!!
!! @brief Public interface of the MyUnit unit
!!
!! @details
!! A standard Flash-X Fortran module that encapsulates the interface declarations
!! of all routine's in the MyUnit unit that are part of this unit's public
!! interface.
module MyUnit_interface

    implicit none

    ! Since a unit's public interface is defined by the stubs, this module, whose
    ! existence is motivated by Flash-X's internal design rules, does not need
    ! much documentation.  In particular, do *not* use doxygen to document the
    ! different interfaces below.
    interface
        subroutine MyUnit_myRoutine(a, b, c)
            implicit none
            ...
        end subroutine MyUnit_myRoutine
    end interface

    interface MyUnit_simpleGeneric
        ! These two share the same contract aside from slight differences that
        ! don't need explicit documentation/explanation.
        subroutine MyUnit_simpleGeneric_int(a, b)
            implicit none
            integer, intent(IN)  :: a
            real,    intent(OUT) :: b
        end subroutine MyUnit_simpleGeneric_int
        subroutine MyUnit_simpleGeneric_real(a, b)
            implicit none
            real, intent(IN)  :: a
            real, intent(OUT) :: b
        end subroutine MyUnit_simpleGeneric_real
    end interface MyUnit_simpleGeneric

    interface MyUnit_complexGeneric
        ! These two share a commonality but require different contracts at the
        ! very least to accommodate significant argument differences
        subroutine MyUnit_complexGeneric_noRuntime(a, b)
            implicit none
            integer, intent(IN)  :: a
            real,    intent(OUT) :: b
        end subroutine MyUnit_complexGeneric_noRuntime
        subroutine MyUnit_complexGeneric_cpuOnly(a1, a2, b)
            implicit none
            integer, intent(IN)  :: a1
            integer, intent(IN)  :: a2
            real,    intent(OUT) :: b
        end subroutine MyUnit_complexGeneric_cpuOnly
    end interface MyUnit_complexGeneric

end module MyUnit_interface
