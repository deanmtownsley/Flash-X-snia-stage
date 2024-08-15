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
! The file documentation block should be at the top of the file for legal reasons
! and should not need modification.  C-preprocessor macros can go here between
! documentation blocks if so desired.  The next documentation block should be
! placed as close to the subroutine/function as possible.
!
! @ingroup should refer to the group name associated with the current folder.
! This is declared as the first item in the @defgroup statement in the .dox file
! in the same folder as this file.
!
! The argument to @stubref must always be the routine's name.
!
! Aside from updating the name of the subroutine, no other changes should be made
! to the @brief line.
!
! Note that we don't include @param here since those are documented in the stubs.

#include "Simulation.h"

!> @ingroup MyUnitImpl
!! @stubref{MyUnit_myRoutine}
!!
!! @todo This should indicate pending work on this implementation only
!!
!! @brief Concrete implementation of MyUnit_myRoutine.
!!
!! @details
!! The general goal is to make this content as minimal as possible.  It should
!! contain only implementation-specific details that are useful for
!! developers/maintainers.  No portion of the contract should be specified here.
subroutine MyUnit_myRoutine(a, b, c)
    implicit none

    integer, intent(OUT)            :: a
    real,    intent(INOUT)          :: b
    real,    intent(IN),   optional :: c

    ...
end subroutine MyUnit_myRoutine
