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
! The file documentation header should be at the top of the file and should not
! need modification.  C-preprocessor macros can go here if so desired.  The next
! documentation block should be placed as close to the subroutine/function as
! possible.  The general rule is to put your documentation (e.g., @details,
! @param) within that block as close as possible to the subroutine/function so
! that we can see the documentation and the programmatic interface/parameter
! declarations at the same time.  This is especially useful if the documentation
! or declarations are long.  This implies that technical directives such as
! (@ingroup and @anchor should be up top).
!
! @ingroup should refer to the group name for this unit's stub group.  This is
! declared as the first item in the @defgroup statement in MyUnit.dox.
! @anchor must always be <subroutine/function name>_stub

#include "constants.h"

!> @ingroup MyUnit
!! @anchor MyUnit_myRoutine_stub
!!
!! @todo This should indicate pending work on the interface (including contract) only
!!
!! @brief <Add in one-line description>
!! 
!! @details
!! The documentation in this section should be the general, universal contract that
!! the code establishes with its users.  Developers of concrete implementations of
!! this code need to make certain that their implementations satisfy this contract.
!! If their implementation requires alteration of this contract, then the stub
!! documentation should be updated to reflect the change and developers of all
!! other implementations notified of the change.
!! 
!! This contract can state, for example,
!! * if you give me these inputs, I will give you these outputs.
!! * an input of X is exceptional and is handled by ...
!! * It is assumed that X is true before calling this.
!! * Upon termination, Y is true.
!! * Once this is called, then Y must also ...
!!
!! The contract should be light on implementation-specific details unless, for
!! example, the contract needs to specify a high-level order of operations or
!! algorithmic structure that all implementations must satisfy.
!!
!! @param a   My a argument
!! @param b   My b argument
!! @param c   My c argument
subroutine MyUnit_myRoutine(a, b, c)
    implicit none

    ! NOTE: The full argument details are *not* included in the @param above
    ! because Doxygen can parse the code to autmatically determine such details.
    integer, intent(OUT)            :: a
    real,    intent(INOUT)          :: b
    real,    intent(IN),   optional :: c

    ...
end subroutine MyUnit_myRoutine
