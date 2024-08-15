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
! Analogous to those for the MyUnit_myRoutine stub.

!> @ingroup MyUnitPrivate
!! @anchor mu_memMyRoutine_stub
!!
!! @brief <Add in one-line description>
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the MyUnit_myRoutine stub.
!!
!! @param a   My a argument
!! @param b   My b argument
!! @param c   My c argument
subroutine mu_memMyRoutine(a, b, c)
    implicit none

    ! NOTE: The full argument details are *not* included in the @param above
    ! because Doxygen can parse the code to autmatically determine such details.
    integer, intent(OUT)            :: a
    real,    intent(INOUT)          :: b
    real,    intent(IN),   optional :: c

    ...
end subroutine mu_memMyRoutine
