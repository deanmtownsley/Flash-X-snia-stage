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

!> @ingroup MyUnit
!! @anchor MyUnit_complexGeneric_noRuntime_stub
!!
!! @brief No-runtime version of @ref myunit_interface::myunit_complexgeneric "MyUnit_complexGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the MyUnit_myRoutine stub and should establish the contract for
!! this particular version of the generic interface.
!!
!! @param a   My a argument
!! @param b   My b argument
subroutine MyUnit_complexGeneric_noRuntime(a, b)
    implicit none

    ! NOTE: The full argument details are *not* included in the @param above
    ! because Doxygen can parse the code to autmatically determine such details.
    integer, intent(IN)  :: a
    real,    intent(OUT) :: b

    ...
end subroutine MyUnit_complexGeneric_noRuntime

!> @ingroup MyUnit
!! @anchor MyUnit_complexGeneric_cpuOnly_stub
!!
!! @brief CPU-only version of @ref myunit_interface::myunit_complexgeneric "MyUnit_complexGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the MyUnit_myRoutine stub and should establish the contract for
!! this particular version of the generic interface.  Documentation writers
!! should determine how to relate this content to the details in the other versions
!! so that documentation creation & maintenance is easy and successful.
!!
!! @param a1  My a1 argument
!! @param a2  My a2 argument
!! @param b   My b argument
subroutine MyUnit_complexGeneric_cpuOnly(a1, a2, b)
    implicit none

    integer, intent(IN)  :: a1
    integer, intent(IN)  :: a2
    real,    intent(OUT) :: b

    ...
end subroutine MyUnit_complexGeneric_cpuOnly
