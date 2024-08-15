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
! Analogous to those for the concrete implementation of MyUnit_myRoutine.

!> @ingroup MyUnitImpl
!! @stubref{MyUnit_complexGeneric_noRuntime}
!!
!! @brief Concrete no-runtime version of @ref myunit_interface::myunit_complexgeneric "MyUnit_complexGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the concrete implementation of MyUnit_myRoutine.
subroutine MyUnit_complexGeneric_noRuntime(a, b)
    implicit none

    integer, intent(IN)  :: a
    real,    intent(OUT) :: b

    ...
end subroutine MyUnit_complexGeneric_noRuntime

!> @ingroup MyUnitImpl
!! @stubref{MyUnit_complexGeneric_cpuOnly}
!!
!! @brief Concrete CPU-only version of @ref myunit_interface::myunit_complexgeneric "MyUnit_complexGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the concrete implementation of MyUnit_myRoutine.  Note that the
!! cpu-only implementation might require some dev/maintainer docs that
!! are different from the no-runtime version.
subroutine MyUnit_complexGeneric_cpuOnly(a1, a2, b)
    implicit none

    integer, intent(IN)  :: a1
    integer, intent(IN)  :: a2
    real,    intent(OUT) :: b

    ...
end subroutine MyUnit_complexGeneric_cpuOnly
