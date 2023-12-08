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

!> @ingroup MyUnitMain
!! @stubref{MyUnit_simpleGeneric_int}
!!
!! @brief Concrete integer version of @ref myunit_interface::myunit_simplegeneric "MyUnit_simpleGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the concrete implementation of MyUnit_myRoutine.
subroutine MyUnit_simpleGeneric_int(a, b)
    implicit none

    integer, intent(IN)  :: a
    real,    intent(OUT) :: b

    ...
end subroutine MyUnit_simpleGeneric_int

!> @ingroup MyUnitMain
!! @stubref{MyUnit_simpleGeneric_real}
!!
!! @brief Concrete real version of @ref myunit_interface::myunit_simplegeneric "MyUnit_simpleGeneric generic interface"
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the concrete implementation of MyUnit_myRoutine.  Note that the
!! implementation of the real version might require some dev/maintainer docs that
!! are different from the integer version.
subroutine MyUnit_simpleGeneric_real(a, b)
    implicit none

    real, intent(IN)  :: a
    real, intent(OUT) :: b

    ...
end subroutine MyUnit_simpleGeneric_real
