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

!> @ingroup MyUnitMemory
!! @stubref{mu_memMyRoutine}
!!
!! @brief Concrete implementation of mu_memMyRoutine.
!! 
!! @details
!! This section should be written following the same rules for the details
!! section of the concrete implementation of MyUnit_myRoutine.
subroutine mu_memMyRoutine(a, b, c)
    implicit none

    integer, intent(OUT)            :: a
    real,    intent(INOUT)          :: b
    real,    intent(IN),   optional :: c

    ...
end subroutine mu_memMyRoutine
