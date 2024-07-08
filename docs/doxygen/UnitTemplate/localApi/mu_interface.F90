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
! As for MyUnit_interface.

!> @ingroup MyUnitPrivate
!!
!! @brief Private interface of the MyUnit unit
!!
!! @details
!! A standard Flash-X Fortran module that encapsulates the interface declarations
!! of all routine's in the MyUnit unit that are part of the unit's private
!! interface.
module mu_interface

    implicit none

    ! Since a unit's private interface is defined by the stubs, this module, whose
    ! existence is motivated by Flash-X's internal design rules, does not need
    ! much documentation.  In particular, do *not* use doxygen to document the
    ! different interfaces below.
    interface
        ...
    end interface

end module mu_interface
