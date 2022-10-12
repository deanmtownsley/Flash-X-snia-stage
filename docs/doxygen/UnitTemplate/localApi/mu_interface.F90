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

As for MyUnit_interface.  The documentation of stubs and concrete
implementations in the private interface should follow the rules for documenting
stubs and concrete implementations in the public interface.

!> @ingroup MyUnitPrivate
!!
!! @brief Private interface of the MyUnit unit
!!
!! @details
!! A standard Flash-X Fortran module that encapsulates the interface declarations
!! of all routine's in the MyUnit unit that are part of this unit's private
!! interface.
module mu_interface

    implicit none

    interface
        ...
    end interface

end module mu_interface
