!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Public interfaces for the turbulence driving unit Stir
!!
!! @details This is the header file for the turbulence driving module Stir
!!          that defines its public interfaces.

!> @ingroup Stir
!! Interfaces to Stir public procedures
Module Stir_interface

#include "constants.h"

  interface Stir
    subroutine Stir(dt, pass)
      real,intent(IN) :: dt
      integer, intent(in), optional :: pass
    end subroutine Stir
  end interface

  interface Stir_finalize
    subroutine Stir_finalize()
    end subroutine Stir_finalize
  end interface

  interface Stir_init
    subroutine Stir_init(restart)
      logical, intent(IN) :: restart
    end subroutine Stir_init
  end interface

end Module Stir_interface
