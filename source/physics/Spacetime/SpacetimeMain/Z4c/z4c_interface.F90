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
!! @brief Private interfaces for the Z4c solver
!!
!! @details This is the header file for the Z4c solver that defines
!!          its private interfaces.

!> @ingroup Z4c
!! Interfaces to Z4c private procedures
module z4c_interface

#include "constants.h"

   implicit none

   interface
      subroutine z4c_calculateADM(vars, lim)
         implicit none
         real, dimension(:, :, :, :), pointer :: vars
         integer, dimension(LOW:HIGH, MDIM) :: lim
      end subroutine z4c_calculateADM
   end interface

   interface
      subroutine z4c_calculateConstraintViolation(vars, lim, del)
         implicit none
         real, dimension(:, :, :, :), pointer :: vars
         integer, dimension(LOW:HIGH, MDIM) :: lim
         real, dimension(MDIM) :: del
      end subroutine z4c_calculateConstraintViolation
   end interface

end module z4c_interface
