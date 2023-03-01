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
!! @brief ml_memZero implementation

!> @ingroup MoLMemory
!!
!! @brief Implements ml_memZero
!!
!! @note This implementation is specific to UG and Paramesh
!!
!! @stubref{ml_memZero}
subroutine ml_memZero(dataStruct)
   use ml_memData, only: ml_scratch_data

#include "MoL.h"

   implicit none

   integer, intent(in) :: dataStruct

   if (dataStruct .ge. MOL_INITIAL) ml_scratch_data(:, :, :, :, :, dataStruct) = 0.0
end subroutine ml_memZero
