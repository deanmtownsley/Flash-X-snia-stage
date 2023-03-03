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
!! @brief ml_memReleaseActiveRHS implementation

!> @ingroup MoLMemory
!!
!! @brief Implements ml_memReleaseActiveRHS
!!
!! @stubref{ml_memReleaseActiveRHS}
subroutine ml_memReleaseActiveRHS(irhs)
   use ml_memData, only: ml_activeRHS

#include "MoL.h"

   implicit none

   integer, intent(in) :: irhs

   ml_activeRHS = MOL_INVALID
end subroutine ml_memReleaseActiveRHS
