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
!! @brief Main data module for MoLMemory

!> @ingroup MoLMemory
!! MoL scratch memory structures and trackers
module ml_memData

#include "MoL.h"

   implicit none

   !> Scratch memory for intermediate stage storage, indexed as (var,i,j,k,block,dataStruct)
   real, dimension(:, :, :, :, :, :), allocatable, target, save :: ml_scratch_data

   !> Active RHS state to use when MOL_RHS is requested
   !! @hideinitializer
   integer, save :: ml_activeRHS = MOL_INVALID

end module ml_memData
