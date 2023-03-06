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
!! @brief MoL evolved variable tracking

!> @ingroup MoLMain
!! Provides tracking of registered evolved variables
module ml_variables

#include "Simulation.h"
#include "MoL.h"

   implicit none

   !> Associates an evolved variable index in UNK with an RHS index
   !! in MoL scratch memory
   type :: ml_variable_t
      character(len=:), allocatable :: name !< Name of the variable
      integer :: evolIndex                  !< Index of the variable in UNK
      integer :: rhsIndex                   !< Index of the variable in MoL's scratch memory
   end type ml_variable_t

   !> @name Variable Tracking
   !! @{

   !> Registered evolved variables
   type(ml_variable_t), allocatable, save :: ml_vars(:)

   !> Total number of registered evolved variables
   !! @hideinitializer
   integer, save :: ml_nvars = 0
   !> @}

   !> @name Indexing Mapping
   !! @{
   !! Maps variable indices between UNK and MoL's scratch memory
   integer, allocatable, save :: ml_unk_mask(:), ml_scratch_mask(:)
   integer, save :: ml_unk_to_scratch(NUNK_VARS) = MOL_INVALID
   !> @}

end module ml_variables
