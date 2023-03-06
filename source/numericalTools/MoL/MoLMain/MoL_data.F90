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
!! @brief Main data module for MoL

!> @ingroup MoLMain
!! Variables storing data required in MoL top-level implmentations
module MoL_data

   implicit none

   !> @name Scratch Memory Totals
   !! @{

   !> Number of memory-levels (e.g. intermediate stages) required by an integrator
   integer, save :: ml_nscratch

   !> Number of memory-levels (e.g. intermediate stages) required by an
   !! integrator, plus those provided by the base implementation
   integer, save :: ml_nscratch_total
   !> @}

   !> @name Verbosity Settings
   !! @{

   !> Verbosity level for messaging as defined in @ref MoL.h
   integer, save :: ml_verbosity

   !> Whether to abort when a warning is issued
   logical, save :: ml_abortOnWarn

   !> Current MPI rank (only outputs messages on the master process)
   integer, save :: ml_mpiRank
   !> @}

end module MoL_data
