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

!> @ingroup Orchestration
!! @anchor Orchestration_finalize_stub
!!
!! @brief Finalize the Orchestration unit
!!
!! @details
!! This routine should be called at most once during program execution.  It should
!! only be called if Orchestration_init has been called, in which case it must be
!! called.  If this routine is called, it must be called before finalizing the Grid
!! unit.   Called, called, called...
!!
!! It is intended that only the Driver unit call this routine.
subroutine Orchestration_finalize()
    implicit none
    RETURN
end subroutine Orchestration_finalize

