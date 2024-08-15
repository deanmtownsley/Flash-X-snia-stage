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
!! @anchor Orchestration_init_stub
!!
!! @brief Initialize the Orchestration unit
!!
!! @details
!! This routine should be called at most once during program execution and should
!! only be called if the Grid unit has already been initialized.  Upon terminating,
!! this routine leaves the Orchestration unit ready for full use.  If this routine
!! is called, then calling code must eventually call Orchestration_finalize.
!!
!! It is intended that only the Driver unit call this routine.
subroutine Orchestration_init()
    implicit none
    RETURN
end subroutine Orchestration_init

