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
!! @brief Eos_init stub

!> @ingroup physics_Eos
!!
!! @brief Initialize the Eos unit
!!
!! @details
!! @anchor Eos_init_stub
!!
!! This procedure must be called once and only once during startup. Tasks
!! performed by this procedure will include initializationi of various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants f A call to this procedure must be matched by
!! a call to @ref Eos_finalize


subroutine Eos_init()
    
    implicit none
    ! stub for eos initialization.  A non-stub implementation of his routine
    ! will be supplied, if required,  by the main subunit of Eos, normally
    ! located under Eos/EosMain.

    return
end subroutine Eos_init


