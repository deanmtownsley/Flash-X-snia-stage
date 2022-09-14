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
!! A Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface
!! information.
!!
!! This subroutine is presently a noop and as a result, it is one factor that
!! limits the Milhoja Grid unit implementation to pseudo-UG.
!!
!! @todo Implement
subroutine Grid_updateRefinement(nstep, time, gridChanged)
    implicit none

    integer, intent(IN)            :: nstep
    real,    intent(IN)            :: time
    logical, intent(OUT), OPTIONAL :: gridChanged

    RETURN
end subroutine Grid_updateRefinement

