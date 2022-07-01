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
!! This subroutine has not yet been implemented and aborts if called.
!!
!! @todo Implement
subroutine gr_tile2txtytz(tileID, tx, ty, tz)
    use Driver_interface, ONLY : Driver_abort

    implicit none

    integer, intent(IN)  :: tileID
    integer, intent(OUT) :: tx, ty, tz

    tx = -1
    ty = -1
    tz = -1
    CALL Driver_abort("[gr_tile2txtytz] Not implemented yet")
end subroutine gr_tile2txtytz

