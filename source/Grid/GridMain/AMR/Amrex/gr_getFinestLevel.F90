subroutine gr_getFinestLevel(level)
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
    use amrex_amrcore_module, ONLY : amrex_get_finest_level

    implicit none

    integer, intent(OUT) :: level

    ! AMReX uses 0-based level indexing/FLASH uses 1-based
    level = amrex_get_finest_level() + 1
end subroutine gr_getFinestLevel

