!!****if* source/Grid/GridMain/AMR/Amrex/Grid_finalize
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
!!
!! NAME
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  call Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory that has been allocated in the Grid Unit
!!
!!***


subroutine Grid_finalize()
    use gr_bcInterface,    ONLY : gr_bcFinalize
    use gr_ptInterface,    ONLY : gr_ptFinalize
    use gr_amrexInterface, ONLY : gr_amrexFinalize

    implicit none

    call gr_solversFinalize()
    call gr_ptFinalize()
    call gr_bcFinalize()
    call gr_amrexFinalize()
end subroutine Grid_finalize
