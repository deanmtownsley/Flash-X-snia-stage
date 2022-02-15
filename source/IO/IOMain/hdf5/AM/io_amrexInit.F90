!!***if* source/IO/IOMain/hdf5/AM/io_amrexInit
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
!!  io_amrexInit
!!
!! SYNOPSIS
!!
!!  io_amrexInit()
!!
!! DESCRIPTION
!!
!!   Perform initialization for source/IO/IOMain/hdf5/AM 
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! NOTES
!!
!!  Initialize variables need by IO operation involving AMReX
!
!!
!!
!!***

subroutine io_amrexInit()

    use IO_data, ONLY: io_plotfileGridQuantityDP 
    use RuntimeParameters_interface, ONLY : RuntimeParameters_get
    use io_amrexData

    call RuntimeParameters_get('io_plotFileAmrexFormat', io_plotFileAmrexFormat)

    if (.not. io_plotFileAmrexFormat) io_plotfileGridQuantityDP = .TRUE.

end subroutine io_amrexInit
