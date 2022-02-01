!!****f* source/Grid/GridMain/AMR/Amrex/Grid_restrictAllLevels
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
!!  Grid_restrictAllLevels
!!
!! SYNOPSIS
!!  call Grid_restrictAllLevels()
!!  
!! DESCRIPTION 
!!  Restricts the grid-managed, cell-centered data from leaf blocks down
!!  to all refinement levels. Normally FLASH only evolves on the leaf blocks,
!!  calling this routine makes all levels have valid data.  This is mostly for
!!  visualization purposes to be able to look at different levels of resolution.
!!
!! ARGUMENTS
!!  None
!!
!!***

#include "constants.h"

subroutine Grid_restrictAllLevels()
    use Grid_data,         ONLY : gr_convertToConsvdInMeshInterp, &
                                  gr_convertToConsvdForMeshCalls
    use gr_amrexInterface, ONLY : gr_restrictAllLevels

    implicit none

    logical :: needConversion

    needConversion = gr_convertToConsvdInMeshInterp .OR. gr_convertToConsvdInMeshInterp

    call gr_restrictAllLevels(CENTER, &
                              convertPtoC=needConversion, &
                              convertCtoP=needConversion)
end subroutine Grid_restrictAllLevels

