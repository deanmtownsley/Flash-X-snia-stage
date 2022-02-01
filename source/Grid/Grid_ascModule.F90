!!****f* source/Grid/Grid_ascModule
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
!!  Grid_ascModule
!!
!! SYNOPSIS
!!
!!  use Grid_ascModule
!!
!! DESCRIPTION 
!!  
!!  Module with data and some accessor routines for allocatable scratches.
!!
!!   
!!***

Module Grid_ascModule

  implicit none


contains
  subroutine Grid_ascStart()
  end subroutine Grid_ascStart

  subroutine Grid_ascAllocMem(gds,var1,nvars, nGuardCtr,nGuardFaceN,nGuardFaceT, &
                                leafBlocks, arrayRank, highSize)

    integer, intent(IN) :: gds,var1,nvars
    integer, OPTIONAL, intent(IN) :: nGuardCtr,nGuardFaceN,nGuardFaceT
    integer, OPTIONAL, intent(IN), dimension(:) :: leafBlocks
    integer, OPTIONAL, intent(IN) :: arrayRank
    integer, OPTIONAL, intent(IN) :: highSize
  end subroutine Grid_ascAllocMem

  subroutine Grid_ascDeallocMem(gds, arrayRank)

    integer, OPTIONAL, intent(IN) :: gds
    integer, OPTIONAL, intent(IN) :: arrayRank

  end subroutine Grid_ascDeallocMem

end Module Grid_ascModule
