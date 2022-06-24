!!****f* source/numericalTools/MoL/MoLMain/MoL_regrid
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
!!  NAME
!!
!!      MoL_regrid
!!
!!  SYNOPSIS
!!
!!      call MoL_regrid()
!!
!!  DESCRIPTION
!!
!!      Regrid internally-managed MoL memory  for intermediate stages
!!
!!  ARGUMENTS
!!
!!
!!  NOTES
!!
!!      Required information to regrid will be queried from Grid
!!
!!***
subroutine MoL_regrid()
    use ml_memInterface, only: ml_memAlloc, ml_memFree

    implicit none

    call ml_memFree
    call ml_memAlloc
end subroutine MoL_regrid
