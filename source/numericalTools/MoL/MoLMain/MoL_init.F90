!!****if* source/numericalTools/MoL/MoLMain/MoL_init
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
!!      MoL_init
!!
!!  SYNOPSIS
!!
!!      call MoL_init()
!!
!!  DESCRIPTION
!!
!!      Initialize the method of lines unit
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine MoL_init()
    use ml_interface, only: ml_init
    use MoL_data

    implicit none

    MoL_nscratch = 0
    MoL_nvars    = 0

    call ml_init()

    ! +2 for MOL_INITIAL & MOL_RHS
    MoL_nscratch_total = 2 + MoL_nscratch
end subroutine MoL_init
