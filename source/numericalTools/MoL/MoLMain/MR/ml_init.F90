!!****if* source/numericalTools/MoL/MoLMain/MR/ml_init
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
!!      ml_init
!!
!!  SYNOPSIS
!!
!!      call ml_init()
!!
!!  DESCRIPTION
!!
!!      Initialize a method of lines unit implementation
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_init()
    use MoL_data, only: MoL_nscratch
    use mr_data
    use erk, only: erk_init
    use gark, only: gark_init

    use RuntimeParameters_interface, only: RuntimeParameters_get

#include "Simulation.h"
#include "constants.h"
#include "MoL.h"

    implicit none

    character(len=MAX_STRING_LENGTH) :: slowMethod_str, fastMethod_str

    integer :: i

    call RuntimeParameters_get("mr_slowMethod", slowMethod_str)
    call RuntimeParameters_get("mr_fastMethod", fastMethod_str)

    mr_slowMethod = trim(slowMethod_str)
    mr_fastMethod = trim(fastMethod_str)

    call RuntimeParameters_get("mr_nsubcycle", mr_nsubcycle)

    call gark_init
    call erk_init

    ! Setup RHS indexing
    allocate(FE(mr_nstages_slow))
    allocate(FI(mr_nstages_slow))
    allocate(FF(mr_nstages_fast))

    FE = MOL_INVALID; FI = MOL_INVALID; FF = MOL_INVALID

    ! Uses MOL_RHS as first index
    do i = 1, mr_nstages_slow, 2
        FE(i) = MOL_RHS + i - 1
        FI(i) = MOL_RHS + i
    end do ! i

    FAST_INITIAL = FI(mr_nstages_slow-1) + 1

    do i = 1, mr_nstages_fast
        FF(i) = FAST_INITIAL + i
    end do ! i

    MoL_nscratch = mr_nstages_slow + mr_nstages_fast
end subroutine ml_init
