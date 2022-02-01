!!****if* source/flashUtilities/rng/mt_rng/ut_randomNumber
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
!!  ut_randomNumber
!!
!! SYNOPSIS
!!
!!  call ut_randomNumber(real, intent(OUT)  :: x)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   x : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

!! generates a random number on [0,1) with 53-bit resolution
  subroutine ut_randomNumber(x)
    implicit none
    real, intent(OUT) :: x
    integer,parameter :: realEightKind = selected_real_kind(15)
    real(kind=realEightKind)::y
    call ut_rand(y)
    x=y
  end subroutine ut_randomNumber
