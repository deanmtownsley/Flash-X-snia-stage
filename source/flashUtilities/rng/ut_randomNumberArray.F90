!!****if* source/flashUtilities/rng/ut_randomNumberArray
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
!!  ut_randomNumberArray
!!
!! SYNOPSIS
!!
!!  call ut_randomNumberArray(real, intent(OUT)  :: r1)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   r1 : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

  subroutine ut_randomNumberArray(r1)
    implicit none
    real, intent(OUT) :: r1(:)
    call random_number(r1)
  end subroutine ut_randomNumberArray
