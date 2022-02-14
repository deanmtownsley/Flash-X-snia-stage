!!****if* source/flashUtilities/general/ut_constRealZeroFn1
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
!!  ut_constRealZeroFn1
!!  ut_constRealOneFn1
!!
!! SYNOPSIS
!!
!!  real a0 = ut_constRealZeroFn1 (real(in)  :: x)
!!
!!  real a1 = ut_constRealOneFn1  (real(in)  :: x)
!!
!! DESCRIPTION
!!
!!  Functions of one real argument that always returns
!!  the same value.
!!
!! ARGUMENTS
!!
!!  x         : ignored argument
!!
!! NOTES
!!
!!***

real function ut_constRealZeroFn1 (x)
  implicit none
  real,    intent (in)  :: x

  ut_constRealZeroFn1 = 0.0

end function ut_constRealZeroFn1

real function ut_constRealOneFn1 (x)
  implicit none
  real,    intent (in)  :: x

  ut_constRealOneFn1 = 1.0

end function ut_constRealOneFn1
