!***************************************************************************************************
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
! xnet_types.f90 10/18/17
! This file contains the module defining different fortran types.
!***************************************************************************************************

Module xnet_types
  !-------------------------------------------------------------------------------------------------
  ! These are explicitly defined type parameters for floating-point (real) number precision, with
  ! double precision declared as Real(dp) instead of Real(8).
  !-------------------------------------------------------------------------------------------------
  Use, Intrinsic :: iso_fortran_env, Only: int32, int64, real32, real64, real128
  Implicit None
  Integer, Parameter :: i4 = int64   ! 32-bit integer
  Integer, Parameter :: i8 = int64   ! 64-bit integer
  Integer, Parameter :: sp = real32  ! Single precision
  Integer, Parameter :: dp = real64  ! Double precision
  Integer, Parameter :: qp = real128 ! Quad percision
End Module xnet_types