!!****if* source/IO/IOMain/io_getNumScalars
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
!!  io_getNumScalars
!!
!! SYNOPSIS
!!
!!  io_getNumScalars(integer(out) :: nscalars)
!!                      
!!
!! DESCRIPTION
!!
!! This function returns the number of scalars of a certain type
!! (real/int/str/log) from a 
!! linked list implemented under
!! the hood.  This function is not overladed.  Call the specific routine
!! directly.  
!!
!! ARGUMENTS
!!
!! nscalars:     number of scalars
!!
!!
!! EXAMPLE
!!
!!  call io_getNumScalarsReal(numRealScalars)
!!
!!***


   
subroutine io_getNumScalarsReal (nscalars)

  use IO_data, ONLY : io_scalar

implicit none
  integer, intent(inout)                      :: nscalars

  call NameValueLL_getNumReal(io_scalar, nscalars)

  return
  
end subroutine io_getNumScalarsReal


   
subroutine io_getNumScalarsInt (nscalars)

  use IO_data, ONLY : io_scalar

implicit none
  integer, intent(inout)                :: nscalars

  call NameValueLL_getNumInt(io_scalar, nscalars)

  return
  
end subroutine io_getNumScalarsInt


   
subroutine io_getNumScalarsStr (nscalars)

  use IO_data, ONLY : io_scalar

implicit none
  integer, intent(inout)          :: nscalars

  call NameValueLL_getNumStr(io_scalar, nscalars)

  return
  
end subroutine io_getNumScalarsStr


   
subroutine io_getNumScalarsLog (nscalars)

  use IO_data, ONLY : io_scalar

implicit none
  integer, intent(inout)                      :: nscalars

  call NameValueLL_getNumLog(io_scalar, nscalars)

  return
  
end subroutine io_getNumScalarsLog



