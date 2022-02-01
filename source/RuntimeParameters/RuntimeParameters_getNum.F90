!!****f* source/RuntimeParameters/RuntimeParameters_getNum
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
!!  RuntimeParameters_getNum
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_getNum(integer(out) :: nparms)
!!                      
!!
!! DESCRIPTION
!!
!!  Returns the number of parameters of a given type.  The
!!  routines below are implemented; they are not overloaded,
!!  so call them directly.
!!
!!  RuntimeParameter_getNumReal(integer(out) :: nparms)
!!  RuntimeParameter_getNumInt(integer(out) :: nparms)
!!  RuntimeParameter_getNumStr(integer(out) :: nparms)
!!  RuntimeParameter_getNumLog(integer(out) :: nparms)
!!
!! ARGUMENTS
!!
!!  nparms:     number of parameters
!!
!!
!!
!!***


   
subroutine RuntimeParameters_getNumReal (nparms)
  implicit none
  integer, intent(out)                :: nparms
end subroutine RuntimeParameters_getNumReal


   
subroutine RuntimeParameters_getNumInt (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumInt


   
subroutine RuntimeParameters_getNumStr (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumStr


   
subroutine RuntimeParameters_getNumLog (nparms)
  implicit none
  integer, intent(out)               :: nparms
end subroutine RuntimeParameters_getNumLog



