!!****if* source/Driver/DriverMain/Driver_envGetScalar
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
!!  Driver_envGetScalar
!!
!! SYNOPSIS
!!
!!  call Driver_envGetScalar(character(len=*)(in) :: name,
!!                           integer(out) :: value)
!!  call Driver_envGetScalarInt(character(len=*)(in) :: name,
!!                           integer(out) :: value)
!!
!! DESCRIPTION
!!
!!  Get an integer value from a named environment variable.
!!
!! ARGUMENTS
!!
!!   name : the name as a Fortran character variable
!!
!!   value : the integer value returned, or -1 if the
!!           variable did not exist in the environment,
!!           is not a recognized representation of
!!           an integer, or some other error occurred.
!!
!! NOTES
!!
!!   The specific implementation Driver_envGetScalarInt is
!!   made available under the generic name Driver_envGetScalar
!!   by a declaration in Driver_interface.
!!
!! SEE ALSO
!!
!!   Driver_interface
!!
!!***

subroutine Driver_envGetScalarInt(name, value)

  implicit none
#include "constants.h"

  character(len=*), intent(in)          :: name
  integer, intent(out)                  :: value

  character(len=MAX_STRING_LENGTH)  :: buff1
  integer :: status,vlen,res

  call get_environment_variable(name, buff1, length=vlen, status=status)

  if (vlen .LE. 0 .OR. status .NE. 0) then
     res = -1
  else
1    format(I0)
     read(buff1,1,iostat=status) res
     if (status .NE. 0) then
        res = -1
     end if
  end if

  value = res

end subroutine Driver_envGetScalarInt
