!!****if* source/RuntimeParameters/RuntimeParametersMain/rp_getArgument
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
!!    rp_getArgument
!!
!! SYNOPSIS
!!    call rp_getArgument(integer(in) :: pos,
!!                        character(len=MAX_STRING_LENGTH)(out) :: arg)
!!
!! DESCRIPTION
!!
!!     Plug-in routine to return a command-line argument.
!!     Specify the position of the argument on the command
!!     line (beginning with 1) and a string variable to
!!     receive the argument value.
!!
!!     This is a wrapper which can be replaced when running
!!     on systems that do not have this extension to the standard.
!!
!! ARGUMENTS
!!
!!     pos :   integer argument position
!!     arg :  the returned command-line argument
!!
!!
!!
!!***
 
 
subroutine rp_getArgument (pos, arg)

#ifdef NAGF95
  use f90_unix_env, ONLY: getarg
#endif

  implicit none
 
#include "constants.h"

  integer, intent(in) ::          pos
  character(len=MAX_STRING_LENGTH), intent(out) :: arg
 
  call getarg (pos, arg)
 
  return
end subroutine rp_getArgument
