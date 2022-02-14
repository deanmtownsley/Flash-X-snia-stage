!!****f* source/RuntimeParameters/RuntimeParameters_read
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
!!  RuntimeParameters_read
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_read( character(in), (len=MAX_STRING_LEN) :: parmfile)
!!
!!
!! DESCRIPTION
!!
!!
!!         Parses the parameter file parmfile, which contains
!!         job-dependent parameter definitions.  Syntax of
!!         parameter file lines is:
!!
!!                       # comment               rest of line is a comment
!!                       variable = value        set variable to value
!!                       strvar   = "word"     set string variable strvar
!!                                                to word
!!
!!         Some error checking is performed:  if a variable is
!!         unrecognized, it is ignored.  Syntax errors are signaled
!!         along with the offending lines.  However, type mismatch
!!         errors force termination.  Also, no checking is done to
!!         determine whether all of the variables have received some
!!         value.  In case of repeated definitions of the same
!!         variable, the last definition overrides the others.
!!
!! ARGUMENTS
!!
!!   parmfile :       the name of the parameter file to read
!! 
!! NOTES
!!
!!   This routine is called during FLASH initialization to read the flash.par file.
!!   In general it would not be used by users.  Instead, use the restart = .true.
!!   capability within a flash.par to restart a run from a checkpoint file,
!!   using a different flash.par configuration.
!!
!!
!!
!!***



subroutine RuntimeParameters_read (parmfile)

#include "constants.h"

  implicit none
  

  character(len=MAX_STRING_LENGTH), intent(in)    :: parmfile


  return
end subroutine RuntimeParameters_read



