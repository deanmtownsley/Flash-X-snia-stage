!!****f* source/Grid/Grid_formatNonRep
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
!!  Grid_formatNonRep
!!
!!
!! SYNOPSIS
!!
!!  call Grid_formatNonRep(integer(IN) :: nonrep,
!!                    integer(IN) :: idx,
!!                    character(out) :: str(*))
!!
!!
!! DESCRIPTION
!!
!!  Given the ID of a nonreplicated variable array and the index into that array,
!!  returns a string name suitable for IO.
!!
!!
!! ARGUMENTS
!!  
!!   nonrep - array ID
!!   idx - index into array
!!   str - receives string name
!!
!!***
subroutine Grid_formatNonRep(nonrep, idx, str)
    implicit none
    integer, intent(in) :: nonrep, idx
    character(*), intent(out) :: str
    str = ''
end subroutine
