!!****f* source/Grid/Grid_formatNonRep
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!
!!  Grid_formatNonRep
!!
!!
!! SYNOPSIS
!!
!!  Grid_formatNonRep(integer(IN) :: nonrep,
!!                    integer(IN) :: idx,
!!                    character(out) :: str(*))
!!
!!
!! DESCRIPTION
!!
!!  Given a nonreplicated variable array id and index into that array, returns a string name suitable for IO
!!
!!
!! ARGUMENTS
!!  
!!   nonrep - array id
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
