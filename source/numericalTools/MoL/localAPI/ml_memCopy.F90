!!****f* source/numericalTools/MoL/localAPI/ml_memCopy
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
!!  NAME
!!
!!      ml_memCopy
!!
!!  SYNOPSIS
!!
!!      call ml_memCopy(integer, intent(in) :: dst
!!                      integer, intent(in) :: src)
!!
!!  DESCRIPTION
!!
!!      Copy the variables evolved by MoL from provided source to destination
!!      memory locations.  These memory locations either be UNK or MoL-specific
!!      scratch memory.
!!
!!      Valid locations include (defined in MoL.h):
!!          - MOL_EVOLVED : Evolved variables in UNK
!!          - MOL_INITIAL : Copy of the evolved variables at the start of a timestep
!!          - MOL_RHS     : The currently-being-calculated RHS terms
!!          - other       : Each integrator may specify some additional number of
!!                          of scratch-memory for intermediate stages/RHS terms
!!
!!  ARGUMENTS
!!
!!      dst : Index of the destination location to copy to
!!      src : Index of the source location to copy from
!!
!!***
subroutine ml_memCopy(dst, src)
    implicit none

    integer, intent(in) :: dst, src

    return
end subroutine ml_memCopy
