!!****f* source/Grid/Grid_computeVarNorm
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
!!  Grid_computeVarNorm
!!
!! SYNOPSIS
!!  
!!  call Grid_computeVarNorm(integer(in)  :: level,
!!                           integer(in)  :: normType,
!!                           integer(in)  :: ivar,
!!                           real(out)    :: norm,
!!                           integer(in)  :: leafOnly)
!!
!! DESCRIPTION
!!
!!  Computes the L1 or L2 norm of the variable specified by ivar.  This
!!  can be done per-level, or on leaf or all nodes.  For multigrid, the
!!  L2 norm is used for convergence, but the L1 norm is incredibly useful
!!  for debugging purposes.
!!
!! ARGUMENTS
!!
!!  level     - If the norm is restricted to a given level; 0 is all
!!  normType - p in the Lp norm where choices of p are 1 or 2
!!  ivar      - the grid variable being normed; -1 for work
!!  norm      - the variable with which to return the norm
!!  leafOnly - if this isn't 0, compute the norm only on leaf nodes
!!
!! RESULT
!!
!!  The norm of ivar is in norm.
!!
!! NOTES
!!
!!  DEV: Currently only implemented for Paramesh4 and UG!
!!
!! EXAMPLE
!!  
!!  call gr_restrictTree()
!!  do i = 1, lrefine_max
!!    call Grid_computeVarNorm(i, 1, pdens, norm(i), 0)
!!  enddo
!!  do i = 1, lrefine_max
!!    if (norm(0) - norm(i) > 0.0000001) then
!!    call Driver_abort("restriction is highly nonconservatory!")
!!    endif
!!  enddo
!!
!!***

subroutine Grid_computeVarNorm (level, normType, ivar, norm, leafOnly)

  implicit none

  integer, intent(IN)  :: normType, level, ivar, leafOnly
  real, intent(OUT)    :: norm

  norm = 0.0
end subroutine Grid_computeVarNorm
