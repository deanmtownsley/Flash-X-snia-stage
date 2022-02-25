!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_getIntCoords.F90
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
!!   gr_getIntCoords.F90
!!
!! SYNOPSIS
!!
!!   gr_getIntCoords(integer(IN) :: lblock,
!!                   integer(IN) :: lcoord(MDIM) )
!!
!! DESCRIPTION
!!   Returns integer coordinates (0-based) for any local block.
!!
!! ARGUMENTS
!!   lblock - local block number
!!   lcoord - stores output
!!
!! NOTES
!!
!!   This subroutine is used in the Paramesh4 Grid implementation with
!!   Bittree. It should also be usable, and may be useful, when using
!!   the Paramesh4 Grid implementation without Bittree, and does not
!!   rely on any data items that are specific to Bittree; therefore
!!   it was moved to a more general directory level. - KW 2022-02-17
!!***
#include "constants.h"
      subroutine gr_getIntCoords(lblock,lcoord)

      use paramesh_dimensions, only: ndim
      use tree, only: bsize,coord,grid_xmin,grid_ymin,grid_zmin 

      implicit none

      integer, intent(in):: lblock
      integer, dimension(MDIM),intent(out):: lcoord


      lcoord(:) = 0
      if(ndim >= 1) lcoord(1) = int((coord(1,lblock)-grid_xmin)&
                                    /bsize(1,lblock))
      if(ndim >= 2) lcoord(2) = int((coord(2,lblock)-grid_ymin)&
                                    /bsize(2,lblock))
      if(ndim >= 3) lcoord(3) = int((coord(3,lblock)-grid_zmin)&
                                    /bsize(3,lblock))

      return

      end subroutine
