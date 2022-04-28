!!****if* source/Grid/GridMain/paramesh/bittree/source/gr_getNeighIntCoords.F90
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
!!   gr_getNeighIntCoords.F90
!!
!! SYNOPSIS
!!
!!   gr_getNeighIntCoords(integer(IN)    :: lblock,
!!                        integer(IN)(3) :: gCell,
!!                        integer(OUT)   :: neighCoord(MDIM),
!!                        logical(in)    :: force (optional) )
!!
!! DESCRIPTION
!!   Returns integer coordinates (0-based) for neighbor of lblock 
!!   in a given direction. If periodic BCs, applies modulo to
!!   appropriate coordinate, allowing neighCoord to loop to other side.
!!
!! ARGUMENTS
!!   lblock - local number of block
!!   gCell - direction of desired neighbor (range [-1:1])
!!   neighCoord - integer coordinates of neighbor block
!!   force - logical that determines whether to return out-of-bounds case
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
      subroutine gr_getNeighIntCoords(lblock,gCell,neighCoord,force)

      use paramesh_dimensions, only: ndim
      use Grid_data, only: gr_domainBC
      use tree, only : grid_xmax,grid_xmin,grid_ymax, &
                       grid_ymin,grid_zmax,grid_zmin, &
                       bsize

      use Driver_interface,only: Driver_abort

      implicit none

      integer, intent(in):: lblock
      integer, dimension(MDIM),intent(in) :: gCell
      integer, dimension(MDIM),intent(out):: neighCoord
      logical, optional,intent(in) :: force
      
      
      integer, dimension(MDIM) :: lcoord
      integer :: i,maxcoord

      call gr_getIntCoords(lblock,lcoord)

!-----Calculate integer coordinates of neighbor in direction
      neighCoord = lcoord + gCell

!-----Make sure not out-of-bounds. If periodic BCs, apply modulo
      do i = 1,ndim
        if (i.eq.1) maxcoord = (grid_xmax-grid_xmin)/bsize(1,lblock) 
        if (i.eq.2) maxcoord = (grid_ymax-grid_ymin)/bsize(2,lblock) 
        if (i.eq.3) maxcoord = (grid_zmax-grid_zmin)/bsize(3,lblock) 

        if (neighCoord(i).lt.0 ) then
          if ( gr_domainBC(1,i).eq.PERIODIC ) then
            neighCoord(i) = modulo(neighCoord(i),maxcoord)
          else
            if (present(force)) then
              if (force) then
                neighCoord(i) = gr_domainBC(1,i)
              else
                call Driver_abort('Error in gr_getNeighIntCoords, &
                 &trying to compute out-of-bounds integer coordinates.')
              end if
            else
              call Driver_abort('Error in gr_getNeighIntCoords, &
                 &trying to compute out-of-bounds integer coordinates.')
            end if
          end if
        end if
        
        if (neighCoord(i).ge.maxcoord) then
          if ( gr_domainBC(2,i).eq.PERIODIC ) then
            neighCoord(i) = modulo(neighCoord(i),maxcoord)
          else
            if (present(force)) then
              if (force) then
                neighCoord(i) = gr_domainBC(2,i)
              else
                call Driver_abort('Error in gr_getNeighIntCoords, &
                 &trying to compute out-of-bounds integer coordinates.')
              end if
            else 
              call Driver_abort('Error in gr_getNeighIntCoords, &
                 &trying to compute out-of-bounds integer coordinates.')
            end if
          end if
        end if
      
      end do !i=1,3

      return

      end subroutine
