!!****if* source/Grid/GridMain/Grid_getCellCoords
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
!!  Grid_getCellCoords
!!
!! SYNOPSIS
!!
!!  call Grid_getCellCoords(integer(IN)  :: axis,
!!                      integer(IN):: edge,
!!                      integer(IN):: level,
!!                      integer(IN):: lo(1:MDIM),
!!                      integer(IN):: hi(1:MDIM),
!!                      real(OUT)  :: coordinates(:))
!!
!!
!! DESCRIPTION
!!
!!    This subroutine is an accessor function that gets the coordinates of
!!    the cells in a given block.
!!    Coordinates are retrieved one axis at a time, 
!!    meaning you can get the i, j, _or_ k coordinates with one call.  
!!    If you want all the coordinates, all axes, you
!!    need to call Grid_getCellCoords 3 times, one for each axis.
!!    The code carries coordinates at cell centers as well as faces.
!!    It is possible to get coordinates for CENTER, only LEFT_EDGE,
!!    only RIGHT_EDGE or for all FACES along a dimension.
!!
!!
!!
!!
!! ARGUMENTS
!!            
!!   axis - specifies the integer index coordinates of the cells being retrieved.
!!          axis can have one of three different values, IAXIS, JAXIS or KAXIS 
!!          (defined in constants.h as 1,2 and 3)
!!
!!   edge - integer value with one of four values, 
!!          LEFT_EDGE, RIGHT_EDGE, CENTER or FACES
!!          The edge argument specifies what side of the zone to get, 
!!          the CENTER point, the LEFT_EDGE  or the RIGHT_EDGE of the zone.
!!          FACES gets the left and right face of each cell, but since 
!!          two adjacent cells have a common face, there are only N+1
!!          unique values if N is the number of cells.
!!
!!   level - refinement level.
!!           This is 1-based, i.e., the root level is numbered 1.
!!          
!!   lo   - Indices of a low point, in the level-wide integer index space.
!!          Only the component indicated by the axis argument is used.
!!          The coordinates returned pertain to cells in the range
!!          lo(axis) ... hi(axis).
!!
!!   hi   - Indices of a high point, in the level-wide integer index space.
!!          Only the component indicated by the axis argument is used.
!!          The coordinates returned pertain to cells in the range
!!          lo(axis) ... hi(axis).
!!
!!   coordinates - The array holding the data returning the coordinate values.
!!                 The array must be at large enough to hold the number of
!!                 coordinate values that are requested occording to lo and hi
!!                 arguments. That is, the following must be true:
!!           
!!          if edge = CENTER/LEFT_EDGE/RIGHT_EDGE then
!!                size(coordinates)  >=  hi(axis) - lo(axis) + 1 ;
!!
!!          If edge=FACES then
!!                size(coordinates)  >=  hi(axis) - lo(axis) + 2 .
!!
!!               
!!  EXAMPLE 
!!
!!  1. Getting cell centered values
!!
!!   #include "constants.h"
!!
!!      
!!      integer             :: coordSize
!!      integer,allocatable :: xCoord(:)
!!      integer             :: level,lo(3),hi(3)
!!      .....
!!      
!!      do while(itor%isValid())
!!          call itor%currentTile(tileDesc)
!!          level = tileDesc % level
!!          lo    = tileDesc % limits(LOW ,:)
!!          hi    = tileDesc % limits(HIGH,:)
!!
!!          coordSize = hi(IAXIS) - lo(IAXIS) + 1
!!          allocate(xCoord(coordSize)) !sized to be number of coords returned
!!
!!          call Grid_getCellCoords(IAXIS, CENTER, level, lo, hi, xCoord)
!!          .....
!!          deallocate(xCoord)
!!          call itor%next()
!!     end do    
!!
!!  2. Getting face values
!! 
!!   #include "constants.h"
!!
!!      
!!      integer             :: coordSize
!!      integer,allocatable :: xCoord(:)
!!      integer             :: level,lo(3),hi(3)
!!      
!!      do while(itor%isValid())
!!          call itor%currentTile(tileDesc)
!!          level = tileDesc % level
!!          lo    = tileDesc % limits(LOW ,:)
!!          hi    = tileDesc % limits(HIGH,:)
!!
!!          coordSize = hi(IAXIS) - lo(IAXIS) + 2
!!          allocate(xCoord(coordSize)) !sized to be number of coords returned
!!
!!          call Grid_getCellCoords(IAXIS, FACES, level, lo, hi, xCoord)
!!          .....
!!          deallocate(xCoord)
!!          call itor%next()
!!     end do    
!!
!!
!!  NOTES
!!   Variables that start with "gr_" are variables of Grid unit scope
!!   and are stored in the Fortran module Grid_data. Variables that are not
!!   starting with gr_ are local variables or arguments passed to the
!!   routine.
!!
!!***

#ifdef DEBUG
#define DEBUG_GRID
#endif

#include "constants.h"
#include "Simulation.h"

subroutine Grid_getCellCoords(axis, edge, level, lo, hi, coordinates)
  use Grid_interface,   ONLY : Grid_getDeltas
  use Grid_data,        ONLY : gr_globalDomain
  use Driver_interface, ONLY : Driver_abort

  implicit none

  integer, intent(in)  :: axis
  integer, intent(in)  :: edge
  integer, intent(in)  :: level
  integer, intent(in)  :: lo(1:MDIM)
  integer, intent(in)  :: hi(1:MDIM)
  real,    intent(out) :: coordinates(:)

  real    :: shift
  integer :: nElements
  real    :: deltas(1:MDIM)
  integer :: i

#ifdef DEBUG_GRID
  if((axis /= IAXIS) .AND. (axis /= JAXIS) .AND. (axis /= KAXIS)) then
     call Driver_abort("[Grid_getCellCoords] invalid axis, must be IAXIS, JAXIS or KAXIS ")
  end if
#endif
 
  ! Calculate number of cells
  nElements = hi(axis) - lo(axis) + 1
  if      (edge == FACES)  then
      shift = 2.0
      ! Get number of faces
      nElements = nElements + 1
  else if (edge == LEFT_EDGE) then
      shift = 2.0
  else if (edge == CENTER) then
      shift = 1.5
  else if (edge == RIGHT_EDGE) then
      shift = 1.0
  else
      call Driver_abort('[Grid_getCellCoords] Invalid edge')
  end if

  if (SIZE(coordinates) < nElements) then
      call Driver_abort("[Grid_getCellCoords] coordinates is too small")
  end if

  call Grid_getDeltas(level, deltas)

  associate (x0 => gr_globalDomain(LOW, axis), &
             dx => deltas(axis))
      do i = 1, nElements
          coordinates(i) = x0 + (lo(axis) + i - shift) * dx
      end do
  end associate
end subroutine Grid_getCellCoords

