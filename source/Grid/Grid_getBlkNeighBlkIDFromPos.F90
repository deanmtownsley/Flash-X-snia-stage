!!****f* source/Grid/Grid_getBlkNeighBlkIDFromPos
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
!!  Grid_getBlkNeighBlkIDFromPos
!!
!! SYNOPSIS
!!
!!  call Grid_getBlkNeighBlkIDFromPos(
!!                            type(Grid_tile_t)(IN) :: blockDesc,
!!                            real(IN)    :: pos(MDIM), 
!!                            integer(IN) :: neghDir(MDIM),
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID)
!!
!!  call Grid_getBlkIDFromPos(type(Grid_tile_t)(IN) :: blockDesc,
!!                            real(IN)    :: pos(MDIM), 
!!                            integer(IN) :: neghDir(MDIM),
!!                            integer(OUT) :: ansBlockID,
!!                            integer(OUT) :: ansProcID)
!!
!! DESCRIPTION
!!
!!   Given the physical coordinates of a point outside the current
!!   block (described by blockDesc) but in its neighborhood, and
!!   the direction in which the neighboring block lies, this routine
!!   finds the processor number and blockID  within that processor
!!   for the neighboring block that contains the point.
!!
!! ARGUMENTS
!!
!!   blockDesc : describes block in current processor
!!
!!   pos :     coordinates of the point of interest
!!
!!   neghDir : the location of the neighbor with respect to the
!!             current block, in other words specification on which
!!             face/edge/point is common between the current block and
!!             neighbor of interest. For example
!!             negh(1:MDIM)=LEFT_EDGE indicates that the lowest
!!             left hand corner of the current block is the same as
!!             the highest right end corner of the neighbor. Similarly
!!             negh(IAXIS)=LEFT_EDGE, negh(JAXIS:KAXIS) =
!!             CENTER implies that the left face of current block is
!!             common with the right face of the neighbor.
!!             This array is in the form returned by Grid_outsideBoundBox.
!!
!!   ansBlockID : identity of the neighbor, the first number 
!!                is the blocknumber within the processor
!!   ansProcID :  identity of the neighbor, the second number is the processor
!!                number where the neighbor is located
!!
!! NOTES
!!
!!   Currently at most implemented for PARAMESH Grid implementations.
!!
!!   The specific subroutine Grid_getBlkNeighBlkIDFromPos is also available
!!   under the generic name Grid_getBlkIDFromPos.
!!
!!   !DEV: There is no proper implementation, Paramesh4 implementation will abort if called.
!!
!! SEE ALSO
!!   Grid_outsideBoundBox
!!***

subroutine Grid_getBlkNeighBlkIDFromPos(blockDesc,pos,neghDir,ansBlockID,ansProcID)
#include "constants.h"
 
  use Grid_tile, ONLY : Grid_tile_t
  
  implicit none
  type(Grid_tile_t), intent(IN) :: blockDesc
  real,dimension(MDIM),intent(IN) :: pos
  integer,dimension(MDIM),intent(IN) :: neghDir
  integer,intent(OUT) :: ansBlockID, ansProcID

  
  ansProcID=NONEXISTENT
  ansBlockID=NONEXISTENT

end subroutine Grid_getBlkNeighBlkIDFromPos
