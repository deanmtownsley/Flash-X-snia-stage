!!****if* source/Grid/GridMain/paramesh/Grid_getDeltas
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
!!  Grid_getDeltas
!!
!! SYNOPSIS
!!
!!  Grid_getDeltas(integer(IN) :: blockId,
!!                 real(OUT)   :: del(MDIM))
!!  
!! DESCRIPTION 
!!  
!!  Gets the grid spacing dx/dy/dz for a given blockId on the Grid.
!!  dx is the size of one cell in the x direction of a block.
!!  
!!  
!! ARGUMENTS 
!!
!!  blockId - local block number
!!  del - array of size MDIM returned holding the dx, dy, and dz values
!!
!!***

subroutine Grid_getDeltas(lev,del)
  use Grid_data, ONLY: gr_delta
  implicit none

#include "constants.h"
  
  integer, intent(IN)   :: lev
  real, dimension(MDIM), intent(out) :: del

  if (lev < 1) then
     print*,'**lev**',lev,'!'
  end if
  del(1:MDIM) = gr_delta(1:MDIM,lev)
  return
end subroutine Grid_getDeltas
