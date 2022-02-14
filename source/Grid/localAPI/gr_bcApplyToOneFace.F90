!!****if* source/Grid/localAPI/gr_bcApplyToOneFace
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
!!  gr_bcApplyToOneFace
!!
!! SYNOPSIS
!!
!!  gr_bcApplyToOneFace(integer(IN)          :: axis,
!!                      integer(IN)          :: bcType,
!!                      integer(IN)          :: gridDataStruct,
!!                      integer(IN)          :: varCount,
!!                      integer(IN)          :: regionType(MDIM)
!!                      Grid_tile_t(IN) :: blockDesc,
!!                      integer(IN)          :: idest)
!!  
!! DESCRIPTION 
!!
!!  
!! 
!! ARGUMENTS
!!  
!!    axis           - the direction for applying BC, one of IAXIS, JAXIS, or KAXIS
!!    bcType         - the type of boundary condition
!!    gridDataStruct - In PM3 and PM4 it can have values (CENTER,FACEX,FACEY,FACEZ,WORK),
!!                     in UG (CENTER,FACEX,FACEY,FACEZ), and in PM2 (CENTER,WORK).
!!    varCount       - the number of variable in the data structure specified in gridDataStruct
!!    regionType     - The part of the block that is involved in the boundary condition. This integer
!!                     array can have values (LEFT_EDGE, CENTER, RIGHT_EDGE, WHOLE_VEC and NO_VEC)
!!                     for each of the three dimensions of the physical data.
!!                     LEFT_EDGE implies guard cells along lower face. If this value is specified
!!                     for the axis along which the BC are being applies, that means that we 
!!                     are applying BC to the lowerface, if it is one of the other axes, then
!!                     the BC are being applied to one of the corners. Same is true of RIGHT_EDGE
!!                     except that implies the upperface. CENTER, WHOLE_VEC and NO_VEC values
!!                     are valid only for dimensions other than the one specified in "axis". 
!!                     CENTER implies only the interior cells, WHOLE_VEC implies all cells in 
!!                     the block and NO_VEC implies that the correspoding dimension is not
!!                     a part of the region. Normally this value is most likely to be used
!!                     along KAXIS in a 2D problems, and JAXIS and KAXIS in a 1D problem
!!    blockDesc      - Derived type that encapsulates metadata that uniquely
!!                     characterizes local block to be operated on
!!    idest          - this is useful in paramesh, when the boundary conditions
!!                    are being applied to the workspace array WORK
!!
!! NOTES
!!  A specific direction is required in axis - no ALLDIR at this point.
!!
!!***

subroutine gr_bcApplyToOneFace(axis,bcType,gridDataStruct,varCount,&
     regionType,blockDesc,idest)
  use Grid_tile, ONLY : Grid_tile_t

  implicit none
#include "constants.h"
  
  integer, intent(in) :: axis,bcType,gridDataStruct,varCount,idest
  integer,dimension(MDIM),intent(IN) :: regionType
  type(Grid_tile_t), intent(IN) :: blockDesc

end subroutine gr_bcApplyToOneFace
