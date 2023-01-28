#include "constants.h"
#include "Simulation.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! Set the initial conditions in interiors and guardcells to simple values
!! that are used for checking correct read access to data via iterator and
!! Grid_tile_t.
!!
!! @todo Use Klaus's header file so that we can declare the pointer as
!!       intent(IN)
!!
!! @param initData  Cell-centered physical data array associated with the
!!                  given tile.  Flash-X only requires that the initial
!!                  conditions be set on the interior.  However, we set in the 
!!                  guardcells as well for enhanced testing.
!! @param tileDesc  The tile whose data is to be set.

!!REORDER(4): initData
subroutine Simulation_initBlock(initData, tileDesc)
    use Grid_tile, ONLY : Grid_tile_t 

    implicit none
    
    real,                         pointer :: initData(:, :, :, :)
    type(Grid_tile_t), intent(IN)         :: tileDesc

    integer :: i, j, k, var

    associate(loGC => tileDesc%grownLimits(LOW,  :), &
              hiGC => tileDesc%grownLimits(HIGH, :))
        do           var = UNK_VARS_BEGIN, UNK_VARS_END
            do         k = loGC(KAXIS), hiGC(KAXIS)
                do     j = loGC(JAXIS), hiGC(JAXIS)
                    do i = loGC(IAXIS), hiGC(IAXIS)
                        initData(var, i, j, k) = 1.1 * var
                    end do
                end do
            end do
        end do
    end associate
end subroutine Simulation_initBlock

