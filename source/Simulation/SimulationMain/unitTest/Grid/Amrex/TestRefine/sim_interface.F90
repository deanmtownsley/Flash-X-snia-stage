module sim_interface
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

    interface
        subroutine sim_writeDataPoints(initData, tileDesc, points, values)
            use Grid_tile, ONLY : Grid_tile_t
            implicit none
            real,              intent(IN), pointer :: initData(:, :, :, :)
            type(Grid_tile_t), intent(IN)          :: tileDesc
            real,              intent(IN)          :: points(:, :)
            real,              intent(IN)          :: values(:)
        end subroutine sim_writeDataPoints
    end interface

    interface
        subroutine sim_collectLeaves
            implicit none
        end subroutine sim_collectLeaves
    end interface 

    interface
        subroutine sim_printLeaves(title)
            implicit none
            character(*), intent(IN)    :: title
        end subroutine sim_printLeaves
    end interface 

    interface
        subroutine sim_advance(step, points, values, set_msg, leaf_msg)
            implicit none
            integer,      intent(IN)    :: step
            real,         intent(IN)    :: points(:, :)
            real,         intent(IN)    :: values(:)
            character(*), intent(IN)    :: set_msg
            character(*), intent(IN)    :: leaf_msg
        end subroutine sim_advance
    end interface 

end module sim_interface

