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

module HeatAD_interface

   implicit none

   interface
      subroutine HeatAD_init(restart)
         implicit none
         logical, intent(in) :: restart
      end subroutine HeatAD_init
   end interface

   interface
      subroutine HeatAD_finalize()
         implicit none
      end subroutine HeatAD_finalize
   end interface

   interface
      subroutine HeatAD_solve(tileDesc, dt)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
         real, INTENT(IN) :: dt
      end subroutine HeatAD_solve
   end interface

   interface
      subroutine HeatAD_advection(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine HeatAD_advection
   end interface

   interface
      subroutine HeatAD_diffusion(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine HeatAD_diffusion
   end interface

   interface
      subroutine HeatAD_indicators()
         implicit none
      end subroutine HeatAD_indicators
   end interface

   interface
      subroutine HeatAD_reInitGridVars(tileDesc)
         use Grid_tile, ONLY: Grid_tile_t
         implicit none
         type(Grid_tile_t), INTENT(IN) :: tileDesc
      end subroutine HeatAD_reInitGridVars
   end interface

   interface
      subroutine HeatAD_getGridVar(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         integer, intent(out)          :: value
      end subroutine HeatAD_getGridVar
   end interface

   interface
      subroutine HeatAD_getScalarProp(name, value)
         implicit none
         character(len=*), intent(in)  :: name
         real, intent(out)             :: value
      end subroutine HeatAD_getScalarProp
   end interface

end module HeatAD_interface
