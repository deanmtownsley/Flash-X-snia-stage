!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!! @brief Public interfaces for MoL
!!
!! @details This is the header file for the method of lines time integration unit
!!          that defines its public interfaces.

!> @ingroup MoL
!! Interfaces to MoL public procedures
module MoL_interface

   implicit none

   interface
      subroutine MoL_init
      end subroutine MoL_init
   end interface

   interface
      subroutine MoL_finalize
      end subroutine MoL_finalize
   end interface

   interface
      subroutine MoL_registerVariable(name, evolIndex, rhsIndex)
         implicit none
         character(len=*), intent(in) :: name
         integer, intent(in) :: evolIndex
         integer, intent(out) :: rhsIndex
      end subroutine MoL_registerVariable
   end interface

   interface
      subroutine MoL_registerRHS(rhsType, rhsFunc)
         use MoL_functionTypes, only: MoL_rhs_t
         implicit none
         integer, intent(in) :: rhsType
         procedure(MoL_rhs_t) :: rhsFunc
      end subroutine MoL_registerRHS
   end interface

   interface
      subroutine MoL_registerUpdate(updateType, updateFunc)
         use MoL_functionTypes, only: MoL_update_t
         implicit none
         integer, intent(in) :: updateType
         procedure(MoL_update_t) :: updateFunc
      end subroutine MoL_registerUpdate
   end interface

   interface
      subroutine MoL_registerPostUpdate(postUpdateType, postUpdateFunc)
         use MoL_functionTypes, only: MoL_postUpdate_t
         implicit none
         integer, intent(in) :: postUpdateType
         procedure(MoL_postUpdate_t) :: postUpdateFunc
      end subroutine MoL_registerPostUpdate
   end interface

   interface
      subroutine MoL_releaseFunctions
      end subroutine MoL_releaseFunctions
   end interface

   interface
      subroutine MoL_advance(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine MoL_advance
   end interface

   interface
      subroutine MoL_regrid
      end subroutine MoL_regrid
   end interface

   interface
      subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
         use Grid_tile, only: Grid_tile_t
         implicit none
         class(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: dataPtr(:, :, :, :)
         integer, intent(in) :: dataStruct
      end subroutine MoL_getDataPtr
   end interface

   interface
      subroutine MoL_releaseDataPtr(tileDesc, dataPtr, dataStruct)
         use Grid_tile, only: Grid_tile_t
         implicit none
         class(Grid_tile_t), intent(in) :: tileDesc
         real, pointer :: dataPtr(:, :, :, :)
         integer, intent(in) :: dataStruct
      end subroutine MoL_releaseDataPtr
   end interface

   interface
      function MoL_getRHSIndex(evolIndex)
         implicit none
         integer :: MoL_getRHSIndex
         integer, intent(in) :: evolIndex
      end function MoL_getRHSIndex
   end interface

   interface
      function MoL_getOrder()
         implicit none
         integer :: MoL_getOrder
      end function MoL_getOrder
   end interface

end module MoL_interface
