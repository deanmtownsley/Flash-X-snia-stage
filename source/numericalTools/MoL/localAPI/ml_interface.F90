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
!! @brief Private interfaces for MoL
!!
!! @details This is the header file for the method of lines time integration unit
!!          that defines its private interfaces.

!> @ingroup MoLPrivate
!! Interfaces to private MoL procedures
module ml_interface

   implicit none

   interface
      subroutine ml_init
      end subroutine ml_init
   end interface

   interface
      subroutine ml_finalize
      end subroutine ml_finalize
   end interface

   interface
      subroutine ml_initTableau
      end subroutine ml_initTableau
   end interface

   interface
      subroutine ml_advance(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine ml_advance
   end interface

   interface
      subroutine ml_calcRHS(rhsType, rhsStruct, t, dtWeight)
         implicit none
         integer, intent(in) :: rhsType, rhsStruct
         real, intent(in) :: t
         real, intent(in) :: dtWeight
      end subroutine ml_calcRHS
   end interface

   interface
      subroutine ml_error(msg)
         implicit none
         character(len=*), intent(in) :: msg
      end subroutine ml_error
   end interface

   interface
      subroutine ml_warn(msg)
         implicit none
         character(len=*), intent(in) :: msg
      end subroutine ml_warn
   end interface

   interface
      subroutine ml_status(msg)
         implicit none
         character(len=*), intent(in) :: msg
      end subroutine ml_status
   end interface

end module ml_interface
