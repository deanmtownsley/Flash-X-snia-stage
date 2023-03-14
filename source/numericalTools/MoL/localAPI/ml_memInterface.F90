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
!! @brief Private interfaces for MoL's scratcy memory system
!!
!! @details This is the header file for the method of lines time integration unit
!!          that defines its private interfaces related to scratch memory.

!> @ingroup MoLPrivate
!! Interfaces to private procedures implementing MoL's memory management
module ml_memInterface

   implicit none

   interface
      subroutine ml_memSetActiveRHS(irhs)
         implicit none
         integer, intent(in) :: irhs
      end subroutine ml_memSetActiveRHS
   end interface

   interface
      subroutine ml_memReleaseActiveRHS()
      end subroutine ml_memReleaseActiveRHS
   end interface

   interface
      subroutine ml_memAlloc
      end subroutine ml_memAlloc
   end interface

   interface
      subroutine ml_memFree
      end subroutine ml_memFree
   end interface

   interface
      subroutine ml_memZero(dst)
         implicit none
         integer, intent(in) :: dst
      end subroutine ml_memZero
   end interface

   interface
      subroutine ml_memCopy(dst, src)
         implicit none
         integer, intent(in) :: dst, src
      end subroutine ml_memCopy
   end interface

   interface
      subroutine ml_memAddToVars(dst, dstFac, nsrcs, srcs, facs)
         implicit none
         integer, intent(in) :: dst
         real, intent(in) :: dstFac
         integer, intent(in) :: nsrcs
         integer, intent(in) :: srcs(nsrcs)
         real, intent(in) :: facs(nsrcs)
      end subroutine ml_memAddToVars
   end interface
end module ml_memInterface
