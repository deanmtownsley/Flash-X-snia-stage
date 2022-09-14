!!****h* source/numericalTools/MoL/localAPI/ml_memInterface
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
!!  NAME
!!
!!      ml_memInterface
!!
!!  SYNOPSIS
!!
!!      use ml_memInterface
!!
!!  DESCRIPTION
!!
!!      This is the header file for the method of lines time integration unit
!!      that defines memory storage and operations procedures
!!
!!***
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

    !! ================================ !!
    !!  Memory allocation/deallocation  !!
    !! ================================ !!

   interface
      subroutine ml_memAlloc
      end subroutine ml_memAlloc
   end interface

   interface
      subroutine ml_memFree
      end subroutine ml_memFree
   end interface

    !! ===================================== !!
    !!  Directly set or copy to/from memory  !!
    !! ===================================== !!

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

    !! ============================== !!
    !!  Linear combination operators  !!
    !! ============================== !!

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
