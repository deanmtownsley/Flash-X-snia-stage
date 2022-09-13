!!****if* source/numericalTools/MoL/MoL_functionTypes
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
!!      MoL_functionTypes
!!
!!  SYNOPSIS
!!
!!      use MoL_functionTypes, only:
!!
!!  DESCRIPTION
!!
!!      Commonly used interfaces for procedures registered and utilized by MoL
!!
!!***
module MoL_functionTypes

   implicit none

   abstract interface
      subroutine MoL_rhs_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine MoL_rhs_t
   end interface

   abstract interface
      subroutine MoL_update_t(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine MoL_update_t
   end interface

   abstract interface
      subroutine MoL_postUpdate_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine MoL_postUpdate_t
   end interface
end module MoL_functionTypes
