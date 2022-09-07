!!****if* source/numericalTools/MoL/MoLMain/MoL_functions
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
!!      MoL_functions
!!
!!  SYNOPSIS
!!
!!      use MoL_functions
!!
!!  DESCRIPTION
!!
!!      Function pointers used internally by MoL
!!
!!***

module MoL_functions

   use Grid_tile, only: Grid_tile_t

   implicit none

   abstract interface
      subroutine MoL_rhs_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine MoL_rhs_t

      subroutine MoL_implicitUpdate_t(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine MoL_implicitUpdate_t

      subroutine MoL_postUpdate_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine MoL_postUpdate_t
   end interface

   procedure(MoL_rhs_t), pointer :: MoL_rhsE => MoL_rhsE_default
   procedure(MoL_rhs_t), pointer :: MoL_rhsI => MoL_rhsI_default
   procedure(MoL_rhs_t), pointer :: MoL_rhsF => MoL_rhsF_default

   procedure(MoL_implicitUpdate_t), pointer :: MoL_implicitUpdate &
      => MoL_implicitUpdate_default

   procedure(MoL_postUpdate_t), pointer :: MoL_postUpdate => MoL_PostUpdate_default
   procedure(MoL_postUpdate_t), pointer :: MoL_postUpdateFast &
      => MoL_PostUpdateFast_default

contains

   subroutine MoL_rhsE_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine MoL_rhsE_default

   subroutine MoL_rhsI_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine MoL_rhsI_default

   subroutine MoL_rhsF_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine MoL_rhsF_default

   subroutine MoL_implicitUpdate_default(t, dt)
      implicit none

      real, intent(in) :: t, dt

      return
   end subroutine MoL_implicitUpdate_default

   subroutine MoL_PostUpdate_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine MoL_PostUpdate_default

   subroutine MoL_PostUpdateFast_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine MoL_PostUpdateFast_default

end module MoL_functions
