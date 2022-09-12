!!****if* source/numericalTools/MoL/MoLMain/ml_functions
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
!!      ml_functions
!!
!!  SYNOPSIS
!!
!!      use ml_functions, only:
!!
!!  DESCRIPTION
!!
!!      Function pointers used internally by MoL
!!
!!***
module ml_functions

   implicit none

   abstract interface
      subroutine ml_rhs_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine ml_rhs_t

      subroutine ml_implicitUpdate_t(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine ml_implicitUpdate_t

      subroutine ml_postUpdate_t(t)
         implicit none
         real, intent(in) :: t
      end subroutine ml_postUpdate_t
   end interface

   procedure(ml_rhs_t), pointer :: ml_rhsE => ml_rhsE_default
   procedure(ml_rhs_t), pointer :: ml_rhsI => ml_rhsI_default
   procedure(ml_rhs_t), pointer :: ml_rhsF => ml_rhsF_default

   procedure(ml_implicitUpdate_t), pointer :: ml_implicitUpdate &
      => ml_implicitUpdate_default

   procedure(ml_postUpdate_t), pointer :: ml_postUpdate => ml_postUpdate_default
   procedure(ml_postUpdate_t), pointer :: ml_postUpdateFast &
      => ml_postUpdateFast_default

contains

   subroutine ml_rhsE_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_rhsE_default

   subroutine ml_rhsI_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_rhsI_default

   subroutine ml_rhsF_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_rhsF_default

   subroutine ml_implicitUpdate_default(t, dt)
      implicit none

      real, intent(in) :: t, dt

      return
   end subroutine ml_implicitUpdate_default

   subroutine ml_postUpdate_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_postUpdate_default

   subroutine ml_postUpdateFast_default(t)
      implicit none

      real, intent(in) :: t

      return
   end subroutine ml_postUpdateFast_default

end module ml_functions
