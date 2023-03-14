!!****h* source/Driver/DriverMain/MoL/dr_molInterface
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!      dr_molInterface
!!
!!  SYNOPSIS
!!
!!      use dr_molInterface
!!
!!  DESCRIPTION
!!
!!      MoL-specific driver features
!!
!!***
module dr_molInterface

   implicit none

   interface
      subroutine dr_molExplicitRHS(t, activeRHS, dtWeight)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
         real, intent(in) :: dtWeight
      end subroutine dr_molExplicitRHS
   end interface

   interface
      subroutine dr_molImplicitRHS(t, activeRHS, dtWeight)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
         real, intent(in) :: dtWeight
      end subroutine dr_molImplicitRHS
   end interface

   interface
      subroutine dr_molFastRHS(t, activeRHS, dtWeight)
         implicit none
         real, intent(in) :: t
         integer, intent(in) :: activeRHS
         real, intent(in) :: dtWeight
      end subroutine dr_molFastRHS
   end interface

   interface
      subroutine dr_molImplicitUpdate(t, dt)
         implicit none
         real, intent(in) :: t, dt
      end subroutine dr_molImplicitUpdate
   end interface

   interface
      subroutine dr_molPostUpdate(t)
         implicit none
         real, intent(in) :: t
      end subroutine dr_molPostUpdate
   end interface

   interface
      subroutine dr_molPostFastUpdate(t)
         implicit none
         real, intent(in) :: t
      end subroutine dr_molPostFastUpdate
   end interface

   interface
      subroutine dr_molPreEvolve(t)
         implicit none
         real, intent(in) :: t
      end subroutine dr_molPreEvolve
   end interface

   interface
      subroutine dr_molPostTimeStep(t)
         implicit none
         real, intent(in) :: t
      end subroutine dr_molPostTimeStep
   end interface

   interface
      subroutine dr_molPostRegrid(t)
         implicit none
         real, intent(in) :: t
      end subroutine dr_molPostRegrid
   end interface

   interface
      subroutine dr_molRegisterFunctions
      end subroutine dr_molRegisterFunctions
   end interface

end module dr_molInterface
