!!****f* source/numericalTools/MoL/MoL_registerUpdate
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
!!      MoL_registerUpdate
!!
!!  SYNOPSIS
!!
!!      call MoL_registerUpdate(integer,   intent(in) :: updateType
!!                              procedure,            :: updateFunc)
!!
!!  DESCRIPTION
!!
!!      Register a function with MoL to perform an update (e.g. implicit)
!!
!! ARGUMENTS
!!
!!      updateType : One of the following function types defined in MoL.h
!!                       MOL_IMPLICIT_UPDATE   -  Implicit update (slow)
!!      updateFunc : Procedure to register
!!
!! NOTE
!!
!!      There is only one valid option for updateType (MOL_IMPLICIT_UPDATE),
!!      but this is left generic to accomodate for new update-types in the
!!      future, e.g. distinct implcit-updates for both slow- and fast-integration
!!      steps in the multi-rate integrator.
!!
!!***
subroutine MoL_registerUpdate(updateType, updateFunc)
   implicit none

   integer, intent(in) :: updateType

   interface
      subroutine updateFunc(t, dt)
         real, intent(in) :: t, dt
      end subroutine updateFunc
   end interface

   return
end subroutine MoL_registerUpdate
