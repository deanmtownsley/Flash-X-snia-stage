!!****if* source/physics/ImBound/ImBoundMain/ImBound_getBodyPtr
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
!! NAME
!!
!!  ImBound_data
!!
!!
!! SYNOPSIS
!!
!!  MODULE ImBound_data()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  This stores data and limiter functions that are specific to the ImBound module.
!!
!!***
subroutine ImBound_getBodyPtr(bodyPtr, ibd)

   use ImBound_type, only: ImBound_type_t
   use ImBound_data, only: ib_bodyInfo
   use Driver_interface, only: Driver_abort

   type(ImBound_type_t), pointer :: bodyPtr
   integer, intent(in) :: ibd 

   ! Avoid possible memory leaks
   if (associated(bodyPtr)) then
      call Driver_abort("[ImBound_getBodyPtr] Given data pointer must be NULL")
   end if

   bodyPtr => ib_bodyInfo(ibd)

end subroutine ImBound_getBodyPtr
