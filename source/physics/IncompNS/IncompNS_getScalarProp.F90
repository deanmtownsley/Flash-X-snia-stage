!!****
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
!!
!!***

#include "constants.h"

subroutine IncompNS_getScalarPropReal(name, value)

   use Driver_interface, ONLY: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   real, intent(out)             :: value

   character(len=100)            :: errorMessage

   value = 0.
   write (errorMessage, *) '[IncompNS_getScalarProp] Unknown scalar: ', name
   call Driver_abort(errorMessage)

end subroutine IncompNS_getScalarPropReal

subroutine IncompNS_getScalarPropInteger(name, value)

   use Driver_interface, ONLY: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   integer, intent(out)          :: value

   character(len=100)            :: errorMessage

   value = 0
   write (errorMessage, *) '[IncompNS_getScalarProp] Unknown scalar: ', name
   call Driver_abort(errorMessage)

end subroutine IncompNS_getScalarPropInteger

subroutine IncompNS_getScalarPropLogical(name, value)

   use Driver_interface, ONLY: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   logical, intent(out)          :: value

   character(len=100)            :: errorMessage

   value = .false.
   write (errorMessage, *) '[IncompNS_getScalarProp] Unknown scalar: ', name
   call Driver_abort(errorMessage)

end subroutine IncompNS_getScalarPropLogical
