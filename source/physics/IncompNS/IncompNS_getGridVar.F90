!!**if***source/physics/IncompNS
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
subroutine IncompNS_getGridVar(name, value)

   use Driver_interface, ONLY: Driver_abort

   implicit none
   character(len=*), intent(in)  :: name
   integer, intent(out)          :: value

   value = -1
   print *, "Error in setting grid var: ", name
   call Driver_abort("IncompNS_getGridVar: Unknown IncompNS Grid Variable")

   return
end subroutine IncompNS_getGridVar
