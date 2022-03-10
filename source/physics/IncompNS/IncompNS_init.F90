!!****f* source/physics/IncompNS/IncompNS_init
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
!!  IncompNS_init
!!
!!
!! SYNOPSIS
!!
!!  call IncompNS_init(restart)
!!
!!
!! DESCRIPTION
!!
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initAll.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!!
!!***

subroutine IncompNS_init(restart)

   implicit none
   logical, intent(IN) :: restart

end subroutine IncompNS_init

