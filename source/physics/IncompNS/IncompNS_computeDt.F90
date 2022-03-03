!!****f* source/physics/IncompNS/IncompNS_computeDt
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
!!  IncompNS_computeDt
!!
!!
!! SYNOPSIS
!!
!!
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!***

subroutine IncompNS_computeDt(ins_mindt, ins_minloc)
   implicit none
   real, intent(INOUT) :: ins_mindt
   integer, intent(INOUT) :: ins_minloc(5)
end subroutine IncompNS_computeDt
