!!****if* source/Driver/DriverMain/MoL/dr_molPostFastUpdate
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
!!      dr_molPostFastUpdate
!!
!!  SYNOPSIS
!!
!!      call dr_molPostFastUpdate(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Perform any post-fast-update (post-stage/timestep) work
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine dr_molPostFastUpdate(t)
   use Hydro_interface, only: Hydro_molPostFastUpdate
   use RadTrans_interface, only: RadTrans_molPostFastUpdate
   use Simulation_interface, only: Simulation_molPostFastUpdate

   implicit none

   real, intent(in) :: t

   call Hydro_molPostFastUpdate(t)
   call RadTrans_molPostFastUpdate(t)
   call Simulation_molPostFastUpdate(t)
end subroutine dr_molPostFastUpdate
