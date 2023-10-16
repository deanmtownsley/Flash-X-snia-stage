!!****if* source/Driver/DriverMain/MoL/dr_molPostUpdate
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
!!      dr_molPostUpdate
!!
!!  SYNOPSIS
!!
!!      call dr_molPostUpdate(real, intent(in) :: t)
!!
!!  DESCRIPTION
!!
!!      Perform any post-update (post-stage/timestep) work
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine dr_molPostUpdate(t)
   use Hydro_interface, only: Hydro_molPostUpdate
   use RadTrans_interface, only: RadTrans_molPostUpdate
   use Simulation_interface, only: Simulation_molPostUpdate

   implicit none

   real, intent(in) :: t

   call Hydro_molPostUpdate(t)
   call RadTrans_molPostUpdate(t)
   call Simulation_molPostUpdate(t)
end subroutine dr_molPostUpdate
