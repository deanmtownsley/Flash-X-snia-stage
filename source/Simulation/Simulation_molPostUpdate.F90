!!****f* source/Simulation/Simulation_molPostUpdate
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
!!      Simulation_molPostUpdate
!!
!!  SYNOPSIS
!!
!!      call Simulation_molPostUpdate(real, intent(in) :: t)
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
subroutine Simulation_molPostUpdate(t)
   implicit none

   real, intent(in) :: t

   return
end subroutine Simulation_molPostUpdate
