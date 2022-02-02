!!****f* source/monitors/Timers/Timers_reset
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
!!  Timers_reset
!!
!! SYNOPSIS
!!
!!  Timers_reset()
!!  
!! DESCRIPTION 
!!  
!!   Reset the accumulated data in timers.  This does not delete
!!   timers, so integer mappings to named timers will remain intact,
!!   but all information previously gathered about timers will be
!!   lost.  In the middle of a run, for instance, this could be called
!!   once per timestep along with Timers_getSummary to get timer
!!   summary information for each timestep.
!!  
!! ARGUMENTS 
!!
!!***

subroutine Timers_reset()
  implicit none

  return
end subroutine Timers_reset

