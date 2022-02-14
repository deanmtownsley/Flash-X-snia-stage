!!****f* source/monitors/Timers/Timers_getSummary
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
!! Timers_getSummary
!!
!! SYNOPSIS
!!
!!  call Timers_getSummary(integer(in) :: nIntervals)
!!
!! DESCRIPTION
!!  Write out the Timers summary to the logfile.
!!  Timers collects the performance data and then calls
!!  Logfile_writeSummary to do the actual formatting.
!!
!!  It should be safe to call this at any point in the 
!!  simulation-- an updated summary will appear in the 
!!  logfile.  
!!
!! ARGUMENTS
!!
!!  nIntervals - number of subintervals timed, which is determined 
!!               by the caller of this code, but this will 
!!               typically be the number of timesteps taken
!!               since the last time Timers_init was called. 
!!
!!***


subroutine Timers_getSummary( nIntervals)

  implicit none

  integer, intent(in) :: nIntervals

  return

end subroutine Timers_getSummary



