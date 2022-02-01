!!****f* source/monitors/Profiler/Profiler_getSummary
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
!! NAME
!!   
!! Profiler_getSummary
!!
!! SYNOPSIS
!!  Profiler_getSummary(integer(in) :: nIntervals)
!!
!! DESCRIPTION
!!   Write out the Profiler summary.  This can be done in whatever
!!   format the particular profiler chooses to do it.  
!!
!! ARGUMENTS
!!
!!  nIntervals - nIntervals    :: number of subintervals timed
!!                                 typically timsteps
!! 
!!
!!
!!***


subroutine Profiler_getSummary(nIntervals)

  implicit none
  integer, intent(in) :: nIntervals

end subroutine Profiler_getSummary

