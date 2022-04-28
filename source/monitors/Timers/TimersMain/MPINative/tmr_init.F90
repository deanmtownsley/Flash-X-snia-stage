!!****if* source/monitors/Timers/TimersMain/MPINative/tmr_init
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
!!   tmr_init
!!
!! SYNOPSIS
!!   tmr_init()
!!
!! DESCRIPTION
!!   Initialize the timer data structures.  This will
!!   essentially delete all information previously gathered by all timers
!!   and make it safe to start timers from scratch.  In the middle of 
!!   a run, for instance, this could be called once per timestep along with
!!   Timers_getSummary to get timer summary information for each timestep. 
!!
!! ARGUMENTS
!!
!! PARAMETERS
!!
!!***

subroutine tmr_init()

  use Timers_data, ONLY: tmr_acctSegs, tmr_callStack, tmr_maxSegments, &
       tmr_maxTimerParents, tmr_numSegments, tmr_initDate, tmr_initTime
  implicit none
  integer   :: i, j


#ifdef NOOP
  return
#endif
  tmr_numSegments = 0
  call tmr_stackZero(tmr_callStack)
  do i = 1, tmr_maxSegments
     tmr_acctSegs(i)%name     = "***************"
     call tmr_stackListZero(tmr_acctSegs(i)%stacks)
     do j = 1, tmr_maxTimerParents
        tmr_acctSegs(i)%time(j)     = 0.
        tmr_acctSegs(i)%dtime(j)    = 0.
        tmr_acctSegs(i)%pctTime(j)    = 0.
        tmr_acctSegs(i)%timesCalled(j) = 0
        tmr_acctSegs(i)%isTimed(j) = .false.
     enddo
  enddo
  tmr_acctSegs(tmr_maxSegments)%name = "Excess"
  
  call current_date_time(tmr_initDate)
  call tmr_etime(tmr_initTime)
!!$  ! initialize the profiling library
!!$    call profile_initialize()

  return

end subroutine tmr_init
