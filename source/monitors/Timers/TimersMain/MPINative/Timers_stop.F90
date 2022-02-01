!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_stop
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
!!  Timers_stopIndex - stop a timer given an integer key
!!
!! SYNOPSIS
!!
!!  Timers_stopIndex(integer(IN) :: i)
!!
!! DESCRIPTION
!!  Stop timing a timer specified by a supplied integer key.  This
!!  implementation also calls Profiler_stop with the corresponding
!!  timer name.
!!
!! ARGUMENTS
!!  i --   an integer key specifiying the timer to stop
!!
!! PARAMETERS
!!
!! NOTES
!!
!!  Multithreaded timers
!!
!!  * The following concerns can be ignored if the runtime parameter
!!    tmr_suppressInParallel is set TRUE.  In that case, all
!!    Timers_start() and Timers_stop() calls within active parallel
!!    regions will simply be ignored.  Otherwise:
!!
!!  * Do not call Timers_start and Timers_stop within an OpenMP single
!!    section because timing information is only recorded for the
!!    master thread.  It is possible that the single thread is the
!!    master thread on some MPI ranks but a non-master thread on other
!!    MPI ranks.  When this happens the timer stacks will be different
!!    across MPI ranks and we will get the 'whaddayamakeofthat?'
!!    message.
!!
!!  * Do not replace omp master with an omp single in the Timers unit
!!    because it may cause a deadlock.
!!
!!  * Do not replace omp master with a nowait omp single in the Timers
!!    unit because the Timers unit is of unknown thread-safety.  This
!!    would be a problem if different threads are in different nowait
!!    single sections.
!!
!!***

recursive subroutine Timers_stopIndex (i)

  use Timers_data, ONLY:  tmr_callStack, tmr_acctSegs
  use Timers_data, ONLY: tmr_suppressInParallel
  use Timers_interface, ONLY : Timers_startIndex
!$ use omp_lib

  implicit none  

  integer, intent(in) :: i
  integer            :: result
  integer            :: j
  real :: temp_time

#ifdef NOOP
  return
#endif

  !$ if (tmr_suppressInParallel .AND. omp_get_num_threads() > 1) RETURN

  !$omp master
  call tmr_stackTop(tmr_callStack, result)
  if (i == result) then
     call tmr_stackPop(tmr_callStack, j)
     if (j >= 0) then
!!$       call profile_end(tmr_acctSegs(i)%name, i)      
        call tmr_stackListIndex(tmr_acctSegs(i)%stacks, tmr_callStack, j)
        call tmr_etime(temp_time)
        tmr_acctSegs(i)%time(j) = temp_time - tmr_acctSegs(i)%dtime(j) + tmr_acctSegs(i)%time(j)
        tmr_acctSegs(i)%isTimed(j) = .false.
     end if
  else
     call tmr_stackTop(tmr_callStack, j)
     if (j >= 0) then
        call Timers_stopIndex(j)
        call Timers_stopIndex(i)
        call Timers_startIndex(j)
     end if
  end if
  !$omp end master
  return

end subroutine Timers_stopIndex

!!****if* source/Timers/TimersMain/Timers_stopString
!!
!! NAME
!!   Timers_stopString
!!
!! SYNOPSIS
!!
!!   Timers_stopString(character(IN) :: name(:))
!!
!! DESCRIPTION
!!
!!   Stop timing a timer specified by a supplied name. This
!!   implementation also calls Profiler_stop.
!!
!! ARGUMENTS
!!   name --   a name specifiying the timer to stop
!!
!! PARAMETERS
!!
!!***
subroutine Timers_stopString(name)
  use Timers_interface, ONLY : Timers_stopIndex

  implicit none  

  character(len=*), intent(IN)   :: name
  integer            :: i
  
  !$omp master
  call tmr_findTimerIndex (name,.FALSE., i)
  !$omp end master
  call Timers_stopIndex(i)
  
  return
end subroutine Timers_stopString
