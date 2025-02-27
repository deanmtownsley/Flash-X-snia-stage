!!****if* source/monitors/Timers/TimersMain/MPINative/Timers_data
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
!!  Timers_data
!!
!! SYNOPSIS
!!
!!  use Timers_data
!!
!! DESCRIPTION
!!
!!  Holds the data needed by the Timers Unit
!!
!!***


module Timers_data

#include "Timers.h"

  integer, parameter :: tmr_timerInvalid = TIMERS_INVALID
  integer, parameter :: tmr_maxTimerParents = 30
  integer, parameter :: tmr_maxCallStackDepth = 30
  integer, parameter :: tmr_nameSize = 30
  integer, parameter :: tmr_maxSegments = 100

  type tmr_stack
     integer                 :: stackPointer
     integer                 :: stack(tmr_maxCallStackDepth)
  end type tmr_stack
  
  type tmr_stackList       
     type(tmr_stack)         :: stacks(tmr_maxTimerParents)
     integer                 :: lastItem
  end type tmr_stackList

  type tmr_acctSegType
     character(len=tmr_nameSize) :: name
     real                        :: time(tmr_maxTimerParents)
     real                        :: dtime(tmr_maxTimerParents)
     real                        :: pctTime(tmr_maxTimerParents)
     logical                     :: isTimed(tmr_maxTimerParents)
     integer                     :: timesCalled(tmr_maxTimerParents)
     type(tmr_stackList)         :: stacks
  end type tmr_acctSegType

  type (tmr_stack), save         :: tmr_CallStack
  type (tmr_acctSegType), save   :: tmr_acctSegs(tmr_maxSegments)

  
  integer, save                :: tmr_numSegments
  integer, save                :: tmr_globalComm, tmr_globalMe, tmr_globalNumProcs
  real, save                   :: tmr_initTime
  
  integer, parameter :: tmr_bigInt = SELECTED_INT_KIND(15)
  integer (kind=tmr_bigInt), save :: tmr_cu_block_count = 0
  integer (kind=tmr_bigInt), save :: tmr_global_block_count = 0
  
  character(len=40), save      :: tmr_initDate

  logical, save                :: tmr_writeStatSummary
  logical, save                :: tmr_eachProcWritesSummary
  logical, save                :: tmr_csvSummaryAllProcs
  logical, save                :: tmr_logSummaryAllProcs  
  logical, save                :: tmr_suppressInParallel

end module Timers_data
