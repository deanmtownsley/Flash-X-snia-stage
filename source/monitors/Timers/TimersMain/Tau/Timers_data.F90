!!****if* source/monitors/Timers/TimersMain/Tau/Timers_data
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

#include "constants.h"

  integer, parameter :: tmr_MAX_CUSTOM_TIMERS = 500
  character(len=*), parameter :: tmr_customPrefix = "*** custom:"

  type tauTimerObj
     integer, dimension(2) :: tauSavedData
     character (len=MAX_STRING_LENGTH) :: tauString
     logical :: timerStarted
  end type tauTimerObj

  type (tauTimerObj), save, &
       dimension(tmr_MAX_CUSTOM_TIMERS) :: tmr_tauList

  integer, save :: tmr_freeSlot, tmr_globalMe, tmr_globalNumProcs, tmr_prefixLen

end module Timers_data
