!!****ih* source/monitors/Timers/TimersMain/Tau/tmr_interface
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
!!***
Module tmr_interface
#include "Simulation.h"
#include "constants.h"

  interface tmr_findTimerIndex
     subroutine tmr_findTimerIndex(name, createIfNone, index) 
       character(len=*), intent(in) :: name
       logical, intent(in)          :: createIfNone
       integer, intent(out) :: index
     end subroutine tmr_findTimerIndex
  end interface

  interface tmr_etime
     subroutine tmr_etime(time)
       real, intent(OUT) :: time
     end subroutine tmr_etime
  end interface

end Module tmr_interface
