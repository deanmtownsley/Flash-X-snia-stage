!!****if* source/monitors/Debugger/DebuggerMain/Debugger_init
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
!!  Debugger_init
!!
!! SYNOPSIS
!!
!!  Debugger_init()
!!                   
!!  
!! DESCRIPTION 
!!  
!!  Initialize the debugger unit
!!  
!! ARGUMENTS 
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine Debugger_init()
  use Debugger_data, ONLY : dbg_globalMe, dbg_doHeapCheck
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype
  implicit none

  call Driver_getMype(GLOBAL_COMM,dbg_globalMe)

  call RuntimeParameters_get("doHeapCheck", dbg_doHeapCheck)
#if defined(DMALLOC) || defined(DDT)
  if (dbg_doHeapCheck) then
     call dbg_enable_heap_check()
  else
     call dbg_disable_heap_check()
  end if
#else
  if (dbg_globalMe == MASTER_PE) then
     print *, "*** DMALLOC / DDT support not compiled into this application ***"
  end if
#endif
end subroutine Debugger_init
