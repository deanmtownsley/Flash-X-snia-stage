!!****ih* source/flashUtilities/system/memoryUsage/ut_sysMemInterface
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
!!  ut_sysMemInterface 
!!
!! SYNOPSIS
!!  use ut_sysMemInterface
!!
!!***

#include "constants.h"

module ut_sysMemInterface
  implicit none
  interface
     subroutine ut_sysMemSummaryStats(comm, verbosity, memorySampler, &
          memsummary, numStats)
       use ut_sysMemData, ONLY : memsummary_t, meminfo_t
       implicit none
       integer, intent(IN) :: comm, verbosity, memorySampler
       type(memsummary_t), dimension(:), intent(OUT) :: memsummary
       integer, intent(OUT) :: numStats
     end subroutine ut_sysMemSummaryStats
  end interface

  interface
     subroutine ut_sysMemStats(verbosity, memorySamplerArg, meminfo, numStats)
       use ut_sysMemData, ONLY : meminfo_t
       implicit none
       integer, intent(IN) :: verbosity, memorySamplerArg
       type(meminfo_t), dimension(:), intent(OUT) :: meminfo
       integer, intent(OUT) :: numStats
     end subroutine ut_sysMemStats
  end interface
end module ut_sysMemInterface
