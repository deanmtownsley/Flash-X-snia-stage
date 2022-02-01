!!****if* source/IO/IOMain/io_memoryReport
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
!!  io_memoryReport
!!
!!
!! SYNOPSIS
!!
!!  call io_memoryReport()
!!
!!
!! ARGUMENTS
!!
!!
!! DESCRIPTION
!!
!!  Report on the memory usage of FLASH.  
!!
!!  The implementation simply calls ut_sysMemSummaryStats, for which an
!!  implementation is provided in flashUtitilies, which in turn calls
!!  calls getrusage or whatever else is necessary to access the information
!!  on resident memory size on the current processor and other statistics.
!!  One or several single line reports, each giving the average memory/processor,
!!  the maximum amount of memory used on a processors, and the minimum amount
!!  of memory used on a processor, are written to the logfile.
!!
!!  Normally this subroutine is called every memory_stat_freq timestep.
!!
!!
!! NOTES
!!
!!  This routine is a no-op if the preprocessor FLASH_USE_MEMORYUSAGE is undefined.
!!
!!***

#include "Simulation.h"

#ifdef FLASH_USE_MEMORYUSAGE
#include "ut_sysMem.h"
#include "constants.h"
#endif

Subroutine io_memoryReport()
#ifdef FLASH_USE_MEMORYUSAGE
  use IO_data, ONLY : io_globalComm, io_globalMe, io_measRSS
  use ut_sysMemData, ONLY : memsummary_t
  use ut_sysMemInterface, ONLY : ut_sysMemSummaryStats
  use Logfile_interface, ONLY: Logfile_stamp
  implicit none
  integer, parameter :: maxStats = 20, verbosity = 0, &
       memorySampler = UT_SYSMEM_AUTO
  character (len=*), parameter :: rssStr = 'rss'
  type(memsummary_t), dimension(maxStats) :: memSummary
  integer :: i, rssPos, numStats
  logical :: rssFound
  character (len=12) :: minMem_string, maxMem_string, avgMem_string

  call ut_sysMemSummaryStats(io_globalComm, verbosity, memorySampler, &
       memSummary, numStats)
  
  if (io_globalMe == MASTER_PE) then
     io_measRSS = -1.0
     rssFound = .false.

     do i = 1, numStats
        write (minMem_string, '(f12.2)') memSummary(i) % min
        write (maxMem_string, '(f12.2)') memSummary(i) % max
        write (avgMem_string, '(f12.2)') memSummary(i) % avg

        call Logfile_stamp(trim(memSummary(i) % description) // & 
             trim(minMem_string) // ' (min)  ' // &
             trim(maxMem_string) // ' (max)  ' // &
             trim(avgMem_string) // ' (avg) ', tag="memory")

        if (.not.rssFound) then
           rssPos = index(memSummary(i) % description, rssStr)
           if (rssPos > 0) then
              rssFound = .true.
              io_measRSS = memSummary(i) % max
           end if
        end if
     end do
  end if

#endif

end Subroutine io_memoryReport
