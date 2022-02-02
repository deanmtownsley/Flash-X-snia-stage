!!****f* source/flashUtilities/system/memoryUsage/UnitTest1/ut_sysMemUnitTest1
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
!!  ut_sysMemUnitTest1
!!
!! SYNOPSIS
!!  ./ut_sysMemUnitTest1 (standalone program)
!!
!! DESCRIPTION
!!  This is a standalone unit test independent of FLASH.
!!
!!***

#include "ut_sysMem.h"

Program ut_sysMemUnitTest1
  implicit none
  integer, parameter :: XSIZE=10000
  integer, parameter :: YSIZE=10000
  real, dimension(:,:), allocatable  :: matrix

  call print_memory_stats("START")

  allocate(matrix(XSIZE,YSIZE))
  matrix = -1.0
  call print_memory_stats("Just allocated a big heap array")

  deallocate(matrix)
  call print_memory_stats("Just deallocated a big heap array")
End Program ut_sysMemUnitTest1


subroutine print_memory_stats(msg)
  use ut_sysMemData, ONLY : meminfo_t
  use ut_sysMemInterface, ONLY  : ut_sysMemStats
  implicit none
  character(len=*), intent(IN) :: msg

  integer, parameter :: verbosity = 1, meminfoSize = 20
  type(meminfo_t), dimension(meminfoSize) :: meminfo
  integer :: numStats, i
  integer, parameter :: memorySampler = UT_SYSMEM_ALL
  
  call ut_sysMemStats(verbosity, memorySampler, meminfo, numStats)

  write(*,*) msg
  do i = 1, numStats
     write (*,'(a,f12.2)') trim(meminfo(i) % description), &
          meminfo(i) % measurement
  end do
  write(*,*) ""
end subroutine print_memory_stats
