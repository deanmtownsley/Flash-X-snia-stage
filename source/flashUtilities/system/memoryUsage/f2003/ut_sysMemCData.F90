!!****if* source/flashUtilities/system/memoryUsage/f2003/ut_sysMemCData
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
!!  ut_sysMemCData
!!
!! SYNOPSIS
!!  use ut_sysMemCData
!!   
!!***

#include "constants.h"

module ut_sysMemCData
  use iso_c_binding, ONLY : c_double, c_ptr
  implicit none

  type, bind(c) :: c_meminfo_t
     real(c_double) :: measurement
     type(c_ptr) :: description
  end type c_meminfo_t
end module ut_sysMemCData
