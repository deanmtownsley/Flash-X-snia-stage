!!****if* source/IO/IOMain/IO_finalize
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
!!  IO_finalize
!!
!! SYNOPSIS
!!
!!  IO_finalize()
!!
!! DESCRIPTION
!!
!!  This function cleans up the IO unit, deallocates memory, etc.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine IO_finalize()
  implicit none

! For async case it signals to wait until all async io is complete
#ifdef FLASH_IO_ASYNC_HDF5
  call io_h5_close_async()
#endif

#ifdef FLASH_IO_EXPERIMENTAL
  call io_typeFinalize()
#endif
end subroutine IO_finalize
