!!****if* source/IO/IOTypes/io_typeFinalize
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
!!  io_typeFinalize
!!
!! SYNOPSIS
!!
!!  io_typeFinalize()
!!
!! DESCRIPTION
!!
!!  Cleans up IOTypes sub unit.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine io_typeFinalize()
#ifdef USE_IO_C_INTERFACE
  use io_c_type_interface, ONLY : io_free_grid_mpi_types
#ifdef FLASH_IO_PNETCDF
  use io_c_type_interface, ONLY : io_ncmpi_nonblocking_finalize
#endif
#endif
  implicit none

  call io_free_grid_mpi_types()

#ifdef FLASH_IO_PNETCDF
  call io_ncmpi_nonblocking_finalize()
#endif
end subroutine io_typeFinalize
