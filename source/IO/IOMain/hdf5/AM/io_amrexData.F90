!!****if* source/IO/IOMain/hdf5/AM/io_amrexData
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
!!  io_amrexData
!!
!! SYNOPSIS
!!
!!  use io_amrexData
!!
!! DESCRIPTION 
!!  
!!  Holds all the data need by source/IO/IOMain/hdf5/AM
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

module io_amrexData

#ifdef USE_IO_C_INTERFACE
  use iso_c_binding
#else
  !Older compilers do not support Fortran 2003.
  integer, parameter :: c_char = KIND('A')
#endif

  logical, save :: io_plotFileAmrexFormat

end module io_amrexData
