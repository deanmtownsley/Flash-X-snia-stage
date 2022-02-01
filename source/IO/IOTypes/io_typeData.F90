!!****if* source/IO/IOTypes/io_typeData
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
!!  io_typeData
!!
!! SYNOPSIS
!!
!!  use io_typeData
!!
!! DESCRIPTION 
!!  
!!  Holds all the IO data that is needed by the IOtype sub unit
!!
!! ARGUMENTS
!!
!!  none    
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

!Used only with derived datatype I/O (names may soon change)
!The names will eventually change to something like io_typeXXX.
module io_typeData
  implicit none
  logical, save :: io_packMeshPlotWriteHDF5
  logical, save :: io_packMeshChkWriteHDF5
  logical, save :: io_packMeshChkReadHDF5
  logical, save :: io_asyncMeshPlotWritePnet
  logical, save :: io_asyncMeshChkWritePnet
  logical, save :: io_asyncMeshChkReadPnet
  logical, save :: io_useLegacyLabels
  integer, parameter :: io_legacyLabelLength = 4
end module io_typeData
