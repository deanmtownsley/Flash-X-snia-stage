!!****if* source/IO/IOMain/hdf5/IO_startRayWrite
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
!!    IO_startRayWrite
!!
!! SYNOPSIS
!!
!!    IO_startRayWrite()
!!
!! DESCRIPTION
!!
!!   This routine reopens the plot file so that laser rays can be
!!   written to it. It also creates the extendible RayData dataset in
!!   the HDF5 file by calling io_h5create_raydset.
!!
!!***
subroutine IO_startRayWrite()
  use IO_data, ONLY: io_wrotePlot,       &
                     io_oldPlotFileName, &
                     io_meshComm,        &
                     io_outputSplitNum,  &
                     io_rayFileID

  use Driver_interface, ONLY: Driver_abort
  implicit none

#include "constants.h"

  integer :: existing

  if(.not. io_wrotePlot) then
     call Driver_abort("[IO_startRayWrite] Rays can only be written after a plot")
  end if

  ! Re-open the HDF5 plot file:
  existing = 1
  io_rayFileID = -1
  call io_h5init_file(io_rayFileID, io_oldPlotFileName, io_meshComm, io_outputSplitNum, existing)
  if(io_rayFileID == -1) then
     call Driver_abort("[IO_writeRays] unable to open hdf5 file: " // &
          trim(io_oldPlotFileName))
  end if

  ! Create an extendible dataset to store ray data:
  call io_h5create_raydset(io_rayFileID)          

end subroutine IO_startRayWrite
