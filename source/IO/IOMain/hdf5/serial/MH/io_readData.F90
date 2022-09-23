!!****if* source/IO/IOMain/hdf5/serial/MH/io_readData
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
!!  io_readData
!!
!!
!! SYNOPSIS
!!
!!  call io_readData()
!!
!!
!! DESCRIPTION
!!
!!  This is the reading counterpart to io_writeData.  It reads an HDF5
!!  file and distributes it to the processors to restart a simulation.
!!
!!  this is the reading counterpart to checkpoint_wr.  It eats the HDF
!!  output and distributes it to the processors.
!!
!!  Subroutine to read checkpoint file using AMR package.
!!  Currently reads are done serially by processor 0 and data sent to
!!  other processors.
!!
!! ARGUMENTS
!!
!!   none
!!
!!***

subroutine io_readData()
    use Driver_interface, ONLY : Driver_abort

    CALL Driver_abort("[io_readData] Not implemented yet")
end subroutine io_readData

