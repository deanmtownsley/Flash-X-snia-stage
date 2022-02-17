!!****if* source/IO/localAPI/io_writeGrid
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
!!
!!
!!
!! NAME
!!
!!  io_writeGrid
!!
!!
!! SYNOPSIS
!!
!!  call io_writeGrid()
!!
!!
!! DESCRIPTION
!!
!! This function writes the grid information to an hdf5 file to store the
!! Paramesh or AMReX Grid cell coordinates (Left, Center, Right)
!! and the cell metrics for later use in post-processing
!!
!! Currently only supports hdf5 IO
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!  This is stub, does nothing
!!
!!***

subroutine io_writeGrid()
   implicit none
   return
end subroutine io_writeGrid
