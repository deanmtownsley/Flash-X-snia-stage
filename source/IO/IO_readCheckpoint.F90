!!****f* source/IO/IO_readCheckpoint
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
!!  IO_readCheckpoint
!!
!!
!! SYNOPSIS
!!
!!  IO_readCheckpoint()
!!
!!
!! DESCRIPTION
!!
!!  IO_readCheckpoint is a generic subroutine that retrieves
!!  the unklabels and then calls io_readData
!!  which is specific to pnetcdf, hdf5 and
!!  the necessary grid, UG, paramesh etc.
!!  io_readData reads a checkpoint file and reinitializes
!!  grid and scalar values to resume the run
!!  
!!
!!
!! ARGUMENTS
!!  
!!
!!
!!***

subroutine IO_readCheckpoint()


implicit none
 
  return
end subroutine IO_readCheckpoint

