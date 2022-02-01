!!****if* source/IO/IOMain/hdf5/parallel/io_closeFile
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
!!  io_closeFile
!!
!! SYNOPSIS
!!
!!  io_closeFile(integer(in) :: fileID)
!!
!! DESCRIPTION
!!  
!!  closes the hdf5 file      
!!
!! ARGUMENTS
!!
!!  fileID - integer file identifier
!!
!!***

subroutine io_closeFile( fileID)

  use io_intfTypesModule, ONLY : io_fileID_t
  implicit none

  integer(io_fileID_t), intent(in) :: fileID
  call io_h5close_file(fileID)

end subroutine io_closeFile
