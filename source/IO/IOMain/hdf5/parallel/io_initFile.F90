!!****if* source/IO/IOMain/hdf5/parallel/io_initFile
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
!!  io_initFile
!!
!! SYNOPSIS
!!
!!  call io_initFile(integer(in)                 :: filenum,
!!                   integer(io_fileID_t)(INOUT) :: fileID,
!!                   character(INOUT)            :: filename(MAX_STRING_LENGTH),
!!                   character(len=*)(in)        :: outputType,
!!                   logical(in)                 :: forced)
!!
!! DESCRIPTION
!!
!!  
!!  Initialized the hdf5 file
!!
!! ARGUMENTS
!!
!!
!!  filenum - number of the output file
!!
!!  fileID - file handle returned from hdf5 init file calls
!!
!!  filename - name of the returned file made up of the outputType, filenum and basename
!!
!!  outputType - string indicating output type, usually "chk" or "plt_cnt" or "part"
!!
!!  forced - .true. if this file is considered "forced."
!!
!!***


subroutine io_initFile( filenum, fileID, filename, outputType, forced)

  use IO_data, ONLY : io_comm, io_outputSplitNum
  use Driver_interface, ONLY : Driver_abort
  use io_intfTypesModule, ONLY : io_fileID_t

  implicit none
  
#include "constants.h"

  integer, intent(IN) :: filenum
  integer(io_fileID_t), intent(INOUT) :: fileID
  character(len=*), intent(in) :: outputType
  character (len=MAX_STRING_LENGTH), intent(inout) :: filename
  logical, intent(IN) :: forced
  integer :: existing

  call io_getOutputName(filenum, "hdf5", outputType, filename, forced)

  existing = 0
  fileID = -1
  call io_h5init_file(fileID, filename, io_comm, io_outputSplitNum, existing)
  if(fileID == -1) then
     print *, "Error: unable to initialize file"
     call Driver_abort("unable to initialize hdf5 file")
  end if

end subroutine io_initFile
