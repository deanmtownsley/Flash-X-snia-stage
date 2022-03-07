!!****if* source/IO/IOMain/io_prepareSimInfo
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
!!  io_prepareSimInfo
!!
!!
!! SYNOPSIS
!!
!!  io_prepareSimInfo()
!!          
!!          
!!
!!
!!
!! DESCRIPTION
!!
!!  This function prepares the simulation info for io output.
!!  This includes things like, the setup call, file creation time,
!!  build directory and other statistics about the run
!!  
!!  
!!
!! ARGUMENTS
!!   none   
!!
!!***


subroutine io_prepareSimInfo() 
  
  use IO_data, ONLY : io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, io_fflags, &
       io_setupTimeStamp, io_buildTimeStamp
  

#include "Flashx_mpi_implicitNone.fh"
#include "constants.h"

  character(len=80)   :: setup_flashRelease


  io_flashRelease = setup_flashRelease()


  call setup_buildstats(io_buildDate, io_buildDir, io_buildMachine, io_setupCall, &
             io_cflags, io_fflags)

  call setup_buildstamp (io_setupTimeStamp, io_buildTimeStamp, MAX_STRING_LENGTH)
  
  call current_date_time(io_fileCreationTime)

  call removeNullChar(io_buildDate)
  call removeNullChar(io_buildDir)
  call removeNullChar(io_buildMachine)
  call removeNullChar(io_setupCall)
  call removeNullChar(io_cflags)
  call removeNullChar(io_fflags)
  call removeNullChar(io_setupTimeStamp)
  call removeNullChar(io_buildTimeStamp)
  call removeNullChar(io_fileCreationTime)
        


end subroutine io_prepareSimInfo

