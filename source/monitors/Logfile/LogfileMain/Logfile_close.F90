!!****if* source/monitors/Logfile/LogfileMain/Logfile_close
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
!! NAME
!!   Logfile_close
!!
!! SYNOPSIS
!!   Logfile_close(logical(IN) :: logUnitLocal) 
!!
!! DESCRIPTION
!!   Close the log file
!!
!! ARGUMENTS
!!   logUnitLocal - indicates whether to close the local or global logfile
!!
!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
!!***

subroutine Logfile_close(logUnitLocal)

  use Logfile_data, ONLY : log_globalMe, log_fileOpen, log_lun,log_fileOpenLocal, log_lunLocal

#include "constants.h"
  
  implicit none

  logical, optional, intent(IN) :: logUnitLocal


  ! Close "local" (PE-specific) logfile instead of the global one if so requested
  if(present(logUnitLocal)) then
     if(logUnitLocal) then
        close(log_lunLocal)
        log_fileOpenLocal=.false.
        return
     end if
  end if
  
  ! We get here only if logUnitLocal was not present or if it was .FALSE.
  if(.not. log_fileOpen) return
   
  if(log_globalMe == MASTER_PE) close(log_lun)

  log_fileOpen = .false.

end subroutine Logfile_close

