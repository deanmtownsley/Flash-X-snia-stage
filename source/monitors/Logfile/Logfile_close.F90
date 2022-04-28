!!****f* source/monitors/Logfile/Logfile_close
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

  implicit none

  logical, optional, intent(IN) :: logUnitLocal

end subroutine Logfile_close

