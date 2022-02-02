!!****f* source/monitors/Logfile/Logfile_create
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
!!     Logfile_create
!! 
!! SYNOPSIS
!!     Logfile_create()
!!
!! DESCRIPTION
!!
!!     Creates the named log file and writes some header information
!!     to it, including the included units, runtime parameters,
!!     physical constants.  Meta data about the run is also stored
!!     like a time stamp, the run dimensionality, compiler flags and
!!     more.  The logfile can be stamped from any unit to store the
!!     simulation's progress.
!!
!!     The name of the parameter file is taken as an input; it is
!!     echoed to the log file.  Only the master processor actually
!!     writes anything.  In order to avoid accidentally overwriting an
!!     important logfile during a science run, the logfile is always
!!     opened in append mode.
!!
!!     
!!
!! ARGUMENTS
!!     
!!
!!***

!! NOTES
!!  variables that begin with "log_" are defined in the fortran 
!!  module Logfile_data.  The prefix "log_" is meant to indicate
!!  that these variables have Logfile unit scope.  Other variables
!!  are local to the individual subroutines
!!
subroutine Logfile_create ()

  implicit none


  return
end subroutine Logfile_create
