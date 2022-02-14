!!****if* source/Driver/DriverMain/Driver_driftSetSrcLoc
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
!!  Driver_driftSetSrcLoc
!!
!! DESCRIPTION
!!
!!  Saves the current source file & line to module vars in Drift_data.
!!  It is recommended that you call this method like so:
!!    call Driver_driftSetSrcLoc(__FILE__,__LINE__)
!!
!! ARGUMENTS
!!
!!  filename: source file location to log in case of changed hash
!!  line: source line location to log in case of changed hash
!!
!!***
subroutine Driver_driftSetSrcLoc(filename, line)
  use Driver_data, only: dr_driftSrcFile, dr_driftSrcLine
  implicit none
  character(len=*), intent(in) :: filename
  integer, intent(in) :: line
  dr_driftSrcFile = filename
  dr_driftSrcLine = line
end subroutine Driver_driftSetSrcLoc
