!!****f* source/Driver/Driver_driftUnk
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
!!  Driver_driftUnk
!!
!! DESCRIPTION
!!
!!  Compute one hash per unk variable over all registered block hashes.  Vars with
!!  hashes that have changed since the previous call will be logged.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  flags: bitmask of options (for now there is only one flag)
!!    - DRIFT_NO_PARENTS: if present then only leaf blocks are included
!!
!!***
subroutine Driver_driftUnk(src_file, src_line, flags)
  implicit none
  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in), optional :: flags
end subroutine Driver_driftUnk
