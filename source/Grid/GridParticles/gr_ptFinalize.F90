!!****if* source/Grid/GridParticles/gr_ptFinalize
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
!!  gr_ptFinalize
!!
!! SYNOPSIS
!!
!!  call gr_ptFinalize()
!!
!! DESCRIPTION
!!
!!  deallocate the scratch buffers used in moving particles data
!!
!! ARGUMENTS
!!
!!  none
!!
!!***


subroutine gr_ptFinalize()
  use Logfile_interface, ONLY : Logfile_stamp
  use gr_ptData, ONLY :   gr_ptDestBuf,gr_ptSourceBuf
  use Grid_data, ONLY :   gr_useParticles
  implicit none

  integer :: istat

  if (allocated(gr_ptDestBuf) .OR. allocated(gr_ptSourceBuf) .OR. gr_useParticles) then
     ! These are allocated in gr_ptInit
     deallocate(gr_ptDestBuf,stat=istat)
     if (istat .NE. 0) then
        call Logfile_stamp('WARNING - failed to deallocate gr_ptDestBuf','[gr_ptFinalize]')
     end if
     deallocate(gr_ptSourceBuf,stat=istat)
     if (istat .NE. 0) then
        call Logfile_stamp('WARNING - failed to deallocate gr_ptSourceBuf','[gr_ptFinalize]')
     end if
  end if

  return
end subroutine gr_ptFinalize
