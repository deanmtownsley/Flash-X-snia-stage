!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_localPostSendMsg
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
!!  pl_localPostSendMsg
!!
!! SYNOPSIS
!! 
!!  call pl_localPostSendMsg (integer, intent (in) :: channel)
!!
!! DESCRIPTION
!!
!!  Posts a sending message for the specified channel on the local processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to post the send
!!
!!***

subroutine pl_localPostSendMsg (channel)

  use Pipeline_data,     ONLY : pl_channelSize,        &
                                pl_comm,               &
                                pl_doLog,              &
                                pl_itemSize,           &
                                pl_logUnit,            &
                                pl_procStatusLocal,    &
                                pl_procList,           &
                                pl_rank,               &
                                pl_sendBuf,            &
                                pl_sendCount,          &
                                pl_sendRequest,        &
                                pl_tag
                                
  use Driver_interface,  ONLY : Driver_checkMPIErrorCode

  implicit none

#include "Pipeline.h"
 include "Flashx_mpi.h"

  integer, intent (in) :: channel

  integer :: error
  integer :: msgSize
  integer :: procID
!
!
!     ...Post the send and write this action to the log file (if requested).
!        Record this event into the global status array.
!
!
  msgSize = pl_sendCount (channel)

  if (msgSize >= 0) then

      procID = pl_procList (channel)

      call MPI_Isend (pl_sendBuf (1,1,channel),     &
                      pl_itemSize * msgSize,        &
                      FLASH_REAL,                   &
                      procID,                       &
                      pl_tag,                       &
                      pl_comm,                      &
                      pl_sendRequest (channel),     &
                      error                         )

      call Driver_checkMPIErrorCode (error)

      if (pl_doLog) then
          write (pl_logUnit,'(a,i6)') ' Posted sending message   to proc ID ', procID
      end if

      pl_procStatusLocal (pl_rank, PL_STATUS_COMM) = pl_procStatusLocal (pl_rank, PL_STATUS_COMM) + 1

  end if
!
!
!    ...Ready!
!
!
  return
end subroutine pl_localPostSendMsg
