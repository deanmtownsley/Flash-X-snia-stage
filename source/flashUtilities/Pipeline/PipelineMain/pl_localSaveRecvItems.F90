!!****if* source/flashUtilities/Pipeline/PipelineMain/pl_localSaveRecvItems
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
!!  pl_localSaveRecvItems
!!
!! SYNOPSIS
!! 
!!  call pl_localSaveRecvItems (integer, intent (in)  :: channel,
!!                             logical, intent (out) :: isSaved)
!!
!! DESCRIPTION
!!
!!  Saves all received items in the specified channel on the local processor.
!!
!! ARGUMENTS
!!
!!  channel : channel index for which to save the received items
!!  isSaved : if true, all items have been successfully saved
!!
!!***

subroutine pl_localSaveRecvItems (channel, isSaved)

  use Pipeline_data,     ONLY : pl_doLog,            &
                                pl_itemBuf,          &
                                pl_itemCount,        &
                                pl_itemSize,         &
                                pl_logUnit,          &
                                pl_maxItems,         &
                                pl_procStatusLocal,  &
                                pl_procList,         &
                                pl_rank,             &
                                pl_recvBuf,          &
                                pl_recvCount

#include "Pipeline.h"

  implicit none

  integer, intent (in)  :: channel
  logical, intent (out) :: isSaved

  integer :: bufferBeg, bufferEnd
  integer :: n, nItems
  integer :: procID
!
!
!     ...Save the received items (if any) and write this action to the log file (if requested).
!
!
  nItems = pl_recvCount (channel)

  if (nItems > 0) then

      bufferBeg = pl_itemCount + 1
      bufferEnd = pl_itemCount + nItems

      if (bufferEnd <= pl_maxItems) then

          n = pl_itemSize

          pl_itemBuf (1:n , bufferBeg : bufferEnd) = pl_recvBuf (1:n , 1:nItems , channel)

          pl_itemCount = bufferEnd
          pl_recvCount (channel) = 0

          isSaved = .true.

          if (pl_doLog) then

              procID = pl_procList (channel)

              write (pl_logUnit,'(a,i6,2(a,2(i6)),a)') ' Saved receive items     from procID ',procID, & 
                                                       ' : copy from msg slice (',1,nItems,            &
                                                       ' ) to buffer slice (', bufferBeg, bufferEnd, ')'
          end if

      else
          isSaved = .false.
      end if

  end if
!
!
!     ...Update item buffer status on local processor.
!
!
  pl_procStatusLocal (pl_rank, PL_STATUS_ITEM) = min (pl_itemCount, 1)
!
!
!    ...Ready!
!
!
  return
end subroutine pl_localSaveRecvItems
