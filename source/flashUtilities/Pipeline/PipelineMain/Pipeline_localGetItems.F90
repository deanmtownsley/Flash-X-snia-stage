!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_localGetItems
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
!!  Pipeline_localGetItems
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localGetItems (real,    intent (inout) :: userArray (:,:),
!!                               integer, intent (out)   :: userCount)
!!
!! DESCRIPTION
!!
!!  The routine copies items currently sitting in the item buffer of the local
!!  processor to the user supplied array. Once the items in the buffer have been
!!  copied, the item buffer counter is reset. This routine constitutes thus the
!!  local retriever of the pipeline.
!!
!! ARGUMENTS
!!
!!  userArray : the user supplied array to hold the copied items
!!  userCount : the number of items copied
!!
!! NOTES
!!
!!  The user supplied array must be able to hold all elements of each item. If this
!!  is not the case, the program aborts.
!!
!!***

subroutine Pipeline_localGetItems (userArray, userCount)

  use Pipeline_data,     ONLY : pl_doLog,            &
                                pl_itemBuf,          &
                                pl_itemCount,        &
                                pl_itemSize,         &
                                pl_logUnit,          &
                                pl_procStatusLocal,  &
                                pl_rank

  use Driver_interface,  ONLY : Driver_abort

#include "Pipeline.h"

  implicit none

  real,    intent (inout) :: userArray (:,:)
  integer, intent (out)   :: userCount

  integer :: firstItem
  integer :: itemsToCopy
!
!
!     ...Check, if number of elements per item fits into the user array.
!
!
  if (size (userArray,1) < pl_itemSize) then
      call Driver_abort ('[Pipeline_localGetItems] ERROR: User array cannot hold item elements!')
  end if
!
!
!     ...Decide how many items to copy. Always copy from the end of the items buffer, to
!        avoid necessary reshuffling of remaining items in buffer.
!
!
  itemsToCopy = min (pl_itemCount, size (userArray,2))

  if (itemsToCopy > 0) then
      firstItem = pl_itemCount - itemsToCopy + 1
      userArray (:,1:itemsToCopy) = pl_itemBuf (:,firstItem:pl_itemCount)
      pl_itemCount = pl_itemCount - itemsToCopy
  end if

  userCount = itemsToCopy

  if (pl_doLog) then
      write (pl_logUnit,'(a,i6,a)') ' Copied                              ',userCount, &
                                    ' buffer items to user supplied array '
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
end subroutine Pipeline_localGetItems
