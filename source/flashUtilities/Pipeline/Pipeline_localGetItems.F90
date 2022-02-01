!!****f* source/flashUtilities/Pipeline/Pipeline_localGetItems
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

  implicit none

  real,    intent (inout) :: userArray (:,:)
  integer, intent (out)   :: userCount

  userCount = 0
  userArray (:,:) = 0.0

  return
end subroutine Pipeline_localGetItems
