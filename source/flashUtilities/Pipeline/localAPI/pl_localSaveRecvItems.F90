!!****if* source/flashUtilities/Pipeline/localAPI/pl_localSaveRecvItems
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
!!                              logical, intent (out) :: isSaved)
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

  implicit none

  integer, intent (in)  :: channel
  logical, intent (out) :: isSaved

  isSaved = .false.

  return
end subroutine pl_localSaveRecvItems
