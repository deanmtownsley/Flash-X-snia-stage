!!****if* source/flashUtilities/Pipeline/PipelineMain/Pipeline_localFlush
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
!!  Pipeline_localFlush
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localFlush (logical, intent (in) :: fullestChannelOnly)
!!
!! DESCRIPTION
!!
!!  Flushes the local pipeline processor. When calling this routine, either
!!  only the fullest channel will be sent (fullestChannelOnly = .true.) or
!!  all channels with non-empty sending buffers are forced to send.
!!
!! ARGUMENTS
!!
!!  fullestChannelOnly : if true, only the fullest channel will send
!!
!! NOTES
!!
!!  All processors need to call this routine in order to move items through the
!!  pipeline, especially near the end of the application.
!!
!!***

subroutine Pipeline_localFlush (fullestChannelOnly)

  use Pipeline_data, ONLY : pl_numChannels, &
                            pl_sendCount,   &
                            pl_sendRequest

  use pl_interface,  ONLY : pl_localPostSendMsg

  implicit none

  include "Flashx_mpi.h"

  logical, intent (in) :: fullestChannelOnly

  logical :: checkChannel
  logical :: sendChannel

  integer :: bufSize
  integer :: channel
  integer :: fullestChannel
  integer :: largestBufSize
!
!
!     ...Take action only if channels are present.
!
!
  if (pl_numChannels > 0) then

      if (fullestChannelOnly) then
!
!
!     ...Send only fullest non-empty buffer.
!
!
          if (any (pl_sendCount (:) > 0)) then

              largestBufSize = -1
              fullestChannel = -1

              do channel = 1, pl_numChannels
                 bufSize = pl_sendCount (channel)

                 checkChannel = pl_sendRequest (channel) == MPI_REQUEST_NULL .and. bufSize  > 0

                 if (checkChannel .and. bufSize > largestBufSize) then
                     fullestChannel = channel
                     largestBufSize = bufSize
                 end if
              end do

              if (fullestChannel >= 1 .and. fullestChannel <= pl_numChannels) then
                  call pl_localPostSendMsg (fullestChannel)
              end if

          end if

      else
!
!
!     ...Force sending non-empty buffers on all local channels.
!
!
          do channel = 1, pl_numChannels
             sendChannel = pl_sendRequest (channel) == MPI_REQUEST_NULL .and. pl_sendCount (channel)  > 0
             if (sendChannel) then
                 call pl_localPostSendMsg (channel)
             end if
          end do
      end if
  end if
!
!
!    ...Ready!
!
!
  return
end subroutine Pipeline_localFlush
