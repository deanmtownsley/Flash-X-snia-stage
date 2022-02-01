!!****f* source/flashUtilities/Pipeline/Pipeline_globalCheckStatus
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
!!  Pipeline_globalCheckStatus
!!
!! SYNOPSIS
!! 
!!  call Pipeline_globalCheckStatus (logical, intent (out) :: finishedComm,
!!                                   logical, intent (out) :: emptyRecvBuffers,
!!                                   logical, intent (out) :: emptyItemBuffers)
!!
!! DESCRIPTION
!!
!!  Returns the global status of the pipeline. This is a blocking operation that
!!  currently evaluates the following statuses: 1) has all communication inside
!!  the pipeline been performed?, 2) are all receive buffers on all channels on
!!  all processors empty? and 3) are all item buffers on all channels empty?.
!!
!!  Procedure:
!!
!!  The status vector on all processors are all summation reduced and redistributed
!!  on all processors. On each processor the communication entry contains a positive
!!  integer +n for all n sends posted and a collection of negative integers {-m} for
!!  all channel processor ID's on which -m receives completed. Reduction under summation
!!  will give the overall active sends that have not completed. All entries in the
!!  reduced status vector must hence be >= 0. A negative entry indicates a bug someplace.
!!
!!  The receive and item buffer information in the status vector is non-overlapping
!!  and specific to each processor (with entries of +1 or 0). Hence reduction under
!!  summation and redistribution will not change the info contained in the status vector.
!!
!! ARGUMENTS
!!
!!  finishedComm     : is set true, if no more communication is scheduled globally
!!  emptyRecvBuffers : is set true, if all receive buffers on all processors are empty
!!  emptyItemBuffers : is set true, if all item buffers on all processors are empty
!!
!! NOTES
!!
!!  This is a collective operation, hence all processors must call it.
!!
!!***

subroutine Pipeline_globalCheckStatus (finishedComm, emptyRecvBuffers, emptyItemBuffers)

  implicit none

  logical, intent (out) :: finishedComm
  logical, intent (out) :: emptyRecvBuffers
  logical, intent (out) :: emptyItemBuffers

  finishedComm     = .false.
  emptyRecvBuffers = .false.
  emptyItemBuffers = .false.

  return
end subroutine Pipeline_globalCheckStatus
