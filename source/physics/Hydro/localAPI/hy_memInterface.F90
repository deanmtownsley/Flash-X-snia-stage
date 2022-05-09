!!****ih* source/physics/Hydro/localAPI/hy_memInterface
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
!!  hy_memInterface
!!
!! SYNOPSIS
!!  use hy_memInterface
!!
!! DESCRIPTION
!!  This is a local interface for memory management routines.
!!
!!***


Module hy_memInterface

  use hy_memScratchData, ONLY: hy_memScratchInitialize, hy_memAllocScratch, hy_memDeallocScratch

  implicit none

#include "FortranLangFeatures.fh"

  interface hy_memGetBlkPtr
     subroutine hy_memGetBlkPtr_blkid(blockID,dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlkPtr_blkid
     subroutine hy_memGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlk5Ptr
     subroutine hy_memGetBlkPtr_desc(tileDesc,dataPtr, gridDataStruct)
       use Grid_tile, ONLY : Grid_tile_t 
       implicit none
       type(Grid_tile_t), intent(IN) :: tileDesc
       real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memGetBlkPtr_desc
  end interface

  interface hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlkPtr(blockId, dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlkPtr
     subroutine hy_memReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlk5Ptr
     subroutine hy_memReleaseBlkPtr_desc(tileDesc, dataPtr, gridDataStruct)
       use Grid_tile, ONLY : Grid_tile_t 
       implicit none
       type(Grid_tile_t), intent(IN) :: tileDesc
       real, POINTER_INTENT_OUT :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine hy_memReleaseBlkPtr_desc
  end interface

End Module hy_memInterface
