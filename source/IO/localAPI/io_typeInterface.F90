!!****ih* source/IO/localAPI/io_typeInterface
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
!!  io_typeInterface
!!
!! SYNOPSIS
!!  use io_typeInterface
!!
!! DESCRIPTION
!!
!! This is an interface module for Fortran subroutines that
!! are used in derived datatype I/O.
!!
!!***

#include "constants.h"
#include "Simulation.h"
#include "io_flash.h"

module io_typeInterface

  interface
     subroutine io_xfer_mesh_data(fileID, fileFmt, fileType, &
     libType, xferType, localNumBlocks, globalBlockOffset)
       use io_intfTypesModule, ONLY : io_fileID_t
       implicit none
       integer(io_fileID_t),intent(IN) :: fileID
       integer, intent(IN) :: fileFmt, fileType, libType, &
            xferType, localNumBlocks, globalBlockOffset
     end subroutine io_xfer_mesh_data
  end interface

  interface
     subroutine io_xfer_tree_data(tree_data, fileID, libType, xferType, &
          localNumBlocksIn, localOffsetIn, presentDims, numFileBlks)
       use IO_data, ONLY : tree_data_t
       use io_intfTypesModule, ONLY : io_fileID_t
       implicit none
       type(tree_data_t), intent(INOUT) :: tree_data
       integer(io_fileID_t),intent(IN) :: fileID
       integer, intent(IN) :: libType, xferType, &
            localNumBlocksIn, localOffsetIn, presentDims, numFileBlks
     end subroutine io_xfer_tree_data
  end interface

  interface
     subroutine io_create_grid_header(myPE, fileID, fileFmt, fileType, &
          libType, dataFloatingPointType, metadataFloatingPointType)
       use io_intfTypesModule, ONLY : io_fileID_t
       implicit none
       integer(io_fileID_t),intent(IN) :: fileID
       integer, intent(IN) :: myPE, fileFmt, fileType, &
            libType, dataFloatingPointType, metadataFloatingPointType
     end subroutine io_create_grid_header
  end interface

  interface
     subroutine io_getZeroBasedVarInfo(fileType, gridDataStruct, numGridVars, &
          numOutputGridVars, gridVarOffsets, gridVarLabels)
       implicit none
       integer, intent(IN) :: fileType, gridDataStruct
       integer, intent(OUT) :: numGridVars, numOutputGridVars 
       integer, dimension(MAX_MESH_VAR), intent(OUT) :: gridVarOffsets
       character (len=MAX_STRING_LENGTH), dimension(MAX_MESH_VAR), intent(OUT) :: &
            gridVarLabels
     end subroutine io_getZeroBasedVarInfo
  end interface

  interface
     subroutine io_getZeroBasedBlkSubarray(gridDataStruct, blockInnerSize, &
          blockOuterSize, blockInnerOffset)
       implicit none
       integer, intent(IN) :: gridDataStruct
       integer, dimension(MDIM), intent(OUT) :: blockInnerSize, &
            blockOuterSize, blockInnerOffset
     end subroutine io_getZeroBasedBlkSubarray
  end interface

  interface
     subroutine io_do_xfer(xferType, gridStruct, dataset, doXfer)
       implicit none
       integer, intent(IN) :: xferType, gridStruct
       character(len=*), intent(IN) :: dataset
       logical, intent(OUT) :: doXfer
     end subroutine io_do_xfer
  end interface
end module io_typeInterface
