!!****ih* source/Grid/localAPI/gr_interface
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
!!  gr_interface
!!
!! SYNOPSIS
!!
!!  use gr_interface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms private to the GridMain subunit.
!!
!!  Currently, these are mostly subprograms used in the PARAMESH3/4 Grid implementations.
!!
!!***

! Modification history:
!     Created gr_flashWithPM3_interfaces          February 2007  KW
!     Renamed gr_pm3Interface                     January  2008  AD
!     Renamed gr_interface and moved to localAPI  June     2008  KW
!     Added gr_findMean, alphabetized functions   July     2008  LBR
!     Added gr_findWhichChildren,gr_findAllNeghID, 
!           gr_checkGridState                     Nov      2008  CD
!     Added gr_updateRefinement                   December 2008  KW
!     Modified gr_setGcFillNLayers                June     2009  KW
!     Added gr_createBlock                        October  2012  KW

module gr_interface
#include "constants.h"
#include "Simulation.h"
  implicit none

  interface
     subroutine gr_createBlock(blockImin, blockImax, &
          blockJmin, blockJmax, blockKmin, blockKmax,blockID)
       implicit none
       real,intent(IN) :: blockImin, blockImax, blockJmin, blockJmax, blockKmin, blockKmax
       integer,intent(IN) :: blockID
     end subroutine gr_createBlock
  end interface

  interface
     integer function gr_extractBCForDirection(packedBCs,axis,leftOrRight)
     ! implementation in GridMain/paramesh/gr_packBCs.F90
       implicit none
       integer,intent(IN) :: packedBCs,axis,leftOrRight
     end function gr_extractBCForDirection
  end interface

  interface
     subroutine gr_findBlock(blkList,blkCount,pos,blockID)
       implicit none
       integer,intent(IN) :: blkCount
       integer,dimension(blkCount),intent(IN) :: blkList
       real,dimension(MDIM),intent(IN) :: pos
       integer,intent(INOUT) :: blockID
     end subroutine gr_findBlock
  end interface

  interface
     subroutine gr_findMean(iSrc, iType, bGuardcell, mean)
       implicit none
       integer, intent(in) :: iSrc, iType
       logical, intent(in) :: bGuardcell
       real, intent(out) :: mean
     end subroutine gr_findMean
  end interface

  interface
     subroutine gr_findWhichChild(pos,bndBox,negh,whichChild)
       implicit none
       real,dimension(MDIM), intent(IN) :: pos
       real,dimension(LOW:HIGH,MDIM),intent(IN) :: bndBox
       integer, dimension(MDIM),intent(IN) :: negh
       integer, intent(OUT) :: whichChild
     end subroutine gr_findWhichChild
  end interface

  interface gr_findNeghID
     subroutine gr_findNeghID(blockID,pos,negh,neghID)
       implicit none
       integer, intent(IN) :: blockID
       real,dimension(MDIM), intent(IN) :: pos
       integer,dimension(MDIM),intent(IN) :: negh
       integer,dimension(BLKNO:PROCNO),intent(OUT) :: neghID
     end subroutine gr_findNeghID
  end interface gr_findNeghID

  interface
     subroutine gr_getBlkHandle(remoteBlockID, pe, blockHandle)
     ! implementation in GridMain/paramesh
       implicit none
       integer, intent(in) :: remoteBlockID, pe
       integer, intent(INOUT) :: blockHandle
     end subroutine gr_getBlkHandle
  end interface

  interface
     integer function gr_packBCs(bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight)
     ! implementation in GridMain/AMR/Paramesh4
       implicit none
       integer,intent(IN) :: bcILeft, bcIRight, bcJLeft, bcJRight, bcKLeft, bcKRight
     end function gr_packBCs
  end interface

  interface
     logical function gr_pmIoTreeMetadataIsValid()
     ! implementation in GridMain/AMR/Paramesh4
       implicit none
     end function gr_pmIoTreeMetadataIsValid
  end interface
  
  interface
     subroutine gr_setGcFillNLayers(layers, idir, guard, minLayers, returnLayers)
     ! implementation in GridMain/AMR
       implicit none
       integer,dimension(MDIM), intent(OUT) :: layers
       integer, intent(IN)  :: idir, guard
       integer, intent(IN),OPTIONAL  :: minLayers
       integer,intent(OUT),OPTIONAL  :: returnLayers(MDIM)
     end subroutine gr_setGcFillNLayers
  end interface

  interface
     subroutine gr_setMasks_gen(gridDataStruct,maskSize,mask, gcell_on_cc,gcell_on_fc,enableMaskedGCFill)
       implicit none
       integer, intent(in) :: gridDataStruct
       integer, intent(in) :: maskSize
       logical,dimension(maskSize),intent(in) :: mask
       logical, intent(INOUT)       :: gcell_on_cc(NUNK_VARS)
       logical, intent(inout),OPTIONAL :: gcell_on_fc(MDIM,NFACE_VARS)
       logical, intent(in),OPTIONAL :: enableMaskedGCFill
     end subroutine gr_setMasks_gen
  end interface

  interface
     subroutine gr_makeMaskConsistent_gen(gridDataStruct,eosMode,needEos,gcell_on_cc)
       implicit none
       integer,intent(IN) :: gridDataStruct
       integer,intent(IN) :: eosMode
       logical,intent(INOUT) :: needEos
       logical,intent(INOUT) :: gcell_on_cc(NUNK_VARS)
     end subroutine gr_makeMaskConsistent_gen
  end interface

  interface
     subroutine gr_neghAtSameLevel(blockID,atSameLevel)
       integer,intent(IN) :: blockID
       logical,dimension(LEFT_EDGE:RIGHT_EDGE,&
            LEFT_EDGE:RIGHT_EDGE,&
            LEFT_EDGE:RIGHT_EDGE),intent(OUT) :: atSameLevel
       
     end subroutine gr_neghAtSameLevel
  end interface

  interface
     subroutine gr_sanitizeDataAfterInterp(ntype, info, layers)
       implicit none
       integer, intent(IN) :: ntype
       character(len=*), intent(IN) :: info
       integer, dimension(MDIM), intent(IN):: layers
     end subroutine gr_sanitizeDataAfterInterp
     subroutine gr_sanitizeDataAfterInterp_blklst(blkList,count, info, layers)
       integer,intent(IN) :: count
       integer, dimension(count), intent(IN) :: blkList
       character(len=*), intent(IN) :: info
       integer, dimension(MDIM), intent(IN):: layers
     end subroutine gr_sanitizeDataAfterInterp_blklst
  end interface


  interface
     subroutine gr_findWhichChildren(numNegh,Negh,whichChildren)
       integer,intent(IN) :: numNegh
       integer, dimension(MDIM),intent(IN) :: Negh
       integer, intent(OUT) :: whichChildren(numNegh)
     end subroutine gr_findWhichChildren
  end interface


  interface
     subroutine gr_findAllNeghID(blockID,surrBlksSummary)
       use gr_interfaceTypeDecl, ONLY: AllBlockRegions_t
       integer, intent(IN) :: blockID
       type (AllBlockRegions_t), intent(OUT) :: surrBlksSummary
     end subroutine gr_findAllNeghID
  end interface


  interface
     subroutine gr_checkGridState()
       implicit none
     end subroutine gr_checkGridState
  end interface

  interface
     subroutine gr_estimateBlkError(error, tileDesc, iref, refine_filter)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       real,intent(INOUT) :: error
       type(Grid_tile_t),intent(IN) :: tileDesc
       integer, intent(IN) :: iref
       real, intent(IN) ::  refine_filter
     end subroutine gr_estimateBlkError
  end interface
  
  interface
     subroutine gr_markRefineDerefine(error, refine_cutoff,derefine_cutoff)
       implicit none
       real, intent(IN) :: refine_cutoff, derefine_cutoff
       real, intent(IN) :: error(MAXBLOCKS)
     end subroutine gr_markRefineDerefine
  end interface
  
  interface
     subroutine gr_updateRefinement( gridChanged,force_rebalance)
       implicit none
       logical, intent(out),OPTIONAL :: gridChanged
       logical, intent(in),OPTIONAL  :: force_rebalance
     end subroutine gr_updateRefinement
  end interface
  
  interface
     subroutine gr_GCAllocScratch(gridDataStruct,blkCnt,blkList,&
          indCnt,indList, gcCnt)
       integer, intent(IN) :: gridDataStruct, blkCnt, indCnt
       integer, dimension(blkCnt), intent(IN) :: blkList
       integer, dimension(indCnt), intent(IN) :: indList
       integer, dimension(NDIM), intent(IN) :: gcCnt
     end subroutine gr_GCAllocScratch
  end interface

  interface
     subroutine gr_GCReleaseScratch(gridDataStruct)
       integer, intent(IN) :: gridDataStruct
     end subroutine gr_GCReleaseScratch
  end interface

  interface
     subroutine gr_GCTransferOneBlk(mode,indCnt,indList,offset,&
          blkLimits,blkLimitsGC,&
          flatArray,blkArray)
       integer, intent(IN) :: indCnt,offset
       logical, intent(IN) :: mode
       integer,dimension(indCnt),intent(IN) :: indList
       integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
       real, pointer, dimension(:) :: flatArray
       real, pointer, dimension(:,:,:,:) :: blkArray
     end subroutine gr_GCTransferOneBlk
  end interface

  interface
     subroutine gr_setBlockType(blockID,type)
       integer, intent(IN) :: blockID, type
     end subroutine gr_setBlockType
  end interface

  interface
     subroutine gr_xyzToBlockLevel(lev, xyz, ijk)
       integer, intent(in) :: lev
       real, intent(in) :: xyz(NDIM)
       integer, intent(out) :: ijk(NDIM)
     end subroutine gr_xyzToBlockLevel
  end interface

  interface
     Subroutine gr_xyzToBlock(xyz, procID, blkID)
       real, dimension(MDIM),intent(IN) :: xyz
       integer, intent(OUT) :: procID
       integer, intent(OUT) :: blkID
     End Subroutine gr_xyzToBlock
  end interface


  interface
    subroutine gr_getRegionDataCoordinates(level, gridDataStruct, &
                                           axis, axis2, axis3, &
                                           regionSize, endPoints, &
                                           coordinates)
      implicit none
      integer, intent(IN)  :: level
      integer, intent(IN)  :: gridDataStruct
      integer, intent(IN)  :: axis
      integer, intent(IN)  :: axis2
      integer, intent(IN)  :: axis3
      integer, intent(IN)  :: regionSize(REGION_DIM)
      integer, intent(IN)  :: endPoints(LOW:HIGH, 1:MDIM)
      real,    intent(OUT) :: coordinates(regionSize(BC_DIR),     &
                                          regionSize(SECOND_DIR), &
                                          regionSize(THIRD_DIR),  &
                                          1:MDIM)
    end subroutine gr_getRegionDataCoordinates
  end interface

end module gr_interface
