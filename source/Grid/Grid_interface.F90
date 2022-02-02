!!****h* source/Grid/Grid_interface
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
!! This is the header file for the Grid unit that defines its
!! public interfaces.
!!***

!!REORDER(4): fluxBuf[XYZ], fluxCorr[XYZ], fluxOld[XYZ]

Module Grid_interface

  use Grid_ascModule, ONLY: Grid_ascStart,Grid_ascAllocMem,Grid_ascDeallocMem
  use Grid_getBlkIndexLimits_mod, ONLY: Grid_getBlkIndexLimits
  use Grid_ptDataTypes, ONLY: Grid_particle_t

  implicit none

#include "Simulation.h"
#ifdef Grid_releaseBlkPtr
! disabling drift macro expansion because it doesn't apply here, only to call sites. see: drift
#undef Grid_releaseBlkPtr
#endif

#include "constants.h"
#include "Particles.h"
!#include "GridParticles.h"

  integer,parameter :: GRID_PDE_BND_ISOLATED  = 0
  integer,parameter :: GRID_PDE_BND_PERIODIC  = 1
  integer,parameter :: GRID_PDE_BND_DIRICHLET = 2
  integer,parameter :: GRID_PDE_BND_NEUMANN   = 3
  integer,parameter :: GRID_PDE_BND_GIVENVAL  = 4
  integer,parameter :: GRID_PDE_BND_GIVENGRAD = 5

  integer,parameter :: GRID_COPYDIR_TO_VECT   = 1
  integer,parameter :: GRID_COPYDIR_FROM_VECT = 2

#include "FortranLangFeatures.fh"

  interface Grid_ascGetBlkPtr
     subroutine Grid_ascGetBlkPtr(blockID,dataPtr, gridDataStruct)
       implicit none
       integer, intent(in) :: blockID
       real, dimension(:,:,:,:), pointer :: dataPtr
       integer, optional,intent(in) :: gridDataStruct
     end subroutine Grid_ascGetBlkPtr
     subroutine Grid_ascGetBlk5Ptr(blockID,data5Ptr, gridDataStruct)
       implicit none
       integer, intent(in) :: blockID
       real, dimension(:,:,:,:,:), pointer :: data5Ptr
       integer, optional,intent(in) :: gridDataStruct
     end subroutine Grid_ascGetBlk5Ptr
  end interface Grid_ascGetBlkPtr

  interface
     subroutine Grid_ascReleaseBlkPtr(blockId, dataPtr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine Grid_ascReleaseBlkPtr
     subroutine Grid_ascReleaseBlk5Ptr(blockId, data5Ptr, gridDataStruct)
       implicit none
       integer,intent(in) :: blockId
       real, POINTER_INTENT_OUT :: data5Ptr(:,:,:,:,:)
       integer,optional, intent(in) :: gridDataStruct
     end subroutine Grid_ascReleaseBlk5Ptr
  end interface

  interface
     subroutine Grid_genGetBlkPtr(blockID,dataPtr, varDesc,dataPtr2,dataPtr3)
       implicit none
       integer, intent(in) :: blockID
       real, dimension(:,:,:,:), POINTER_INTENT_OUT :: dataPtr
       integer,intent(in),OPTIONAL :: varDesc(:)
       real, dimension(:,:,:,:), POINTER_INTENT_OUT,OPTIONAL :: dataPtr2, dataPtr3
     end subroutine Grid_genGetBlkPtr
  end interface

  interface
     subroutine Grid_genReleaseBlkPtr(blockID, dataPtr, varDesc,dataPtr2,dataPtr3)
       implicit none
       integer,intent(in) :: blockId
       real, pointer :: dataPtr(:,:,:,:)
       integer,intent(in),OPTIONAL :: varDesc(:)
       real, dimension(:,:,:,:), pointer,OPTIONAL :: dataPtr2, dataPtr3
     end subroutine Grid_genReleaseBlkPtr
  end interface

  interface
     subroutine Grid_applyBCEdge(bcType,bcDir,guard,var,dataRow,face,&
          gridDataStruct, secondCoord, thirdCoord)
       integer,intent(IN):: bcType,bcDir,guard,var,face,gridDataStruct
       real,dimension(:),intent(INOUT)::dataRow
       real,intent(IN),OPTIONAL :: secondCoord,thirdCoord
     end subroutine Grid_applyBCEdge
  end interface

  interface
     subroutine Grid_applyBCEdgeAllUnkVars(bcType,bcDir,guard,dataRow,face,&
          cellCenterSweepCoord, secondCoord,thirdCoord)
       integer,intent(IN):: bcType,bcDir,guard,face
       real,dimension(2*guard,NUNK_VARS),intent(INOUT)::dataRow
       real,intent(IN):: cellCenterSweepCoord(*), secondCoord,thirdCoord
     end subroutine Grid_applyBCEdgeAllUnkVars
  end interface

  interface Grid_computeUserVars
     subroutine Grid_computeUserVars()
     end subroutine Grid_computeUserVars
  end interface

  interface
     subroutine Grid_computeVarMean(iUnk, mean)
       implicit none
       integer, intent(in) :: iUnk
       real, intent(out) :: mean
     end subroutine Grid_computeVarMean
  end interface

  interface
     subroutine Grid_computeVarNorm (level, normType, ivar, norm, leafOnly)
       implicit none
       integer, intent(IN)  :: normType, level, ivar, leafOnly
       real, intent(OUT)    :: norm
     end subroutine Grid_computeVarNorm
  end interface

  interface
     subroutine Grid_computeVarDiff(level, gr_iRefSoln, gr_iSoln, ires)
       implicit none
       integer, intent(in)          :: level, gr_iRefSoln, gr_iSoln, ires
     end subroutine Grid_computeVarDiff
  end interface

  interface
     subroutine Grid_communicateFluxes(axis, coarse_level)
       implicit none
       integer, intent(IN)                   :: axis
       integer, intent(IN)                   :: coarse_level
     end subroutine Grid_communicateFluxes
  end interface

  interface Grid_conserveFluxes
     subroutine Grid_conserveFluxes(axis, coarse_level)
       integer, intent(IN)                   :: axis
       integer, intent(IN)                   :: coarse_level
     end subroutine Grid_conserveFluxes
  end interface

  interface Grid_dump
     subroutine Grid_dump(var,num,solnData,blockDesc,gcell)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer, intent(IN) :: num
       integer, dimension(num), intent(IN) :: var
       real,dimension(:,:,:,:),pointer     :: solnData
       type(Grid_tile_t),target, intent(in) :: blockDesc
       logical, intent(IN) :: gcell
     end subroutine Grid_dump
  end interface


  interface
     subroutine Grid_fillGuardCells( gridDataStruct,idir,&
          minLayers,eosMode,doEos,maskSize,mask,makeMaskConsistent,&
          doLogMask,selectBlockType,unitReadsMeshDataOnly)
       integer, intent(in) :: gridDataStruct
       integer, intent(in) :: idir
       integer,optional,intent(in) :: minLayers
       integer,optional,intent(in) :: eosMode
       logical,optional,intent(IN) :: doEos
       integer,optional,intent(in) :: maskSize
       logical,optional,dimension(:), intent(IN) :: mask
       logical,optional, intent(IN) :: makeMaskConsistent
       logical,optional,intent(IN) :: doLogMask
       integer,optional,intent(in) :: selectBlockType
       logical, optional, intent(IN) :: unitReadsMeshDataOnly
     end subroutine Grid_fillGuardCells
  end interface

  interface Grid_finalize
     subroutine Grid_finalize()
     end subroutine Grid_finalize
  end interface

  interface Grid_getBlkCenterCoords
     subroutine Grid_getBlkCenterCoords(blockDesc,blockCenter)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       real, dimension(MDIM), intent(out) :: blockCenter
     end subroutine Grid_getBlkCenterCoords
  end interface

  interface Grid_getBlkCornerID
     subroutine Grid_getBlkCornerID_desc(blockDesc, cornerID, stride,cornerIDHigh,inRegion)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,dimension(MDIM), intent(OUT) :: cornerID, stride
       integer,dimension(MDIM),optional,intent(OUT) :: cornerIDHigh
       logical, optional, intent(IN) :: inRegion
     end subroutine Grid_getBlkCornerID_desc
     subroutine Grid_getBlkCornerID(blockID, cornerID, stride,cornerIDHigh,inRegion)
       implicit none
       integer, intent(in) :: blockID
       integer,dimension(MDIM), intent(OUT) :: cornerID, stride
       integer,dimension(MDIM),optional,intent(OUT) :: cornerIDHigh
       logical, optional, intent(IN) :: inRegion
     end subroutine Grid_getBlkCornerID
  end interface


  interface
     subroutine Grid_getBlkIndexLimits_STANDALONE(blockId, blkLimits, blkLimitsGC,gridDataStruct)
       integer,intent(IN)                     :: blockId
       integer,dimension(2,MDIM), intent(OUT) :: blkLimits, blkLimitsGC
       integer,optional,intent(IN)            :: gridDataStruct
     end subroutine Grid_getBlkIndexLimits_STANDALONE
  end interface

  interface Grid_getCellCoords
     subroutine Grid_getCellCoords(axis, edge, level, lo, hi, coordinates)
        integer, intent(in)  :: axis
        integer, intent(in)  :: edge
        integer, intent(in)  :: level
        integer, intent(in)  :: lo(1:MDIM)
        integer, intent(in)  :: hi(1:MDIM)
        real,    intent(out) :: coordinates(:)
     end subroutine Grid_getCellCoords
  end interface

  interface Grid_getCellFaceAreas
     subroutine Grid_getCellFaceAreas(axis, level, lo, hi, areas)
        integer, intent(in)  :: axis
        integer, intent(in)  :: level
        integer, intent(in)  :: lo(1:MDIM)
        integer, intent(in)  :: hi(1:MDIM)
        real,    intent(out) :: areas(lo(IAXIS):hi(IAXIS), &
                                      lo(JAXIS):hi(JAXIS), &
                                      lo(KAXIS):hi(KAXIS))
     end subroutine Grid_getCellFaceAreas
  end interface

  interface Grid_getCellVolumes
     subroutine Grid_getCellVolumes(level, lo, hi, volumes)
        integer, intent(in)  :: level
        integer, intent(in)  :: lo(1:MDIM)
        integer, intent(in)  :: hi(1:MDIM)
        real,    intent(out) :: volumes(lo(IAXIS):hi(IAXIS), &
                                        lo(JAXIS):hi(JAXIS), &
                                        lo(KAXIS):hi(KAXIS))
     end subroutine Grid_getCellVolumes
  end interface

  interface Grid_getDeltas
     subroutine Grid_getDeltas(lev, del)
       integer, intent(in) :: lev
       real, dimension(MDIM), intent(out) :: del
     end subroutine Grid_getDeltas
  end interface

  interface
     subroutine Grid_getGlobalIndexLimits(globalIndexLimits)
       integer, dimension(MDIM), intent(out) :: globalIndexLimits
     end subroutine Grid_getGlobalIndexLimits
  end interface

  interface
     subroutine Grid_getListOfBlocks(blockType, listOfBlocks,count,refinementLevel,region_bndBox,includePartialBlocks)
       integer, intent(in) :: blockType
       integer,dimension(MAXBLOCKS),intent(out) :: listOfBlocks
       integer,intent(out) :: count
       integer,intent(IN),optional :: refinementLevel
       real, dimension(LOW:HIGH,MDIM), intent(IN), optional :: region_bndBox
       logical, intent(IN), optional :: includePartialBlocks
     end subroutine Grid_getListOfBlocks
  end interface

  interface
     subroutine Grid_getLocalNumBlks(numBlocks)
       integer,intent(out) :: numBlocks
     end subroutine Grid_getLocalNumBlks
  end interface
  
  interface
     subroutine Grid_getNumBlksFromType(blockType,numBlocks)
       integer,intent(in)  :: blockType
       integer,intent(out) :: numBlocks
     end subroutine Grid_getNumBlksFromType
  end interface

  interface Grid_getMinCellSize
     subroutine Grid_getMinCellSize(minCellSize)
       real, intent(OUT) :: minCellSize
     end subroutine Grid_getMinCellSize
  end interface

  interface
     subroutine Grid_getMinCellSizes(minCellSizes)
       real, dimension(MDIM), intent(OUT) :: minCellSizes
     end subroutine Grid_getMinCellSizes
  end interface

  interface
     subroutine Grid_subcellGeometry(nsubI,nsubJ,nsubK, &
          dvCell, dvSub, xL,xR, yL,yR, pos, blockID)
       implicit none
       integer, VALUE_INTENT(IN) :: nsubI, nsubJ, nsubK
       real,    intent(in)  :: dvCell
       real,    intent(OUT) :: dvSub(nsubI, nsubJ)
       real,OPTIONAL,intent(in) :: xL, xR
       real,OPTIONAL,intent(in) :: yL, yR
       integer,OPTIONAL, intent(in) :: blockID
       integer,OPTIONAL, intent(in) :: pos(*)
     end subroutine Grid_subcellGeometry
  end interface


  interface Grid_getSingleCellCoords
     subroutine Grid_getSingleCellCoords(ind, level,edge, coords)
       implicit none
       integer,dimension(MDIM), intent(in) :: ind
       integer, intent(in) :: level, edge
       real, dimension(MDIM), intent(out) :: coords
     end subroutine Grid_getSingleCellCoords
  end interface

  interface Grid_getSingleCellVol
     subroutine Grid_getSingleCellVol(point, level, cellvolume)
       integer, intent(in)  :: point(1:MDIM)
       integer, intent(in)  :: level
       real,    intent(out) :: cellvolume
     end subroutine Grid_getSingleCellVol
  end interface Grid_getSingleCellVol

  interface
     subroutine Grid_guardCellMaskHook(ccMask, needEos)
       implicit none
       logical,intent(INOUT) :: ccMask(*)
       logical,intent(IN)    :: needEos
     end subroutine Grid_guardCellMaskHook
  end interface

  interface Grid_init
     subroutine Grid_init()
     end subroutine Grid_init
  end interface

  interface Grid_initDomain
     subroutine Grid_initDomain( restart,particlesInitialized)
       logical, intent(IN) :: restart
       logical,intent(INOUT) :: particlesInitialized
     end subroutine Grid_initDomain
  end interface

  interface
     subroutine Grid_makeVector(vecLen,numVars,newVec,numVec,vecLastFree,copyDirection,gridDataStruct)
       implicit none
       integer, intent(in) :: vecLen
       integer, intent(in) :: numVars
       integer,intent(INOUT) :: numVec
       real, dimension(vecLen,numVars,numVec),intent(INOUT) :: newVec
       integer, optional,intent(OUT):: vecLastFree
       integer, optional,intent(in) :: copyDirection
       integer, optional,intent(in) :: gridDataStruct
     end subroutine Grid_makeVector
  end interface

  interface Grid_markRefineDerefine
     subroutine Grid_markRefineDerefine()
     end subroutine Grid_markRefineDerefine
  end interface

  interface Grid_markRefineSpecialized
     subroutine Grid_markRefineSpecialized(criterion,size,specs,lref)
       integer, intent(IN) :: criterion
       integer, intent(IN) :: size
       real,dimension(size),intent(IN) :: specs
       integer, intent(IN) ::  lref
     end subroutine Grid_markRefineSpecialized
  end interface

  interface Grid_moveParticles
     subroutine Grid_moveParticles(dataBuf, propCount, maxCount, localCount, &
       index_list, indexCount,&
       coords_in_blk)
       
       
       integer,intent(IN) :: maxCount, propCount, indexCount
       integer,intent(INOUT) :: localCount
       
       real, dimension(propCount, maxCount),intent(INOUT) :: dataBuf
       integer, dimension(indexCount),intent(IN) :: index_list
       logical, intent(IN) :: coords_in_blk
       
     end subroutine Grid_moveParticles
  end interface

  interface Grid_notifySolnDataUpdate
     subroutine Grid_notifySolnDataUpdate(gds,mask)
       implicit none
       integer,OPTIONAL,intent(in) :: gds     !"grid data struct", i.e., CENTER, CENTER_FACES, FACES
       logical,OPTIONAL,intent(in) :: mask(*) !optional mask, as for Grid_fillGuardCells
     end subroutine Grid_notifySolnDataUpdate
     subroutine Grid_notifySolnDataUpdateVlist(varList,gds)
       implicit none
       integer,intent(in)          :: varList(:)     ! list of UNK (or other?) variables
       integer,OPTIONAL,intent(in) :: gds     !"grid data struct", i.e., CENTER, CENTER_FACES, FACES
     end subroutine Grid_notifySolnDataUpdateVlist
  end interface


  interface Grid_putFluxData
     subroutine Grid_putFluxData(level,axis, pressureSlots)
       implicit none
       integer, intent(IN) :: level
       integer, intent(IN),optional :: axis
       integer, intent(IN), OPTIONAL,target :: pressureSlots(:)
     end subroutine Grid_putFluxData
     subroutine Grid_putFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo,isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,CONTIGUOUS_INTENT(in),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX,fluxBufY,fluxBufZ
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_putFluxData_block
  end interface Grid_putFluxData

  interface Grid_getFluxData ! get (corrected) flux data
     subroutine Grid_getFluxData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, axis, isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,intent(  OUT),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
       integer, intent(IN),optional :: axis
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_getFluxData_block
  end interface

  interface Grid_getFluxCorrData
     subroutine Grid_getFluxCorrData_xtra(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, fluxCorrX,fluxCorrY,fluxCorrZ, isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,CONTIGUOUS_INTENT(in) ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
       real,CONTIGUOUS_INTENT(OUT),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxCorrX,fluxCorrY,fluxCorrZ
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_getFluxCorrData_xtra
     subroutine Grid_getFluxCorrData_block(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,intent(OUT)  ,dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_getFluxCorrData_block
  end interface

  interface Grid_correctFluxData
     subroutine Grid_correctFluxData(blockDesc,fluxBufX,fluxBufY,fluxBufZ, lo, isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,CONTIGUOUS_INTENT(INOUT),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxBufX, fluxBufY, fluxBufZ
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_correctFluxData

     subroutine Grid_correctFluxData_xtra(blockDesc,scaleF,fluxBufX,fluxBufY,fluxBufZ, lo, &
                                                    scaleC,fluxOldX,fluxOldY,fluxOldZ,     &
                                                    isFluxDensity)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(in) :: blockDesc
       integer,intent(in) :: lo(3)
       real,intent(in)    :: scaleF,scaleC
       real,CONTIGUOUS_INTENT(in   ),dimension(: ,lo(1): ,lo(2): ,lo(3): ),TARGET :: fluxOldX, fluxOldY, fluxOldZ
       real,INTENT(INOUT),dimension(size(fluxOldX,1),size(fluxOldX,2),size(fluxOldX,3),size(fluxOldX,4)),TARGET :: FLUXBUFX !CAPITALIZED TO PREVENT INDEX REORDERING!
       real,INTENT(INOUT),dimension(size(fluxOldY,1),size(fluxOldY,2),size(fluxOldY,3),size(fluxOldY,4)),TARGET :: FLUXBUFY !CAPITALIZATION INTENTIONAL!
       real,INTENT(INOUT),dimension(size(fluxOldZ,1),size(fluxOldZ,2),size(fluxOldZ,3),size(fluxOldZ,4)),TARGET :: FLUXBUFZ !CAPITALIZATION INTENTIONAL!
       logical, intent(IN), OPTIONAL :: isFluxDensity(:) !maybe eliminate
     end subroutine Grid_correctFluxData_xtra
  end interface Grid_correctFluxData


  interface Grid_putLocalNumBlks
     subroutine Grid_putLocalNumBlks(numBlocks)
       integer,intent(in) :: numBlocks
     end subroutine Grid_putLocalNumBlks
  end interface

  interface Grid_restrictAllLevels
     subroutine Grid_restrictAllLevels()
     end subroutine Grid_restrictAllLevels
  end interface

  interface
     subroutine Grid_restrictByLevels( gridDataStruct, fromLevel, toLevel, checkFinestLevel,&
          maskSize,mask)
       integer, intent(in) :: gridDataStruct
       integer, intent(in) :: fromLevel, toLevel
       logical, optional,intent(in) :: checkFinestLevel
       integer, optional,intent(in) :: maskSize
       logical,dimension(*),optional,intent(in) :: mask
     end subroutine Grid_restrictByLevels
  end interface

  interface Grid_sendOutputData
     subroutine Grid_sendOutputData()
     end subroutine Grid_sendOutputData
  end interface

  interface
     subroutine Grid_setFluxHandling(handling, status)
       implicit none
       character(len=*),intent(IN) :: handling
       integer,intent(OUT),OPTIONAL :: status
     end subroutine Grid_setFluxHandling
  end interface

  interface
     subroutine Grid_setInterpValsGcell(setval)
       implicit none
       logical, intent(IN) :: setval
     end subroutine Grid_setInterpValsGcell
  end interface

  interface
     subroutine Grid_setWork(tileDesc,work,mode)
       use Grid_tile,   ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(IN) :: tileDesc
       real, intent(IN)    :: work
       integer, intent(IN),optional :: mode
     end subroutine Grid_setWork
  end interface

  interface
     subroutine Grid_setWorkBounds(pwork_bnd,lwork_bnd,bnd_type)
       implicit none
       real, intent(IN)    :: pwork_bnd,lwork_bnd
       integer,intent(IN)  :: bnd_type
     end subroutine Grid_setWorkBounds
  end interface

  interface
     subroutine Grid_setWorkDefault(pwork,lwork)
       implicit none
       real, intent(IN)    :: pwork,lwork
     end subroutine Grid_setWorkDefault
  end interface

  interface Grid_updateRefinement
     subroutine Grid_updateRefinement( nstep,time, gridChanged)
       integer, intent(in) :: nstep
       real, intent(in) :: time
       logical, intent(out), OPTIONAL :: gridChanged
     end subroutine Grid_updateRefinement
  end interface

  interface Grid_unitTest
     subroutine Grid_unitTest(fileUnit,perfect)

       integer, intent(in)           :: fileUnit ! Output to file
       logical, intent(inout)        :: perfect  ! Flag to indicate errors
     end subroutine Grid_unitTest
  end interface

  interface Grid_getGeometry
     subroutine Grid_getGeometry(geometry)
       integer, intent(OUT) :: geometry
     end subroutine Grid_getGeometry
  end interface


  interface Grid_sortParticles
     subroutine Grid_sortParticles(dataBuf,props,localNumCount,&
          elementTypes,maxPerProc,&
          elementsPerBlk, attrib1, attrib2)

       implicit none
       integer,intent(INOUT) :: localNumCount
       integer,intent(IN) :: maxPerProc, props,elementTypes
       real,intent(INOUT),dimension(props,maxPerProc) :: dataBuf
       integer,intent(OUT),dimension(MAXBLOCKS,elementTypes) :: elementsPerBlk
       integer, intent(IN) :: attrib1
       integer, optional, intent(IN) :: attrib2
     end subroutine Grid_sortParticles
  end interface

  interface Grid_countParticles
      subroutine Grid_countParticles(props, localCount,elementTypes,maxCount)

         implicit none
         integer,intent(INOUT) :: localCount
         integer,intent(IN) :: elementTypes
         integer, optional, intent(IN) :: maxCount
         integer,intent(IN) :: props
      end subroutine Grid_countParticles
  end interface

  interface Grid_countParticlesByBlock
      subroutine Grid_countParticlesByBlock(particlesPerBlock)

         implicit none
         integer,dimension(MAXBLOCKS, NPART_TYPES), intent(OUT) :: particlesPerBlock
      end subroutine Grid_countParticlesByBlock
  end interface


  interface Grid_mapMeshToParticles
     subroutine Grid_mapMeshToParticles_pc (ptContainerPos, part_props,part_blkID,&
                                         posAttrib,&
                                         numAttrib, attrib,&
                                         mapType,gridDataStruct)
       implicit none
       integer, INTENT(in) :: ptContainerPos, part_props, part_blkID
       integer, intent(IN) :: numAttrib
       integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
       integer,dimension(MDIM), intent(IN) :: posAttrib
       integer, INTENT(IN) :: mapType
       integer, optional, intent(IN) :: gridDataStruct
     end subroutine Grid_mapMeshToParticles_pc

      subroutine Grid_mapMeshToParticles (particles, part_props,part_blkID,&
                                         numParticles,&
                                         posAttrib,&
                                         numAttrib, attrib,&
                                         mapType,gridDataStruct)
       implicit none
       integer, INTENT(in) :: part_props, numParticles, part_blkID
       real, INTENT(inout),dimension(part_props,numParticles) :: particles
       integer, intent(IN) :: numAttrib
       integer, dimension(PART_ATTR_DS_SIZE,numAttrib),INTENT(in) :: attrib
       integer,dimension(MDIM), intent(IN) :: posAttrib
       integer, INTENT(IN) :: mapType
       integer, optional, intent(IN) :: gridDataStruct
     end subroutine Grid_mapMeshToParticles
  end interface

  interface
     subroutine Grid_mapParticlesToMesh (particles,part_props,numParticles,&
          maxParticlesPerProc,propPart, varGrid, mode, ptInfo)
  
       integer,intent(IN) :: numParticles, part_props,maxParticlesPerProc
       real,dimension(part_props,numParticles),intent(INOUT) :: particles
       integer, INTENT(in) :: propPart, varGrid
       integer, INTENT(in), optional :: mode
       integer, INTENT(in), optional :: ptInfo
     end subroutine Grid_mapParticlesToMesh
  end interface
  
  
  interface 
     subroutine Grid_solvePoisson (iSoln, iSrc, bcTypes, &
          bcValues, poisfact)
       implicit none
       integer, intent(in)    :: iSoln, iSrc
       integer, intent(in)    :: bcTypes(6)
       real, intent(in)       :: bcValues(2,6)
       real, intent(inout)    :: poisfact
     end subroutine Grid_solvePoisson
  end interface
  
  interface
     subroutine Grid_setSolverDbgContextInfo(component,group)
       implicit none
       integer,intent(in),OPTIONAL :: component, group
     end subroutine Grid_setSolverDbgContextInfo
  end interface

  interface
     subroutine Grid_limitAbundance(blkLimits,solnData)
       integer, dimension(2,MDIM), INTENT(in) :: blkLimits
       real, POINTER :: solnData(:,:,:,:)
     end subroutine Grid_limitAbundance
  end interface


  interface
     subroutine Grid_renormAbundance(tileDesc, tileLimits, solnData)
       use Grid_tile,   ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)         :: tileDesc
       integer,           intent(IN)         :: tileLimits(LOW:HIGH, 1:MDIM)
       real,                         pointer :: solnData(:,:,:,:)
     end subroutine Grid_renormAbundance
  end interface

  interface
     subroutine Grid_renormMassScalars(blkLimits,solnData)
       integer, intent(in), dimension(2,MDIM)::blkLimits
       real,pointer :: solnData(:,:,:,:)
     end subroutine Grid_renormMassScalars
  end interface

  interface
     subroutine Grid_conserveField ()
     end subroutine Grid_conserveField
  end interface


  interface
     subroutine Grid_pfft(direction,inArray,outArray)
       integer, intent(IN) :: direction
       real, dimension(:),intent(IN) :: inArray
       real,dimension(:), intent(OUT) :: outArray
     end subroutine Grid_pfft
  end interface

  interface
     subroutine Grid_pfftInit( ndim, needMap,globalLen, mapSize,&
          transformType, baseDatType, jProcs, kProcs, refinementLevel, region_bndBox)
       integer, intent(IN) :: ndim
       logical, intent(IN) :: needMap
       integer, dimension(MDIM), intent(IN) :: globalLen
       integer,dimension(MDIM),intent(OUT) ::  mapSize
       integer,dimension(MDIM),optional,intent(IN) :: transformType
       integer,dimension(0:MDIM),optional,intent(IN) :: baseDatType
       integer,optional,intent(IN) :: jProcs, kProcs
       integer,optional,intent(IN) :: refinementLevel
       real, dimension(LOW:HIGH,MDIM), optional, intent(IN) :: region_bndBox
     end subroutine Grid_pfftInit
  end interface

  interface
     subroutine Grid_pfftFinalize()
     end subroutine Grid_pfftFinalize
  end interface

  interface
     subroutine Grid_bcApplyToRegionSpecialized(bcType,gridDataStruct,level,&
          guard,axis,face,regionData,regionSize,mask,applied,&
          secondDir,ThirdDir,endPoints,idest)
       implicit none

       integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct, level
       integer,dimension(REGION_DIM),intent(IN) :: regionSize
       real,dimension(regionSize(BC_DIR),&
            regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),&
            regionSize(STRUCTSIZE)),intent(INOUT)::regionData
       logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
       logical, intent(OUT) :: applied
       integer,intent(IN) :: secondDir,thirdDir
       integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints
       integer,intent(IN),OPTIONAL:: idest
     end subroutine Grid_bcApplyToRegionSpecialized
  end interface

  interface
     subroutine Grid_bcApplyToRegion(bcType,gridDataStruct,level,&
          guard,axis,face,regionData,regionSize,mask,applied,&
          secondDir,ThirdDir,endPoints,idest)
       implicit none
       integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct, level
       integer,dimension(REGION_DIM),intent(IN) :: regionSize
       real,dimension(regionSize(BC_DIR),&
            regionSize(SECOND_DIR),&
            regionSize(THIRD_DIR),&
            regionSize(STRUCTSIZE)),intent(INOUT)::regionData
       logical,intent(IN),dimension(regionSize(STRUCTSIZE)):: mask
       logical, intent(OUT) :: applied
       integer,intent(IN) :: secondDir,thirdDir
       integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPoints
       integer,intent(IN),OPTIONAL:: idest
     end subroutine Grid_bcApplyToRegion
  end interface

  interface
     subroutine Grid_bcApplyToRegionMixedGds(bcType,gridDataStruct,level,&
          guard,axis,face,&
          regionDataC,regionDataFN,regionDataFT1,regionDataFT2,&
          regionSizeCtr,&
          applied,&
          tileDesc,secondDir,thirdDir,endPointsCtr,rightHanded,idest)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer, intent(IN) :: bcType,axis,face,guard,gridDataStruct,level
       integer,dimension(REGION_DIM),intent(IN) :: regionSizeCtr
       real,pointer,dimension(:,:,:,:) :: regionDataFN, regionDataFT1, regionDataFT2, regionDataC
       logical, intent(INOUT) :: applied
       type(Grid_tile_t),intent(IN) :: tileDesc
       integer,intent(IN) :: secondDir,thirdDir
       integer,intent(IN),dimension(LOW:HIGH,MDIM) :: endPointsCtr
       logical, intent(IN) :: rightHanded
       integer,intent(IN),OPTIONAL:: idest
     end subroutine Grid_bcApplyToRegionMixedGds
  end interface

  interface
     subroutine Grid_pfftGetIndexLimits(configLimits,phaseLimits)
       integer,dimension(LOW:HIGH,MDIM),intent(OUT) :: configLimits, phaseLimits
     end subroutine Grid_pfftGetIndexLimits
  end interface


  interface
     subroutine Grid_getBlkType(blockID, blkType)
       implicit none
       integer,intent(in) :: blockID
       integer,intent(out) :: blkType
     end subroutine Grid_getBlkType
  end interface

  interface Grid_pfftMapToInput
     subroutine Grid_pfftMapToInput(gridVar, pfft_inArray)
       integer,intent(IN) :: gridVar
       real, dimension(:),intent(OUT) :: pfft_inArray
     end subroutine Grid_pfftMapToInput
     subroutine Grid_pfftMapToInput3DArr(gridVar, pfft_inArray)
       integer,intent(IN) :: gridVar
       real, dimension(:,:,:),intent(OUT) :: pfft_inArray
     end subroutine Grid_pfftMapToInput3DArr
  end interface

  interface Grid_pfftMapFromOutput
     subroutine Grid_pfftMapFromOutput(gridVar, pfft_outArray)
       integer,intent(IN) :: gridVar
       real, dimension(:),intent(IN) :: pfft_outArray
     end subroutine Grid_pfftMapFromOutput
     subroutine Grid_pfftMapFromOutput3DArr(gridVar, pfft_outArray)
       integer,intent(IN) :: gridVar
       real, dimension(:,:,:),intent(IN) :: pfft_outArray
     end subroutine Grid_pfftMapFromOutput3DArr
  end interface

  interface
     subroutine Grid_outsideBoundBox(pos,bndBox,outside,Negh)
       real,dimension(MDIM),intent(IN) :: pos
       real,dimension(LOW:HIGH,MDIM),intent(IN) :: bndBox
       logical, intent(OUT) :: outside
       integer, dimension(MDIM),intent(OUT) :: Negh
     end subroutine Grid_outsideBoundBox
  end interface

  interface
     subroutine Grid_getMaxCommonRefinement(inputComm, maxRefinement)
       implicit none
       integer, intent(IN) :: inputComm
       integer, intent(OUT) :: maxRefinement
     end subroutine Grid_getMaxCommonRefinement
  end interface

  interface
     subroutine Grid_getMaxRefinement(maxRefinement, mode, scope, inputComm)
       implicit none
       integer, intent(IN), OPTIONAL :: mode, scope
       integer, intent(IN), OPTIONAL :: inputComm
       integer, intent(OUT) :: maxRefinement
     end subroutine Grid_getMaxRefinement
  end interface

  interface

     subroutine Grid_GCPutScratch(gridDataStruct,needSetup,releaseSetup,&
          &blkCount,blkList,indCount,indList, gcCnt, putData)

       integer, intent(IN) :: gridDataStruct
       logical, intent(IN) :: needSetup, releaseSetup
       integer, optional,intent(IN) :: blkCount, indCount
       integer, optional, dimension(:), intent(IN) :: blkList
       integer, optional, dimension(:), intent(IN) :: indList
       integer, optional, dimension(NDIM),intent(IN) :: gcCnt  
       logical, optional,intent(IN) :: putData

     end subroutine Grid_GCPutScratch

  end interface

  interface

     subroutine Grid_GCTransferOneBlk(mode,gridDataStruct,blkIndex,blkData)

       logical, intent(IN) :: mode
       integer, intent(IN) :: gridDataStruct, blkIndex
       real, pointer, dimension(:,:,:,:) :: blkData
       
     end subroutine Grid_GCTransferOneBlk

  end interface

  
  interface
     subroutine Grid_getNumVars(gridStruct, nVar)  
       implicit none
       integer, intent(in) :: gridStruct
       integer, intent(out) :: nVar
     end subroutine Grid_getNumVars
  end interface

  interface

     subroutine Grid_addToVar(srcVar, destVar, multFactor, reset)
       integer, intent(in) :: srcVar, destVar
       real,  intent(in) :: multFactor
       logical, intent(in) :: reset
       
     end subroutine Grid_addToVar

  end interface

  interface
     subroutine Grid_smoothVar(ivar, ivarOut, &
          solnData, lbUI,lbUJ,lbUK, &
          blklim, smoothMethod, gcLayers, blockID,&
          useMinSmoothVarVal,&
          minSmoothVarVal,&
          useMaxSmoothVarVal,&
          maxSmoothVarVal,&
          smoothCoeff )
       integer, intent(in) :: ivar
       integer, intent(in) :: ivarOut
       integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
       real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
       integer, intent(in)    :: blklim(2,MDIM)
       integer, intent(in),OPTIONAL :: smoothMethod
       integer, intent(IN),OPTIONAL :: gcLayers
       integer, intent(IN),OPTIONAL :: blockID
       logical, intent(IN),OPTIONAL :: useMinSmoothVarVal,useMaxSmoothVarVal
       real   , intent(IN),OPTIONAL :: minSmoothVarVal,maxSmoothVarVal
       real   , intent(IN),OPTIONAL :: smoothCoeff
     end subroutine Grid_smoothVar
  end interface

  interface
     subroutine Grid_primitiveToConserve(blkList,count,force)
       integer,intent(IN) :: count
       integer,dimension(count),intent(IN) :: blkList 
       logical,intent(IN) :: force
     end subroutine Grid_primitiveToConserve
  end interface

  interface
     subroutine Grid_conserveToPrimitive(blkList,count,allCells,force)
       integer,intent(IN) :: count
       integer,dimension(count),intent(IN) :: blkList 
       logical,intent(IN) :: allCells,force
     end subroutine Grid_conserveToPrimitive
  end interface

  interface

     subroutine Grid_getDomainBoundBox(boundBox)
       real,dimension(LOW:HIGH,MDIM), intent(OUT) :: boundBox
     end subroutine Grid_getDomainBoundBox
  end interface

  interface
     subroutine Grid_getDomainBC(boundary)
       integer,dimension(LOW:HIGH,MDIM), intent(OUT) :: boundary
     end subroutine Grid_getDomainBC
  end interface

  interface
     subroutine Grid_parseNonRep(strlwr, nonrep, idx)
      implicit none
      character(*), intent(in) :: strlwr
      integer, intent(out) :: nonrep, idx
     end subroutine
  end interface
  
  interface
     subroutine Grid_formatNonRep(nonrep, idx, str)
        implicit none
        integer, intent(in) :: nonrep, idx
        character(*), intent(out) :: str
     end subroutine
  end interface
  
  interface
     subroutine Grid_getVarNonRep(mapblock, var, nonrep, locidx)
       implicit none
       integer, intent(in) :: mapblock, var
       integer, intent(out) :: nonrep
       integer, intent(out), optional :: locidx
     end subroutine
  end interface

  interface Grid_getBlkIDFromPos

     subroutine Grid_getBlkNeighBlkIDFromPos(blockDesc, pos, neghDir, ansBlockID, ansProcID)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN) :: blockDesc
       real, dimension(MDIM), intent(IN) :: pos
       integer, dimension(MDIM), intent(IN) :: neghDir
       integer, intent(OUT) :: ansBlockID, ansProcID
     end subroutine Grid_getBlkNeighBlkIDFromPos

     subroutine Grid_getBlkIDFromPos(pos, ansBlockID, ansProcID, comm)
       implicit none
       real, dimension(1:MDIM), intent(IN) :: pos
       integer, intent(OUT) :: ansBlockID, ansProcID
       integer,optional,intent(IN) :: comm
     end subroutine Grid_getBlkIDFromPos
  end interface

  interface Grid_getLocalBlkIDFromPos
     subroutine Grid_getLocalBlkIDFromPosSimple(pos,ansBlockID, ansProcID,blkList, blkCount,blockType)
       implicit none
       real, dimension(1:MDIM), intent(IN) :: pos
       integer, intent(OUT) :: ansBlockID, ansProcID
       integer, OPTIONAL,intent(IN)  :: blkCount
       integer, OPTIONAL,dimension(:),intent(IN),target :: blkList
       integer, OPTIONAL, intent(IN) :: blockType
     end subroutine Grid_getLocalBlkIDFromPosSimple
  end interface

  interface
     subroutine Grid_sbBroadcastParticles()
     end subroutine Grid_sbBroadcastParticles
  end interface

  interface
     subroutine Grid_sbCreateGroups()
     end subroutine Grid_sbCreateGroups
  end interface

  interface
     subroutine Grid_sbSelectMaster()
     end subroutine Grid_sbSelectMaster
  end interface

!  interface
!     subroutine Grid_updateSolidBodyForces()
!     end subroutine Grid_updateSolidBodyForces
!  end interface

  interface
     subroutine Grid_updateSolidBodyForces(blkID,particleData)
      implicit none
      integer, intent(in) :: blkID
      real, intent(inout) :: particleData(NPART_PROPS)
     end subroutine Grid_updateSolidBodyForces
  end interface


  interface
     subroutine Grid_solidBodyUnitTest(fileUnit, perfect)
       integer, intent(IN) :: fileUnit 
       logical, intent(INOUT) :: perfect
     end subroutine Grid_solidBodyUnitTest
  end interface

  interface
     subroutine Grid_receiveInputData(localNumBlocks, alnblocks, xx)
       implicit none
       integer, intent(IN) :: localNumBlocks, alnblocks, xx
     end subroutine Grid_receiveInputData
  end interface

  interface
    subroutine Grid_getNeighProcList(includeMyProc, neighProcList, numNeigh)
      implicit none
      logical, intent(IN) :: includeMyProc
      integer, dimension(:), pointer :: neighProcList
      integer, intent(OUT) :: numNeigh
    end subroutine Grid_getNeighProcList
  end interface

  interface
     subroutine Grid_getBlkNeighLevels(blockID, levels, trackBdry)
       implicit none
       integer,intent(in)  :: blockID
       integer,intent(OUT) :: levels(-1:1, -K2D:K2D , -K3D:K3D)
       logical,intent(in),OPTIONAL :: trackBdry
     end subroutine Grid_getBlkNeighLevels
  end interface

  interface
     subroutine Grid_coordTransfm(x,y,z, xout,yout,zout, geometryIn,geometryOut, ndim, velI,velJ,velK,velIOut,velJOut,velKOut)
       implicit none
       real,intent(IN) :: x,y,z
       real,intent(OUT) :: xout,yout,zout
       integer,OPTIONAL,intent(IN) :: geometryIn, geometryOut
       integer,OPTIONAL,intent(IN) :: ndim
       real,OPTIONAL,intent(IN) :: velI,velJ,velK
       real,OPTIONAL,intent(OUT) :: velIOut,velJOut,velKOut
     end subroutine Grid_coordTransfm
  end interface


  interface
     subroutine Grid_getTileIterator(itor, nodetype, level, tiling, tileSize, nthreads)
       use Grid_iterator, ONLY : Grid_iterator_t
       implicit none
       type(Grid_iterator_t), intent(OUT)          :: itor
       integer,               intent(IN)           :: nodetype
       integer,               intent(IN), optional :: level
       logical,               intent(IN), optional :: tiling
       integer,               intent(IN), optional :: tileSize(1:MDIM)
       integer,               intent(IN), optional :: nthreads
     end subroutine Grid_getTileIterator
  end interface

  interface
     subroutine Grid_releaseTileIterator(itor)
       use Grid_iterator, ONLY : Grid_iterator_t
       implicit none
       type(Grid_iterator_t), intent(INOUT) :: itor
     end subroutine Grid_releaseTileIterator
  end interface

  interface
     logical function Grid_properTilingWanted()
       implicit none
     end function Grid_properTilingWanted
  end interface

  interface
     subroutine Grid_zeroFluxData()
       implicit none
     end subroutine Grid_zeroFluxData
  end interface

  interface
     subroutine Grid_addFineToFluxRegister(level, isDensity, coefficient, &
                                           zeroFullRegister)
       implicit none
       integer, intent(IN)           :: level
       logical, intent(IN), optional :: isDensity(:)
       real,    intent(IN), optional :: coefficient
       logical, intent(IN), optional :: zeroFullRegister
     end subroutine Grid_addFineToFluxRegister
  end interface

  interface
     subroutine Grid_addCoarseToFluxRegister(level, isDensity, coefficient, &
                                             zeroFullRegister)
       implicit none
       integer, intent(IN)           :: level
       logical, intent(IN), optional :: isDensity(:)
       real,    intent(IN), optional :: coefficient
       logical, intent(IN), optional :: zeroFullRegister
     end subroutine Grid_addCoarseToFluxRegister
  end interface

end Module Grid_interface

