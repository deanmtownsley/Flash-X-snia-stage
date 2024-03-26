!!****ih* source/Grid/GridMain/AMR/Paramesh4/Grid_tile
!! NOTICE
!!  Copyright 2024 UChicago Argonne, LLC and contributors
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
!!
!!
!!****

!!REORDER(4): dataPtr

#include "FortranLangFeatures.fh"
#include "constants.h"
#include "Simulation.h"

module Grid_tile
    implicit none

    private

    type, public :: Grid_tile_t
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: grownLimits(LOW:HIGH, MDIM)
        integer :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: enclosingBlock
        procedure, public :: fillTileCInfo
    end type Grid_tile_t

contains

    subroutine deltas(this, dx)
        use Grid_data, ONLY: gr_delta

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        dx(1:MDIM) = gr_delta(1:MDIM, this%level)
    end subroutine deltas

    subroutine boundBox(this, box)
        use Grid_data,        ONLY: gr_delta
        use tree,             ONLY : bnd_box
        use Driver_interface, ONLY : Driver_abort

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)
        integer         :: i
  
        if (this%id <= 0) then
           print *, "blockId = ", this%id
           call Driver_abort("[boundBox] blockId out of bounds")
        end if
  
        box = bnd_box(:, :, this%id)
        do i = 1,NDIM
           box(LOW,i)  = box(LOW,i)  + (this%limits(LOW,i)  - this%blkLimitsGC(LOW,i)  - NGUARD) * gr_delta(i, this%level)
           box(HIGH,i) = box(HIGH,i) + (this%limits(HIGH,i) - this%blkLimitsGC(HIGH,i) + NGUARD) * gr_delta(i, this%level)
        end do
    end subroutine boundBox

    subroutine physicalSize(this, tileSize) 
        use Grid_data, ONLY: gr_delta
        use tree, ONLY : bsize

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM)
        integer         :: i
      
        tileSize = bsize(:, this%id)
        do i = 1,NDIM
           tileSize(i) = tileSize(i) + gr_delta(i, this%level) * &
                ( this%limits(HIGH,i) - this%limits(LOW,i) + this%blkLimitsGC(LOW,i) - this%blkLimitsGC(HIGH,i) + 2*NGUARD )
        end do

    end subroutine physicalSize

    subroutine faceBCs(this, faces, onBoundary)
        use tree,         ONLY : bnd_box
        use Grid_data,    ONLY : gr_domainBC, &
                                 gr_globalDomain, &
                                 gr_delta

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)

        real    :: deltas(1:MDIM)
        integer :: axis, face

        deltas(1:MDIM) = gr_delta(1:MDIM, this%level)

        do    axis = 1, MDIM
           do face = LOW, HIGH
              faces(face, axis) = NOT_BOUNDARY
              if (present (onBoundary)) then
                 onBoundary(face,axis) = NOT_BOUNDARY
              end if

              if (almostEqual(bnd_box(face, axis, this%id), &
                              gr_globalDomain(face, axis), &
                              deltas(axis))) then
                 if (gr_domainBC(face, axis) .NE. PERIODIC) &
                      faces(face,axis) = gr_domainBC(face,axis)
                 if (present (onBoundary)) then
                    onBoundary(face,axis) = gr_domainBC(face,axis)
                 end if
              end if

           end do
        end do

    contains
       logical function almostEqual(x, y, dx)
          real, intent(IN) :: x
          real, intent(IN) :: y
          real, intent(IN) :: dx

          almostEqual = (ABS(x-y) <= (0.01 * dx))
       end function almostEqual
    end subroutine faceBCs

    pure function enclosingBlock(this)

        class(Grid_tile_t), intent(IN)  :: this
        type(Grid_tile_t)               :: enclosingBlock

        integer         :: i
        integer,parameter,dimension(MDIM) :: guardv = (/(NGUARD,i=1,NDIM),(0,i=NDIM+1,MDIM)/)

        select type(this)
        type is (Grid_tile_t)
           enclosingBlock = this
        class default
           enclosingBlock%id          = this%id
           enclosingBlock%cid         = this%cid
           enclosingBlock%stride      = this%stride
           enclosingBlock%level       = this%level
           enclosingBlock%limits      = this%limits
           enclosingBlock%blkLimitsGC = this%blkLimitsGC
        end select

        enclosingBlock%limits(LOW, 1:NDIM) = &
             this%blkLimitsGC(LOW,  1:NDIM) + guardv(1:NDIM)
        enclosingBlock%limits(HIGH,1:NDIM) = &
             this%blkLimitsGC(HIGH, 1:NDIM) - guardv(1:NDIM)

        enclosingBlock%grownLimits(:,:) = this%blkLimitsGC(:,:)

    end function enclosingBlock

    ! DEV: For now, the localFlag parameter might be useful for code that
    !      is ported from FLASH4.  We may want to to get rid of it eventually.
    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use physicaldata,    ONLY : unk, &
                                    facevarx, facevary, facevarz
        use gr_specificData, ONLY : scratch, &
                                    scratch_ctr, &
                                    scratch_facevarx, &
                                    scratch_facevary, & 
                                    scratch_facevarz
#ifdef USE_LEVELWIDE_FLUXES
        use gr_specificData, ONLY : gr_flxx, gr_flxy, gr_flxz
#endif

       class(Grid_tile_t), intent(IN)           :: this
#ifdef DEBUG_GRID
       real,               pointer              :: dataPtr(:, :, :, :)
#else
       real,               POINTER_INTENT_OUT   :: dataPtr(:, :, :, :)
#endif
       integer,            intent(IN)           :: gridDataStruct
       logical,            intent(IN), optional :: localFlag

       integer :: lo(1:MDIM), hi(1:MDIM)

#ifdef DEBUG_GRID
       validGridDataStruct = .false.
       validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_FACEZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
       if(.NOT. validGridDataStruct) then
          print *, "Grid_getBlkPtr: gridDataStruct set to improper value"
          print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ," // &
               "WORK or SCRATCH (defined in constants.h)"
          call Driver_abort("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,SCRATCH (see constants.h)")
       end if

       if ((this%id < 1) .OR. (this%id > MAXBLOCKS)) then
          print *, 'Grid_getBlkPtr:  invalid blockid ',this%id
          call Driver_abort("[getDataPtr] invalid blockid ")
       end if

       ! Avoid possible memory leaks
       if (associated(dataPtr)) then
           call Driver_abort("[getDataPtr] Given data pointer must be NULL")
       end if
#endif

       lo = this%blkLimitsGC(LOW, :)

       ! These multifabs are hardwired at creation so that the FAB data only
       ! exists for the block interiors
       if (     (gridDataStruct == FLUXX) & 
           .OR. (gridDataStruct == FLUXY) &
           .OR. (gridDataStruct == FLUXZ)) then
          lo(1:NDIM) = lo(1:NDIM) + NGUARD
          hi         = this%blkLimitsGC(HIGH, :)
          hi(1:NDIM) = hi(1:NDIM) - NGUARD
          if (present(localFlag)) then
              if (localFlag) then
                  hi(1:NDIM) = NGUARD + 1 + hi(1:NDIM) - lo(1:NDIM)
                  lo(1:NDIM) = NGUARD + 1
              end if
          end if
       else if (present(localFlag)) then
           if (localFlag) then
               lo(:) = 1
           end if
       end if

       select case (gridDataStruct)
       case(CENTER)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => unk(:,:,:,:,this%id)
       case(FACEX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarx(:,:,:,:,this%id)
       case(FACEY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevary(:,:,:,:,this%id)
       case(FACEZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => facevarz(:,:,:,:,this%id)
#ifdef USE_LEVELWIDE_FLUXES
       case(FLUXX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxx(:,:,:,:,this%id)
       case(FLUXY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxy(:,:,:,:,this%id)
       case(FLUXZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => gr_flxz(:,:,:,:,this%id)
#else
       case(FLUXX)
          allocate(dataPtr(NFLUXES,lo(1):hi(1)+1, lo(2):hi(2)  , lo(3):hi(3)  ))
       case(FLUXY)
          allocate(dataPtr(NFLUXES,lo(1):hi(1)  , lo(2):hi(2)+1, lo(3):hi(3)  ))
       case(FLUXZ)
          allocate(dataPtr(NFLUXES,lo(1):hi(1)  , lo(2):hi(2)  , lo(3):hi(3)+1))
#endif
       case(SCRATCH)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch(:,:,:,:,this%id)
       case(SCRATCH_CTR)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_ctr(:,:,:,:,this%id)           
       case(SCRATCH_FACEX)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarx(:,:,:,:,this%id)           
       case(SCRATCH_FACEY)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevary(:,:,:,:,this%id)           
       case(SCRATCH_FACEZ)
          dataPtr(1:, lo(1):, lo(2):, lo(3):) => scratch_facevarz(:,:,:,:,this%id)           
       case DEFAULT
          print *, 'TRIED TO GET SOMETHING OTHER THAN UNK OR SCRATCH OR FACE[XYZ]. NOT YET.'
       case(WORK)
          call Driver_abort("work array cannot be got as pointer")
       end select
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(Grid_tile_t), intent(IN)            :: this
        real,                            pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct

        logical :: myAnon

        myAnon = .FALSE.

#ifndef USE_LEVELWIDE_FLUXES
        select case (gridDataStruct)
        case(FLUXX,FLUXY,FLUXZ)
           myAnon = .TRUE.
        end select
#endif

        if (myAnon) then
           if (associated(dataPtr)) deallocate(dataPtr)
        end if

        nullify(dataPtr)
    end subroutine releaseDataPtr

    subroutine fillTileCInfo(this, cInfo)
        use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t
        use Grid_data, ONLY: gr_useOrchestration
        use Grid_data, ONLY: gr_delta
        use,intrinsic :: iso_c_binding
        class(Grid_tile_t), intent(IN)                :: this
        type(Orchestration_tileCInfo_t),intent(OUT)   :: cInfo
        real,pointer,contiguous :: fBlkPtr(:,:,:,:)

#ifdef FLASHX_ORCHESTRATION
#ifdef FLASHX_ORCHESTRATION_MILHOJA
#include "Milhoja.h"
#ifndef RUNTIME_MUST_USE_TILEITER
        if (gr_useOrchestration) then
           cInfo % CInts % nCcComp      = NUNK_VARS
           cInfo % CInts % nFluxComp    = NFLUXES
           cInfo % CInts % loGC(1:MDIM) = this % blklimitsGC(LOW,:)
           cInfo % CInts % hiGC(1:MDIM) = this % blkLimitsGC(HIGH,:)
           cInfo % CInts % lo(1:MDIM)   = this % limits(LOW,:)
           cInfo % CInts % hi(1:MDIM)   = this % limits(HIGH,:)
           cInfo % CInts % ndim         = NDIM
           cInfo % CInts % level        = this % level
           cInfo % CInts % gridIdxOrBlkId = this % id
           cInfo % CInts % tileIdx      = 0 !not meaningful / not know to Grid_tile
           
           cInfo % CReals % deltas(1:MDIM) = gr_delta(1:MDIM, this%level)
           
           cInfo % CPtrs % ccBlkPtr = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(IAXIS) = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(JAXIS) = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(KAXIS) = C_NULL_PTR

           call this % getDataPtr(fBlkPtr, CENTER)
           if(associated(fBlkPtr)) cInfo % CPtrs % ccBlkPtr = c_loc(fBlkPtr)
#ifdef USE_LEVELWIDE_FLUXES
           call this % getDataPtr(fBlkPtr, FLUXX)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(IAXIS) = c_loc(fBlkPtr)
           call this % getDataPtr(fBlkPtr, FLUXY)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(JAXIS) = c_loc(fBlkPtr)
           call this % getDataPtr(fBlkPtr, FLUXZ)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(KAXIS) = c_loc(fBlkPtr)
#endif
        end if
#endif
#endif
#endif
        
      end subroutine fillTileCInfo
end module Grid_tile

