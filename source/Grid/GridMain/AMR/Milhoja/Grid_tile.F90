#include "constants.h"
#include "Simulation.h"

!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! This is a Milhoja-specific implementation of the Grid_tile_t class.  While
!! its public interface must match that of the Grid_tile_t classes defined by other
!! Grid implementations, its implementation can be significantly different.
!! Note that there is no mechanism to enforce this interface requirement.
module Grid_tile
    use iso_c_binding,     ONLY : C_PTR, &
                                  C_NULL_PTR

    implicit none

    private

    public :: Grid_tile_fromMilhojaTilePtr

    !> At a base level each instance of this derived type can be understood to be
    !! a unique Milhoja-specific tile index that Flash-X code passes around
    !! so that code can later access the tile's metadata and data.  Ideally,
    !! application code should have little need to access or interpret the
    !! content of the index itself (e.g., grid_index and tile_index) aside from
    !! potentially level.  The specific composition of the index and how they are
    !! used are likely library-dependent and, therefore, this interface provides an
    !! abstraction layer.  The procedures attached to the index via the class
    !! definition can be interpreted as syntactic sugar for easing data and
    !! metadata access.  For instance, the same information can be accessed by
    !! using an index's deltas procedure and Grid_getDeltas.
    !!
    !! The set {level, grid_index, tile_index} represents the maximal data set
    !! necessary to form a unique tile index across all Milhoja grid backends.
    !! For example, AMReX needs all three, but another backend might only
    !! require grid_index.
    !!
    !! Refer to the Flash-X users' guide for an explanation of the differences
    !! between limits, grownLimits, and blkLimitsGC as well as their uses.
    !!
    !! @todo Can grid_index and tile_index be made private?  If not, we should
    !!       try for this.
    !! @todo Since physics units should no longer use Grid_tile_t, can and
    !! should we simplify this interface?  In other words, is there a need for
    !! syntactic sugar?
    !! @todo Should blkLimitsGC be renamed to something that reflects more
    !! explicitly its role as data array extent/indexing?
    type, public :: Grid_tile_t
        type(C_PTR),          private :: C_dataPtr                     = C_NULL_PTR
        real,        pointer, private :: F_dataPtr(:,:,:,:)            => null()
        real,        pointer, private :: F_fluxXPtr(:,:,:,:)           => null()
        real,        pointer, private :: F_fluxYPtr(:,:,:,:)           => null()
        real,        pointer, private :: F_fluxZPtr(:,:,:,:)           => null()
        integer,              public  :: level                         = -1
        integer,              public  :: grid_index                    = -1
        integer,              public  :: tile_index                    = -1
        integer,              public  :: limits(LOW:HIGH, 1:MDIM)      =  0
        integer,              public  :: grownLimits(LOW:HIGH, 1:MDIM) =  0
        integer,              public  :: blkLimitsGC(LOW:HIGH, 1:MDIM) =  0
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: fillTileCInfo
    end type Grid_tile_t

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS
    ! The C-to-Fortran interoperability layer
    interface
        !> Fortran interface on routine in C interface of same name.
        function milhoja_tile_get_metadata_C(C_dataItemPtr,              &
                                             C_gId, C_level,             &
                                             C_lo, C_hi, C_loGC, C_hiGC, &
                                             C_nCcVars, C_nFluxVars,     &
                                             C_dataPtr,  C_fluxXPtr,     &
                                             C_fluxYPtr, C_fluxZPtr)     &
                                             result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),          intent(IN), value :: C_dataItemPtr
            integer(MILHOJA_INT), intent(OUT)       :: C_gId
            integer(MILHOJA_INT), intent(OUT)       :: C_level
            type(C_PTR),          intent(IN), value :: C_lo
            type(C_PTR),          intent(IN), value :: C_hi
            type(C_PTR),          intent(IN), value :: C_loGC
            type(C_PTR),          intent(IN), value :: C_hiGC
            integer(MILHOJA_INT), intent(OUT)       :: C_nCcVars
            integer(MILHOJA_INT), intent(OUT)       :: C_nFluxVars
            type(C_PTR),          intent(OUT)       :: C_dataPtr
            type(C_PTR),          intent(OUT)       :: C_fluxXPtr
            type(C_PTR),          intent(OUT)       :: C_fluxYPtr
            type(C_PTR),          intent(OUT)       :: C_fluxZPtr
            integer(MILHOJA_INT)                    :: C_ierr
        end function milhoja_tile_get_metadata_C
    end interface 

contains

    !> Create a Grid_tile_t object and populate it with the metadata associated
    !! with the Milhoja C++ Tile object pointed to by the given C pointer.
    !!
    !! Typically, the pointer will have been obtained from a Milhoja tile iterator,
    !! in which case the Tile object would be the current Tile indexed by the
    !! iterator.  Another typical use is in Fortran task functions.
    !!
    !! It is intended that this subroutine only be used in the Fortran/C++
    !! interoperability layer as opposed to general Flash-X application code.
    !! This justifies the presence of C artifacts in the subroutine's interface.
    !!
    !! The implementation is presently only functional with blocks.
    !!
    !! @todo Code up full implementation included tiling.
    !!
    !! @param MH_tile_Cptr   The C pointer associated with the tile of interest
    !! @param tileDesc       The Fortran Tile object
    subroutine Grid_tile_fromMilhojaTilePtr(MH_tile_Cptr, tileDesc)
        use iso_c_binding,       ONLY : C_ASSOCIATED, &
                                        C_F_POINTER,  &
                                        C_LOC

        use milhoja_types_mod,   ONLY : MILHOJA_INT

        use gr_milhojaInterface, ONLY : gr_checkMilhojaError
        use Driver_interface,    ONLY : Driver_abort

        type(C_PTR),       intent(IN)    :: MH_tile_Cptr
        type(Grid_tile_t), intent(INOUT) :: tileDesc

        integer :: dataShape(1:MDIM+1)
        integer :: i

        integer(MILHOJA_INT), target :: MH_lo(1:MDIM)
        integer(MILHOJA_INT), target :: MH_hi(1:MDIM)
        integer(MILHOJA_INT), target :: MH_loGrown(1:MDIM)
        integer(MILHOJA_INT), target :: MH_hiGrown(1:MDIM)
        type(C_PTR)                  :: lo_Cptr
        type(C_PTR)                  :: hi_Cptr
        type(C_PTR)                  :: loGrown_Cptr
        type(C_PTR)                  :: hiGrown_Cptr
        type(C_PTR)                  :: fluxX_Cptr
        type(C_PTR)                  :: fluxY_Cptr
        type(C_PTR)                  :: fluxZ_Cptr
        integer(MILHOJA_INT)         :: MH_gId
        integer(MILHOJA_INT)         :: MH_level
        integer(MILHOJA_INT)         :: MH_nCcVars
        integer(MILHOJA_INT)         :: MH_nFluxVars
        integer(MILHOJA_INT)         :: MH_ierr

        if      (.NOT. C_ASSOCIATED(MH_tile_Cptr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] Null tile pointer")
        else if (C_ASSOCIATED(tileDesc%C_dataPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] C_dataPtr not NULL")
        else if (ASSOCIATED(tileDesc%F_dataPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_dataPtr not NULL")
        else if (ASSOCIATED(tileDesc%F_fluxXPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxXPtr not NULL")
        else if (ASSOCIATED(tileDesc%F_fluxYPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxYPtr not NULL")
        else if (ASSOCIATED(tileDesc%F_fluxZPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxZPtr not NULL")
        end if

        ! Milhoja only sets the 1:NDIM data.  Preset default value for data above NDIM.
        ! These must be set correctly here so that dataShape is correct.
        MH_lo(:)      = 1
        MH_hi(:)      = 1
        MH_loGrown(:) = 1
        MH_hiGrown(:) = 1

        lo_Cptr      = C_LOC(MH_lo)
        hi_Cptr      = C_LOC(MH_hi)
        loGrown_Cptr = C_LOC(MH_loGrown)
        hiGrown_Cptr = C_LOC(MH_hiGrown)

        fluxX_Cptr = C_NULL_PTR
        fluxY_Cptr = C_NULL_PTR
        fluxZ_Cptr = C_NULL_PTR

        MH_ierr = milhoja_tile_get_metadata_C(MH_tile_Cptr,               &
                                              MH_gId, MH_level,           &
                                              lo_Cptr,   hi_Cptr,         &
                                              loGrown_Cptr, hiGrown_Cptr, &
                                              MH_nCcVars, MH_nFluxVars,   &
                                              tileDesc%C_dataPtr,         &
                                              fluxX_Cptr,                 &
                                              fluxY_Cptr,                 &
                                              fluxZ_Cptr)
        CALL gr_checkMilhojaError("Grid_tile_fromMilhojaTilePtr", MH_ierr)
        if (.NOT. C_ASSOCIATED(tileDesc%C_dataPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] C_dataPtr is NULL")
        end if

        dataShape(1:MDIM) = INT(MH_hiGrown(1:MDIM) - MH_loGrown(1:MDIM)) + 1
        dataShape(MDIM+1) = INT(MH_nCcVars)
        CALL C_F_POINTER(tileDesc%C_dataPtr, tileDesc%F_dataPtr, shape=dataShape)
        if (.NOT. ASSOCIATED(tileDesc%F_dataPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_dataPtr is NULL")
        end if

#if NFLUXES > 0
        if (.NOT. C_ASSOCIATED(fluxX_Cptr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] fluxX_Cptr is NULL")
        end if

        dataShape(1:MDIM) = INT(MH_hi(1:MDIM) - MH_lo(1:MDIM)) + 1
        dataShape(MDIM+1) = INT(MH_nFluxVars)
        dataShape(IAXIS) = dataShape(IAXIS) + 1
        CALL C_F_POINTER(fluxX_Cptr, tileDesc%F_fluxXPtr, shape=dataShape)
        if (.NOT. ASSOCIATED(tileDesc%F_fluxXPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxXPtr is NULL")
        end if
#if NDIM >= 2
        if (.NOT. C_ASSOCIATED(fluxY_Cptr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] fluxY_Cptr is NULL")
        end if

        dataShape(1:MDIM) = INT(MH_hi(1:MDIM) - MH_lo(1:MDIM)) + 1
        dataShape(MDIM+1) = INT(MH_nFluxVars)
        dataShape(JAXIS) = dataShape(JAXIS) + 1
        CALL C_F_POINTER(fluxY_Cptr, tileDesc%F_fluxYPtr, shape=dataShape)
        if (.NOT. ASSOCIATED(tileDesc%F_fluxYPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxYPtr is NULL")
        end if
#endif
#if NDIM == 3
        if (.NOT. C_ASSOCIATED(fluxZ_Cptr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] fluxZ_Cptr is NULL")
        end if

        dataShape(1:MDIM) = INT(MH_hi(1:MDIM) - MH_lo(1:MDIM)) + 1
        dataShape(MDIM+1) = INT(MH_nFluxVars)
        dataShape(KAXIS) = dataShape(KAXIS) + 1
        CALL C_F_POINTER(fluxZ_Cptr, tileDesc%F_fluxZPtr, shape=dataShape)
        if (.NOT. ASSOCIATED(tileDesc%F_fluxZPtr)) then
            CALL Driver_abort("[Grid_tile_fromMilhojaTilePtr] F_fluxZPtr is NULL")
        end if
#endif
#endif

        ! Assuming for C interface that level is 1-based.
        tileDesc%level      = INT(MH_level)
        tileDesc%grid_index = INT(MH_gId)
        tileDesc%tile_index = -1

        ! Assuming for C interface that points are defined w.r.t. a 1-based
        ! global index space.
        do i = 1, MDIM
            tileDesc%limits(LOW,  i)      = INT(MH_lo(i))
            tileDesc%limits(HIGH, i)      = INT(MH_hi(i))
            tileDesc%grownLimits(LOW,  i) = INT(MH_loGrown(i))
            tileDesc%grownLimits(HIGH, i) = INT(MH_hiGrown(i))
            ! Assume no proper tiles for now
            tileDesc%blkLimitsGC(LOW,  i) = INT(MH_loGrown(i))
            tileDesc%blkLimitsGC(HIGH, i) = INT(MH_hiGrown(i))
        end do
    end subroutine Grid_tile_fromMilhojaTilePtr

    !> Obtain  the dx/dy/dz grid spacing for the refinement level of the
    !! associated tile.  dx is the size of one cell in the tile in the x
    !! direction.  Similar statements can be made for dy and dz.
    !!
    !! @param dx  The array of grid spacings.  The spacing values are zero for
    !!            those dimensions above NDIM.
    subroutine deltas(this, dx)
        use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                        MILHOJA_REAL
        use milhoja_grid_mod,    ONLY : milhoja_grid_getDeltas

        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        integer(MILHOJA_INT) :: MH_level
        real(MILHOJA_REAL)   :: MH_deltas(1:MDIM)
        integer(MILHOJA_INT) :: MH_ierr

        integer :: i

        MH_level = INT(this%level, kind=MILHOJA_INT)

        CALL milhoja_grid_getDeltas(MH_level, MH_deltas, MH_ierr)
        CALL gr_checkMilhojaError("Grid_tile_t%deltas", MH_ierr)

        dx(:) = 0.0
        do i = 1, NDIM
            dx(i) = REAL(MH_deltas(i))
        end do
    end subroutine deltas

    !> Obtain the bounding box of the interior of the associated tile in
    !! physical space.
    !!
    !! @todo What value to set for dimensions above NDIM?
    !!
    !! @param box  The bounding box specified by its low and high coordinates.
    !!             A value of 1 is returned for both low and high for
    !!             dimensions above NDIM.
    subroutine boundBox(this, box)
        use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                        MILHOJA_REAL
        use milhoja_grid_mod,    ONLY : milhoja_grid_getDeltas

        use Grid_data,           ONLY : gr_globalDomain
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)

        integer(MILHOJA_INT) :: MH_level
        real(MILHOJA_REAL)   :: MH_dx(1:MDIM)
        integer(MILHOJA_INT) :: MH_ierr

        integer :: i

        MH_level = INT(this%level, kind=MILHOJA_INT)

        CALL milhoja_grid_getDeltas(MH_level, MH_dx, MH_ierr)
        CALL gr_checkMilhojaError("Grid_tile_t%boundBox", MH_ierr)

        box(:, :) = 1.0
        associate(x0 => gr_globalDomain(LOW, :), &
                  lo => this%limits(LOW,  :), &
                  hi => this%limits(HIGH, :))
            do i = 1, NDIM
                box(LOW,  i) = x0(i) + (lo(i) - 1)*REAL(MH_dx(i))
                box(HIGH, i) = x0(i) + (hi(i)    )*REAL(MH_dx(i))
            end do
        end associate
    end subroutine boundBox 

    !> Obtain the extent of the associated tile's interior in physical space.
    !!
    !! This subroutine is not implemented and will abort if called.
    !!
    !! @todo Implement
    !! @todo What should the value be above NDIM?  Should it follow the same
    !!       rules as for cell volumes (See User Manual).
    subroutine physicalSize(this, tileSize) 
        use Driver_interface, ONLY : Driver_abort
    
        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM) 

        tileSize(:) = 0.0
        call Driver_abort("[Grid_tile%physicalSize] Not implemented yet")
    end subroutine physicalSize

    !> For each face of the associated tile, obtain the boundary condition type
    !! if the face is on the domain boundary.  Note that BC values are returned
    !! for all faces regardless of problem dimension.
    !!
    !! This subroutine is not implemented and will abort if called.
    !!
    !! @todo Implement
    !!
    !! @param faces  The BCs for the tile's faces.  Faces that are not on the
    !! boundary or that have periodic BCs are set to NOT_BOUNDARY.
    !! @param onBoundary As for faces except that faces with periodic BC have
    !! their BC value set to PERIODIC.
    subroutine faceBCs(this, faces, onBoundary)
        use Driver_interface, ONLY : Driver_abort

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)

        faces(:,:) = NOT_BOUNDARY
        if (present(onBoundary)) then
            onBoundary(:,:) = NOT_BOUNDARY
        end if
        call Driver_abort("[Grid_tile%faceBCs] Not implemented yet")
    end subroutine faceBCs

    !> Obtain the pointer to the array holding the associated tile's data.
    !! For cell-centered data, the pointer gives access to both the interior of
    !! the tile as well as the tile's guardcells should these exist; for flux
    !! data, interior only.
    !!
    !! If this is called, calling code is obligated to subsequently release the
    !! pointer using the releaseDataPtr subroutine.
    !!
    !! This implementation does not presently grant access to all data types nor
    !! implement the full functionality implied by all arguments.  This routine
    !! aborts if calling code attempts to access unavailable data or use 
    !! unimplemented functionality.
    !!
    !! @todo Code full implementation.
    !!
    !! @param dataPtr  The desired pointer.  This must be null on entry.
    !! @param gridDataStruct The type of data to access via the pointer.  Valid
    !!                       values are CENTER and FLUX[XYZ].  It is a logical
    !!                       error to request flux data if NFLUXES is zero.
    !!                       The pointer is null for flux data whose
    !!                       dimension is above NDIM.
    !! @param localFlag  If given and True then the data array's indices are
    !!                   specified in block-local index space; otherwise, the
    !!                   global index space.  
    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use Driver_interface, ONLY : Driver_abort

        class(Grid_tile_t), intent(IN),  target   :: this
        real,                            pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct
        logical,            intent(IN),  optional :: localFlag

        integer :: lo(1:MDIM)
        integer :: loGC(1:MDIM)

        ! Avoid possible memory leaks
        if (ASSOCIATED(dataPtr)) then
            CALL Driver_abort("[Grid_tile_t%getDataPtr] Given data pointer must be NULL")
        end if

        lo(:)   = this%limits(LOW, :)
        loGC(:) = this%grownLimits(LOW, :)
        if (present(localFlag)) then
            if (gridDataStruct /= CENTER) then
                call Driver_abort("[Grid_tile_t%getDataPtr] not implemented yet!")
            else if (localFlag) then
                loGC(:) = 1
            end if
        end if

        if      (gridDataStruct == CENTER) then
            dataPtr(loGC(IAXIS):, loGC(JAXIS):, loGC(KAXIS):, 1:) => this%F_dataPtr
#if NFLUXES > 0
        else if (gridDataStruct == FLUXX) then
            dataPtr(lo(IAXIS):, lo(JAXIS):, lo(KAXIS):, 1:) => this%F_fluxXPtr
        else if (gridDataStruct == FLUXY) then
#if NDIM >= 2
            dataPtr(lo(IAXIS):, lo(JAXIS):, lo(KAXIS):, 1:) => this%F_fluxYPtr
#else
            nullify(dataPtr)
#endif
        else if (gridDataStruct == FLUXZ) then
#if NDIM == 3
            dataPtr(lo(IAXIS):, lo(JAXIS):, lo(KAXIS):, 1:) => this%F_fluxZPtr
#else
            nullify(dataPtr)
#endif
#endif
        else
            CALL Driver_abort("[Grid_tile_t%getDataPtr] Invalid gridDataStruct")
        end if
    end subroutine getDataPtr

    !> Release the associated tile's given data pointer.
    !!
    !! @param dataPtr  The pointer to release.  It is set to null on exit.
    !! @param gridDataStruct The type of data pointed to.
    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        class(Grid_tile_t), intent(IN)          :: this
        real,                          pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)          :: gridDataStruct

        NULLIFY(dataPtr)
    end subroutine releaseDataPtr

    subroutine fillTileCInfo(this, cInfo)
        use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t
        use Grid_data, ONLY: gr_useOrchestration
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError
        use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                        MILHOJA_REAL
        use milhoja_grid_mod, ONLY : milhoja_grid_getDeltas
        use,intrinsic :: iso_c_binding
        class(Grid_tile_t), intent(IN)                :: this
        type(Orchestration_tileCInfo_t),intent(OUT)   :: cInfo

        integer(MILHOJA_INT) :: MH_level
        real(MILHOJA_REAL)   :: MH_deltas(1:MDIM)
        integer(MILHOJA_INT) :: MH_ierr
        integer :: i
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
           cInfo % CInts % gridIdxOrBlkId = this % grid_index
           cInfo % CInts % tileIdx      = this % tile_index ! may not be meaningful

           MH_level = INT(this%level, kind=MILHOJA_INT)
           CALL milhoja_grid_getDeltas(MH_level, MH_deltas, MH_ierr)
           CALL gr_checkMilhojaError("Grid_tile_t%fillTileCInfo", MH_ierr)
           associate(dx => cInfo % CReals % deltas)
               dx(:) = 0.0
               do i = 1, NDIM
                  dx(i) = REAL(MH_deltas(i))
               end do
           end associate

           cInfo % CPtrs % ccBlkPtr = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(IAXIS) = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(JAXIS) = C_NULL_PTR
           cInfo % CPtrs % fluxBlkPtrs(KAXIS) = C_NULL_PTR

           nullify(fBlkPtr)
           call this % getDataPtr(fBlkPtr, CENTER)
           if(associated(fBlkPtr)) cInfo % CPtrs % ccBlkPtr = c_loc(fBlkPtr)
#ifdef USE_LEVELWIDE_FLUXES
           nullify(fBlkPtr)
           call this % getDataPtr(fBlkPtr, FLUXX)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(IAXIS) = c_loc(fBlkPtr)
           nullify(fBlkPtr)
           call this % getDataPtr(fBlkPtr, FLUXY)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(JAXIS) = c_loc(fBlkPtr)
           nullify(fBlkPtr)
           call this % getDataPtr(fBlkPtr, FLUXZ)
           if(associated(fBlkPtr)) cInfo % CPtrs % fluxBlkPtrs(KAXIS) = c_loc(fBlkPtr)
#endif
        end if
#endif
#endif
#endif

    end subroutine fillTileCInfo
end module Grid_tile

