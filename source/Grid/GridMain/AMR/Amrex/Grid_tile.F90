!!****ih* source/Grid/GridMain/AMR/Amrex/Grid_tile
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
!!  NOTES
!!    It is generally expected that these objects will never be created directly
!!    in FLASH and therefore there are no "constructors" specified.  In the
!!    even that low-level code in the Grid unit should need to create an object
!!    it can be done by directly setting all public fields to appropriate
!!    values.
!!
!!    DEV: TODO Fix this documentation.
!!    The data contained in the fields of this derived type are
!!
!!    level - the refinement level at which the block resides.  Levels are
!!            referenced with the index set {1, 2, ..., gr_lRefineMax}, where
!!            1 refers to the coarsest level and mesh refinement increases with
!!            increasing level index.  The value of gr_lRefineMax is fixed to
!!            the value of the runtime parameter lrefine_max.
!!
!!    grid_index - used internally by AMReX.  Client code should never make
!!                 direct use of the value stored in this field.
!!
!!    tile_index - used internally by AMReX.  Client code should never make
!!                 direct use of the value stored in this field.
!!
!!    limits - a two-dimensional array where limits(LOW, :) is the array
!!             containing the (i, j, k) coordinate in index space of the 
!!             lower-leftmost corner cell of the block.  Similarly, 
!!             limits(HIGH, :) is the array containing the (i, j, k) 
!!             coordinate of the upper-rightmost corner cell of the block.
!!             The two coordinates are specified for the cell-center index space
!!             defined so that 
!!                - i is an integer in {1, 2, ..., n_x_cells},
!!                - j is an integer in {1, 2, ..., n_y_cells, and
!!                - k is an integer in {1, 2, ..., n_z_cells},
!!             where for FLASH, 
!!                  n_[xyz]_cells = 2^(level - 1) * N[XYZ]B * nblock[xyz].
!!             This index space is the global index space so that (1, 1, 1) is
!!             the lower-leftmost cell in the full domain of the problem.
!!
!!    blkLimitsGC - a two-dimensional array identical to limits except that the
!!                  two cell coordinates in index space refer to the lower-leftmost
!!                  and the upper-rightmost guardcells for the tile's subblock.
!!                  Note that if the tile is in fact a block, then blkLimitsGC 
!!                  has the same contents as limitsGC.  The coordinates are
!!                  specified w.r.t. an expanded version of the global
!!                  cell-center index space used for limits so that
!!                   - i is an integer in {1-NGUARD, ..., n_x_cells+NGUARD},
!!                   - j is an integer in {1-NGUARD, ..., n_y_cells+NGUARD},
!!                   - k is an integer in {1-NGUARD, ..., n_z_cells+NGUARD}.
!!                  Hence, (1, 1, 1) still refers to the lower-leftmost guardcell for
!!                  the full domain of the problem.  If gr_[ijk]guard=4, then 
!!                  (-3, -3, -3) indexes the lower-leftmost guardcell in the
!!                  full domain.
!!
!!    limitsGC - a two-dimensional array identical to limits except that the
!!               two cell coordinates in index space refer to the lower-leftmost
!!               and the upper-rightmost guardcells for the block.  The
!!               coordinates are specified w.r.t. an expanded version of the
!!               global cell-center index space used for limits so that
!!                - i is an integer in {1-NGUARD, ..., n_x_cells+NGUARD},
!!                - j is an integer in {1-NGUARD, ..., n_y_cells+NGUARD},
!!                - k is an integer in {1-NGUARD, ..., n_z_cells+NGUARD}.
!!               Hence, (1, 1, 1) still refers to the lower-leftmost guardcell for
!!               the full domain of the problem.  If gr_[ijk]guard=4, then 
!!               (-3, -3, -3) indexes the lower-leftmost guardcell in the
!!               full domain.
!!
!!****

#include "constants.h"
#include "Simulation.h"

module Grid_tile
    implicit none

    private

    type, public :: Grid_tile_t
        integer, public :: level
        integer, public :: grid_index
        integer, public :: tile_index
        integer, public :: limits(LOW:HIGH, MDIM)
        integer, public :: grownLimits(LOW:HIGH, MDIM)
        integer, public :: blkLimitsGC(LOW:HIGH, MDIM)
    contains
        procedure, public :: deltas
        procedure, public :: boundBox
        procedure, public :: physicalSize
        procedure, public :: faceBCs
!        procedure, public :: outsideDomain
        procedure, public :: getDataPtr
        procedure, public :: releaseDataPtr
        procedure, public :: enclosingBlock
    end type Grid_tile_t

contains

    subroutine deltas(this, dx)
        use amrex_amrcore_module, ONLY : amrex_geom

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: dx(1:MDIM)

        ! AMReX uses zero-based level indexing, but FLASH assumes one-based
        dx(:) = 0.0
        dx(1:MDIM) = amrex_geom(this%level - 1)%dx(1:MDIM)
    end subroutine deltas

    subroutine boundBox(this, box)
        use amrex_amrcore_module,  ONLY : amrex_geom
        use amrex_geometry_module, ONLY : amrex_problo

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: box(LOW:HIGH, 1:MDIM)

        ! DEV: FIXME How to manage matching amrex_real to FLASH real
        box(:, :) = 1.0
        associate(x0 => amrex_problo, &
                  dx => amrex_geom(this%level - 1)%dx, &
                  lo => this%limits(LOW,  :), &
                  hi => this%limits(HIGH, :))
!            ! lo is 1-based cell-index of lower-left cell in block 
!            ! hi is 1-based cell-index of upper-right cell in block
            box(LOW,  1:NDIM) = x0(1:NDIM) + (lo(1:NDIM) - 1)*dx(1:NDIM)
            box(HIGH, 1:NDIM) = x0(1:NDIM) + (hi(1:NDIM)    )*dx(1:NDIM)
        end associate
    end subroutine boundBox 

    ! DEV: FIXME What should the size be for the dimensions > NDIM?
    !            Should we follow the same rules used to compute
    !            cell volumes (See User Manual).
    subroutine physicalSize(this, tileSize) 
        use amrex_amrcore_module, ONLY : amrex_geom

        class(Grid_tile_t), intent(IN)  :: this
        real,               intent(OUT) :: tileSize(1:MDIM) 

        tileSize(:) = 0.0

        associate(dx => amrex_geom(this%level - 1)%dx, &
                  lo => this%limits(LOW,  :), &
                  hi => this%limits(HIGH, :))
            tileSize(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) * dx(1:NDIM)
        end associate
    end subroutine physicalSize
    
    subroutine faceBCs(this, faces, onBoundary)
        use amrex_amrcore_module, ONLY : amrex_geom
        use Grid_data,            ONLY : gr_globalDomain, &
                                         gr_domainBC

        class(Grid_tile_t), intent(IN)            :: this
        integer,            intent(OUT)           :: faces(LOW:HIGH, 1:MDIM)
        integer,            intent(OUT), optional :: onBoundary(LOW:HIGH, 1:MDIM)
        
        real    :: bnd_box(LOW:HIGH, 1:MDIM)
        real    :: deltas(1:MDIM)
        integer :: axis, face
 
        deltas(:) = amrex_geom(this%level - 1)%dx(:)

        call this%boundBox(bnd_box)
        do    axis = 1, MDIM
           do face = LOW, HIGH
              faces(face, axis) = NOT_BOUNDARY
              if (present (onBoundary)) then
                 onBoundary(face,axis) = NOT_BOUNDARY
              end if

              if (almostEqual(bnd_box(face, axis), &
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

    function enclosingBlock(this)
        use amrex_fort_module,    ONLY : wp => amrex_real
        use amrex_box_module,     ONLY : amrex_box

        use gr_physicalMultifabs, ONLY : unk

        class(Grid_tile_t), intent(IN)  :: this
        type(Grid_tile_t)               :: enclosingBlock

        type(amrex_box) :: parent_box
        integer         :: idx(1:MDIM+1)

        real(wp), pointer :: parent_dataPtr(:, :, :, :)

        enclosingBlock%level      = this%level
        enclosingBlock%grid_index = this%grid_index
        enclosingBlock%tile_index = 0

        associate(level => this%level - 1, &
                  gId   => this%grid_index)
            ! The boxes in the boxarray describe the interiors of
            ! our blocks
            parent_box = unk(level)%ba%get_box(gId)
            enclosingBlock%limits(:, :) = 1
            enclosingBlock%limits(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%limits(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1

            call parent_box%grow(NGUARD)
            enclosingBlock%grownLimits(:, :) = 1
            enclosingBlock%grownLimits(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%grownLimits(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1

            enclosingBlock%blkLimitsGC(:, :) = 1
            enclosingBlock%blkLimitsGC(LOW,  1:NDIM) = parent_box%lo(1:NDIM) + 1
            enclosingBlock%blkLimitsGC(HIGH, 1:NDIM) = parent_box%hi(1:NDIM) + 1
        end associate
    end function enclosingBlock

    ! DEV: If the client code requests a pointer to data that is not 
    ! included in the problem, this routine will return a null pointer
    ! without indicating an error.
    !
    ! This gives the client code the possibility to use either preprocessor
    ! checks to avoid calling this routine needlessly or to do runtime checks
    ! of pointers.
    !
    ! DEV: For now, the localFlag parameter might be useful as we pull in more
    !      units from FLASH4.4.  Try to get rid of it along the way.
    subroutine getDataPtr(this, dataPtr, gridDataStruct, localFlag)
        use amrex_fort_module,      ONLY : wp => amrex_real

        use gr_physicalMultifabs,   ONLY : unk, &
                                           gr_scratchCtr, &
                                           facevarx, facevary, facevarz
#ifdef USE_LEVELWIDE_FLUXES
        use gr_physicalMultifabs,   ONLY : fluxes
#endif

        class(Grid_tile_t), intent(IN),  target   :: this
        real(wp),                        pointer  :: dataPtr(:, :, :, :)
        integer,            intent(IN)            :: gridDataStruct
        logical,            intent(IN),  optional :: localFlag

        integer :: lo(1:MDIM), hi(1:MDIM)

#ifdef DEBUG_GRID
       validGridDataStruct = .false.
       validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FLUXX).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FLUXY).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == FLUXZ).or.validGridDataStruct
       validGridDataStruct= (gridDataStruct == SCRATCH_CTR).or.validGridDataStruct
       if(.NOT. validGridDataStruct) then
          print *, "Grid_getBlkPtr: gridDataStruct set to improper value"
          print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ," // &
               " or SCRATCH_CTR (defined in constants.h)"
          call Driver_abort("gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,SCRATCH_CTR (see constants.h)")
       end if

#endif

        ! Avoid possible memory leaks
        if (associated(dataPtr)) then
            call Driver_abort("[getDataPtr] Given data pointer must be NULL")
        end if

        lo = this%blkLimitsGC(LOW, :)

        ! These multifabs are hardwired at creation so that the FAB data only
        ! exists for the block interiors
        if (     (gridDataStruct == SCRATCH_CTR) .OR. (gridDataStruct == FLUXX) & 
            .OR. (gridDataStruct == FLUXY)       .OR. (gridDataStruct == FLUXZ)) then
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

        ! Multifab arrays use 0-based level index set (AMReX) instead of 
        ! 1-based set (FLASH/block)
        associate (ilev => this%level - 1, &
                   igrd => this%grid_index)
          select case (gridDataStruct)
          case(CENTER)
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => unk     (ilev)%dataptr(igrd)
          case(FACEX)
#if NFACE_VARS > 0
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarx(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FACEY)
#if NFACE_VARS > 0 && NDIM >= 2
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevary(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
          case(FACEZ)
#if NFACE_VARS > 0 && NDIM == 3
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => facevarz(ilev)%dataptr(igrd)
#else
             nullify(dataPtr)
#endif
#if NFLUXES > 0 && defined(USE_LEVELWIDE_FLUXES)
          case(FLUXX)
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, IAXIS)%dataptr(igrd)
          case(FLUXY)
#  if NDIM >= 2
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, JAXIS)%dataptr(igrd)
#  endif
          case(FLUXZ)
#  if NDIM == 3
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => fluxes(ilev, KAXIS)%dataptr(igrd)
#  endif
#elif NFLUXES > 0
          case(FLUXX)
             allocate(dataPtr(lo(1):hi(1)+1, lo(2):hi(2)  , lo(3):hi(3)  ,NFLUXES))
          case(FLUXY)
#  if NDIM >= 2
             allocate(dataPtr(lo(1):hi(1)  , lo(2):hi(2)+1, lo(3):hi(3)  ,NFLUXES))
#  else
             allocate(dataPtr(      0      ,       0      ,       0      ,NFLUXES))
#  endif
          case(FLUXZ)
#  if NDIM == 3
             allocate(dataPtr(lo(1):hi(1)  , lo(2):hi(2)  , lo(3):hi(3)+1,NFLUXES))
#  else
             allocate(dataPtr(      0      ,       0      ,       0      ,NFLUXES))
#  endif
#else
          case(FLUXX,FLUXY,FLUXZ); nullify(dataPtr)
#endif
          case(SCRATCH_CTR)
             dataPtr(lo(1):, lo(2):, lo(3):, 1:) => gr_scratchCtr(ilev)%dataptr(igrd)
          case DEFAULT
              call Driver_abort("[getDataPtr] Unknown grid data structure")
          end select
        end associate
    end subroutine getDataPtr

    subroutine releaseDataPtr(this, dataPtr, gridDataStruct)
        use amrex_fort_module, ONLY : wp => amrex_real

        class(Grid_tile_t), intent(IN)            :: this
        real(wp),           intent(OUT), pointer  :: dataPtr(:, :, :, :)
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

end module Grid_tile

