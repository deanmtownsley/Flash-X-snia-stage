#include "Simulation.h"
#include "constants.h"

!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!! Refer to the test's design doc for more information.
!!
!! @todo Once the BC unit has been implemented for Milhoja, uncomment the
!! desired BCs in test_milhoja_grid.par and update expected values here.
!! @todo Should we leave this unittest with simple data that does not refine so
!! that testing the block structure is easy?  All blocks presently on coarsest
!! level since error estimate callback still does nothing.
!! @todo Add in check of BCs on individual tile faces (i.e., faceBCs).
!! @todo Add in check of physical size of each tile.
!! @todo Add in check of nrefs?
!!
!! @param fileUnit   Ignored.  All output is written to stdout.
!! @param perfect    True if no errors occurred; False, otherwise.
subroutine Grid_unitTest(fileUnit, perfect)
    use Grid_interface, ONLY : Grid_getGeometry, &
                               Grid_getDomainBoundBox, &
                               Grid_getDomainBC, &
                               Grid_getDeltas, &
                               Grid_getMaxRefinement, &
                               Grid_getTileIterator, &
                               Grid_releaseTileIterator, &
                               Grid_fillGuardCells
    use Grid_iterator,  ONLY : Grid_iterator_t
    use Grid_tile,      ONLY : Grid_tile_t
    use Grid_data,      ONLY : gr_numRefineVarsMax, &
                               gr_numRefineVars, &
                               gr_refine_var, &
                               gr_refine_cutoff, &
                               gr_derefine_cutoff, &
                               gr_refine_filter, &
                               gr_enforceMaxRefinement, &
                               gr_eosMode, &
                               gr_eosModeInit
    use ut_testDriverMod

    implicit none

    integer, intent(in)    :: fileUnit
    logical, intent(inout) :: perfect

    !!!!! EXPECTED RESULTS BASED ON test_milhoja_grid.par
    !            AND SETUP VALUES GIVEN ABOVE
    integer,  parameter :: NXCELL_EX   = 64
    integer,  parameter :: NYCELL_EX   = 64
    integer,  parameter :: NZCELL_EX   =  4
    integer,  parameter :: NXBLK_EX    =  8
    integer,  parameter :: NYBLK_EX    = 16
    integer,  parameter :: NZBLK_EX    =  2
    real,     parameter :: XMIN_EX     = -1.00
    real,     parameter :: XMAX_EX     =  2.00
    real,     parameter :: YMIN_EX     = -1.50
    real,     parameter :: YMAX_EX     =  4.50
    real,     parameter :: ZMIN_EX     =  0.50
    real,     parameter :: ZMAX_EX     =  0.75
    real,     parameter :: XDELTA_EX   = (XMAX_EX-XMIN_EX)/NXCELL_EX
    real,     parameter :: YDELTA_EX   = (YMAX_EX-YMIN_EX)/NYCELL_EX
    real,     parameter :: ZDELTA_EX   = (ZMAX_EX-ZMIN_EX)/NZCELL_EX
    integer,  parameter :: MAXLEVEL_EX =  4
    integer,  parameter :: XL_BC_EX    = PERIODIC
    integer,  parameter :: XH_BC_EX    = PERIODIC
    integer,  parameter :: YL_BC_EX    = PERIODIC
    integer,  parameter :: YH_BC_EX    = PERIODIC
    integer,  parameter :: ZL_BC_EX    = PERIODIC
    integer,  parameter :: ZH_BC_EX    = PERIODIC

    integer :: geometry
    real    :: domain(LOW:HIGH, 1:MDIM)
    integer :: domainBC(LOW:HIGH, 1:MDIM)
    real    :: deltas(1:MDIM)
    real    :: x_expected
    real    :: y_expected
    real    :: z_expected
    integer :: maxLevel

    type(Grid_iterator_t)         :: itor
    type(Grid_tile_t)             :: tileDesc
    real,                 pointer :: solnData(:, :, :, :)
    integer                       :: n_blocks
    integer                       :: blkLimits(LOW:HIGH, 1:MDIM)
    integer                       :: blkLimitsGC(LOW:HIGH, 1:MDIM)
    integer                       :: blkGC(LOW:HIGH, 1:MDIM)
    integer                       :: blkSize(1:MDIM)
    integer                       :: xBlkMin
    integer                       :: xBlkMax
    integer                       :: yBlkMin
    integer                       :: yBlkMax
    integer                       :: zBlkMin
    integer                       :: zBlkMax
    real                          :: xMin
    real                          :: xMax
    real                          :: yMin
    real                          :: yMax
    real                          :: zMin
    real                          :: zMax
    real                          :: boundBox(LOW:HIGH, 1:MDIM)

    integer :: level
    integer :: i, j, k, var

    nullify(solnData)

    CALL start_test_run

    !!!!! CONFIRM PROPER PHYSICAL DOMAIN
    CALL Grid_getGeometry(geometry)
    CALL assertEqual(geometry, CARTESIAN, "Incorrect coordinate system type")

    CALL Grid_getDomainBoundBox(domain)
    CALL assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    CALL assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
#if NDIM == 1
    CALL assertEqual(domain(LOW,  JAXIS), 0.0,    "Incorrect low Y-coordinate")
    CALL assertEqual(domain(HIGH, JAXIS), 0.0,    "Incorrect high Y-coordinate")
    CALL assertEqual(domain(LOW,  KAXIS), 0.0,    "Incorrect low Z-coordinate")
    CALL assertEqual(domain(HIGH, KAXIS), 0.0,    "Incorrect high Z-coordinate")
#elif NDIM == 2
    CALL assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    CALL assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    CALL assertEqual(domain(LOW,  KAXIS), 0.0,    "Incorrect low Z-coordinate")
    CALL assertEqual(domain(HIGH, KAXIS), 0.0,    "Incorrect high Z-coordinate")
#elif NDIM == 3 
    CALL assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    CALL assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    CALL assertEqual(domain(LOW,  KAXIS), ZMIN_EX,"Incorrect low Z-coordinate")
    CALL assertEqual(domain(HIGH, KAXIS), ZMAX_EX,"Incorrect high Z-coordinate")
#endif

    !!!!! CONFIRM PROPER BC
    CALL Grid_getDomainBC(domainBC)
    CALL assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    CALL assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    CALL assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    CALL assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")
    CALL assertEqual(domainBC(LOW,  KAXIS), ZL_BC_EX, "Incorrect Z-left BC")
    CALL assertEqual(domainBC(HIGH, KAXIS), ZH_BC_EX, "Incorrect Z-right BC")

    !!!!! CONFIRM PROPER REFINEMENT
    CALL Grid_getMaxRefinement(maxLevel, mode=4)
    CALL assertEqual(maxLevel, 1, "Incorrect current finest level")
    CALL Grid_getMaxRefinement(maxLevel, mode=1)
    CALL assertEqual(maxLevel, MAXLEVEL_EX, "Incorrect maximum refine level")

    do level = 1, maxLevel
        CALL Grid_getDeltas(level, deltas)

        x_expected = XDELTA_EX / 2.0**(level - 1)
        y_expected = 0.0
        z_expected = 0.0
#if NDIM >= 2
        y_expected = YDELTA_EX / 2.0**(level - 1)
#endif
#if NDIM == 3 
        z_expected = ZDELTA_EX / 2.0**(level - 1)
#endif
        CALL assertEqual(deltas(IAXIS),x_expected,"Incorrect X mesh spacing")
        CALL assertEqual(deltas(JAXIS),y_expected,"Incorrect Y mesh spacing")
        CALL assertEqual(deltas(KAXIS),z_expected,"Incorrect Z mesh spacing")
    end do

    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0

    ! The tests using this iterator were designed specifically such 
    ! that tiling cannot be used.
    CALL Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

    CALL itor%currentTile(tileDesc)
    xBlkMin = tileDesc%limits(LOW,  IAXIS)
    xBlkMax = tileDesc%limits(HIGH, IAXIS)
    yBlkMin = tileDesc%limits(LOW,  JAXIS)
    yBlkMax = tileDesc%limits(HIGH, JAXIS)
    zBlkMin = tileDesc%limits(LOW,  KAXIS)
    zBlkMax = tileDesc%limits(HIGH, KAXIS)
    xMin =  HUGE(1.0)
    xMax = -HUGE(1.0)
    yMin =  HUGE(1.0)
    yMax = -HUGE(1.0)
    zMin =  HUGE(1.0)
    zMax = -HUGE(1.0)
    do while (itor%isValid())
        n_blocks = n_blocks + 1
        CALL itor%currentTile(tileDesc)

        CALL assertEqual(tileDesc%level, 1, "Incorrect block level")

        CALL tileDesc%boundBox(boundBox)
        xMin = MIN(xMin, boundBox(LOW,  IAXIS))
        xMax = MAX(xMax, boundBox(HIGH, IAXIS))
        yMin = MIN(yMin, boundBox(LOW,  JAXIS))
        yMax = MAX(yMax, boundBox(HIGH, JAXIS))
        zMin = MIN(zMin, boundBox(LOW,  KAXIS))
        zMax = MAX(zMax, boundBox(HIGH, KAXIS))

        ! Check guard cells along all directions
        blkLimits   = tileDesc%limits
        blkLimitsGC = tileDesc%blkLimitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
        CALL assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        CALL assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
#if NDIM == 1
        CALL assertEqual(blkGC(LOW,  JAXIS), 0, "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(HIGH, JAXIS), 0, "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        CALL assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")
#elif NDIM == 2
        CALL assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        CALL assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")
#elif NDIM == 3
        CALL assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        CALL assertEqual(blkGC(LOW,  KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
        CALL assertEqual(blkGC(HIGH, KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
#endif

        ! Correct cells per block along each direction
        blkSize = blkLimits(HIGH, :) - blkLimits(LOW, :) + 1
        CALL assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
#if NDIM == 1
        CALL assertEqual(blkSize(JAXIS), 1, "Incorrect cells per block along Y-axis")
        CALL assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 2
        CALL assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        CALL assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 3
        CALL assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        CALL assertEqual(blkSize(KAXIS), NZCELL_EX / NZBLK_EX, &
                         "Incorrect cells per block along Z-axis")
#endif

        xBlkMin = MIN(xBlkMin, blkLimits(LOW,  IAXIS))
        yBlkMin = MIN(yBlkMin, blkLimits(LOW,  JAXIS))
        zBlkMin = MIN(zBlkMin, blkLimits(LOW,  KAXIS))
        xBlkMax = MAX(xBlkMax, blkLimits(HIGH, IAXIS))
        yBlkMax = MAX(yBlkMax, blkLimits(HIGH, JAXIS))
        zBlkMax = MAX(zBlkMax, blkLimits(HIGH, KAXIS))

        CALL itor%next()
    end do

    CALL Grid_releaseTileIterator(itor)

    ! Confirm proper number of blocks and cells
    CALL assertEqual(xBlkMin, 1, "Incorrect origin X-coordinate")
    CALL assertEqual(yBlkMin, 1, "Incorrect origin Y-coordinate")
    CALL assertEqual(zBlkMin, 1, "Incorrect origin Z-coordinate")

    ! FIXME: This only works for one processor.  Need a reduction here.
#if NDIM == 1
    CALL assertEqual(n_blocks, NXBLK_EX, &
                     "Incorrect total number of blocks")

    CALL assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    CALL assertEqual(yBlkMax, 1, "More than one cell along Y-axis")
    CALL assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    CALL assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    CALL assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    CALL assertEqual(yMin, 1.0,     "Incorrect minimum Y-coordinate found")
    CALL assertEqual(yMax, 1.0,     "Incorrect maximum Y-coordinate found")
    CALL assertEqual(zMin, 1.0,     "Incorrect minimum Z-coordinate found")
    CALL assertEqual(zMax, 1.0,     "Incorrect maximum Z-coordinate found")
#elif NDIM == 2
    CALL assertEqual(n_blocks, NXBLK_EX*NYBLK_EX, &
                     "Incorrect total number of blocks")

    CALL assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    CALL assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    CALL assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    CALL assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    CALL assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    CALL assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    CALL assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    CALL assertEqual(zMin, 1.0,     "Incorrect minimum Z-coordinate found")
    CALL assertEqual(zMax, 1.0,     "Incorrect maximum Z-coordinate found")
#elif NDIM == 3
    CALL assertEqual(n_blocks, NXBLK_EX*NYBLK_EX*NZBLK_EX, &
                     "Incorrect total number of blocks")
    
    CALL assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    CALL assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    CALL assertEqual(zBlkMax, NZCELL_EX, &
                     "Incorrect total number of cells along Z-axis")
    
    CALL assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    CALL assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    CALL assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    CALL assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    CALL assertEqual(zMin, ZMIN_EX, "Incorrect minimum Z-coordinate found")
    CALL assertEqual(zMax, ZMAX_EX, "Incorrect maximum Z-coordinate found")
#endif

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    ! Simulation_initBlock is expected to set data in interior and in GCs
    CALL Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        associate(loGC => tileDesc%grownLimits(LOW, :), &
                  hiGC => tileDesc%grownLimits(HIGH, :))
            CALL tileDesc%getDataPtr(solnData, CENTER)
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = loGC(KAXIS), hiGC(KAXIS)
                    do     j = loGC(JAXIS), hiGC(JAXIS)
                        do i = loGC(IAXIS), hiGC(IAXIS)
                            CALL assertEqual(solnData(i, j, k, var), 1.1 * var, &
                                             "Incorrect initial condition in unk")
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(solnData, CENTER)
        end associate

        CALL itor%next()
    end do
    CALL Grid_releaseTileIterator(itor)

    !!!!! CONFIRM WRITING OF DATA
    CALL Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        associate(lo   => tileDesc%limits(LOW, :), &
                  hi   => tileDesc%limits(HIGH, :), &
                  loGC => tileDesc%grownLimits(LOW, :), &
                  hiGC => tileDesc%grownLimits(HIGH, :))
            CALL tileDesc%getDataPtr(solnData, CENTER)

            ! Set obviously different data in GCs so that we can
            ! confirm that GC fill worked.
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = loGC(KAXIS), hiGC(KAXIS)
                    do     j = loGC(JAXIS), hiGC(JAXIS)
                        do i = loGC(IAXIS), hiGC(IAXIS)
                            solnData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do

            ! Set test data only in interiors using values different
            ! from the ICs
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            solnData(i, j, k, var) = -2.2 * var
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(solnData, CENTER)
        end associate

        CALL itor%next()
    end do
    CALL Grid_releaseTileIterator(itor)

    CALL Grid_fillGuardCells(CENTER, ALLDIR)

    ! Confirm that interior & GC values set correctly
    CALL Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        CALL itor%currentTile(tileDesc)

        associate(loGC => tileDesc%grownLimits(LOW, :), &
                  hiGC => tileDesc%grownLimits(HIGH, :))
            CALL tileDesc%getDataPtr(solnData, CENTER)
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = loGC(KAXIS), hiGC(KAXIS)
                    do     j = loGC(JAXIS), hiGC(JAXIS)
                        do i = loGC(IAXIS), hiGC(IAXIS)
                            CALL assertEqual(solnData(i, j, k, var), &
                                             -2.2 * var, &
                                             "Incorrect data in unk")
                        end do
                    end do
                end do
            end do
            CALL tileDesc%releaseDataPtr(solnData, CENTER)
        end associate

        CALL itor%next()
    end do
    CALL Grid_releaseTileIterator(itor)

    !!!!! CONFIRM REFINEMENT SETUP
    ! uses default value
    call assertEqual(gr_numRefineVarsMax, 4, "Incorrect max refinement variables")
    call assertEqual(gr_numRefineVars, 3, "Incorrect max refinement variables")

    call assertEqual(gr_refine_var(1), DENS_VAR, "First refine var not DENS")
    call assertEqual(gr_refine_var(2), TEMP_VAR, "Second refine var not TEMP")
    call assertEqual(gr_refine_var(3), ENER_VAR, "Third refine var not ENER")

    call assertEqual(gr_refine_cutoff(1), 0.8, "Incorrect DENS refine cutoff")
    call assertEqual(gr_refine_cutoff(2), 0.5, "Incorrect TEMP refine cutoff")
    call assertEqual(gr_refine_cutoff(3), 0.6, "Incorrect ENER refine cutoff")

    call assertEqual(gr_derefine_cutoff(1), 0.45, &
                     "Incorrect DENS derefine cutoff")
    call assertEqual(gr_derefine_cutoff(2), 0.325, &
                     "Incorrect TEMP derefine cutoff")
    call assertEqual(gr_derefine_cutoff(3), 0.35, &
                     "Incorrect ENER derefine cutoff")

    call assertEqual(gr_refine_filter(1), 0.05, &
                     "Incorrect DENS derefine cutoff")
    call assertEqual(gr_refine_filter(2), 0.025, &
                     "Incorrect TEMP derefine cutoff")
    call assertEqual(gr_refine_filter(3), 0.035, &
                     "Incorrect ENER derefine cutoff")

    call assertFalse(gr_enforceMaxRefinement, "gr_enforceMaxRefinement True")

    !!!!! CONFIRM EoS SETUP
    call assertEqual(gr_eosMode, MODE_DENS_EI, &
                     "Incorrect eosMode")
    call assertEqual(gr_eosModeInit, MODE_DENS_TEMP, &
                     "Incorrect eosModeInit")

    perfect = finish_test_run()
end subroutine Grid_unitTest

