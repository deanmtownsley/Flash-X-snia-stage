#include "Simulation.h"
#include "constants.h"

!> A subset of simulation configuration data is loaded into AMReX at
!! initialization and is therefore owned by AMReX.  As a result, AMReX is used
!! to provide these data values to client code through the associated public
!! Grid_* and local gr_* interface accessor routines as well as local data
!! members.
!!
!! This code tests that AMReX is properly initialized for a Cartesian domain by
!! verifying correct results as obtained through the accessor routines and data
!! members.
!!
!! This test is not presently MPI-compatible in the sense that its results can
!! only be trusted if run with a single processor.
!!
!! @param fileUnit   Ignored.  All output is written to stdout.
!! @param perfect    True if no errors occurred; False, otherwise.
subroutine Grid_unitTest(fileUnit, perfect)
    use amrex_fort_module,     ONLY : amrex_spacedim
    use amrex_parallel_module, ONLY : amrex_parallel_myproc

    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getCellCoords, &
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getTileIterator, &
                                      Grid_releaseTileIterator
    use Grid_data,             ONLY : gr_meshMe, &
                                      gr_numRefineVarsMax, gr_numRefineVars, &
                                      gr_refine_var, &
                                      gr_refine_cutoff, gr_derefine_cutoff, &
                                      gr_refine_filter, &
                                      gr_enforceMaxRefinement, &
                                      gr_eosMode, &
                                      gr_eosModeInit
    use Grid_iterator,         ONLY : Grid_iterator_t
    use Grid_tile,             ONLY : Grid_tile_t
    use ut_testDriverMod

    implicit none

    integer, intent(in)    :: fileUnit
    logical, intent(inout) :: perfect

    !!!!! EXPECTED RESULTS BASED ON flash.par AND SETUP VALUES GIVEN ABOVE
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
    integer,  parameter :: XL_BC_EX    = OUTFLOW
    integer,  parameter :: XH_BC_EX    = REFLECTING
    integer,  parameter :: YL_BC_EX    = DIRICHLET
    integer,  parameter :: YH_BC_EX    = DIODE
    integer,  parameter :: ZL_BC_EX    = PERIODIC
    integer,  parameter :: ZH_BC_EX    = DIRICHLET

    integer :: geometry
    real    :: domain(LOW:HIGH, MDIM)
    integer :: domainBC(LOW:HIGH, MDIM)
    real    :: deltas(1:MDIM)
    real    :: x_expected
    real    :: y_expected
    real    :: z_expected
    integer :: max_level

    type(Grid_iterator_t) :: itor
    type(Grid_tile_t)     :: tileDesc
    real, pointer         :: solnData(:, :, :, :)
    integer               :: n_blocks
    integer               :: blkLimits(LOW:HIGH, 1:MDIM)
    integer               :: blkLimitsGC(LOW:HIGH, 1:MDIM)
    integer               :: blkGC(LOW:HIGH, 1:MDIM)
    integer               :: blkSize(1:MDIM)
    integer               :: xBlkMin
    integer               :: xBlkMax
    integer               :: yBlkMin
    integer               :: yBlkMax
    integer               :: zBlkMin
    integer               :: zBlkMax
    real                  :: xMin
    real                  :: xMax
    real                  :: yMin
    real                  :: yMax
    real                  :: zMin
    real                  :: zMax
    real                  :: boundBox(LOW:HIGH, 1:MDIM)

    real, allocatable :: x_coords(:)
    real, allocatable :: y_coords(:)
    real, allocatable :: z_coords(:)

    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)
    integer :: offset(1:MDIM)

    integer :: rank
    integer :: ilev
    integer :: i, j, k, var

    nullify(solnData)

    rank = amrex_parallel_myproc()

    !!!!! CONFIRM PROPER COORDINATE SYSTEM
    ! Dimensionality
    write(*,*)
    if (amrex_spacedim /= NDIM) then
        write(*,*) "Wrong dimensionality - ", amrex_spacedim, ' != ', NDIM
        write(*,*) "Recompile AMReX with correct dimensionality"
        write(*,*)
        stop
    end if

    call start_test_run

    !!!!! CONFIRM MPI SETUP
    call assertEqual(rank, gr_meshMe, "AMReX/FLASH ranks are different")
 
    ! Physical domain
    call Grid_getGeometry(geometry)
    call assertEqual(geometry, CARTESIAN, "Incorrect coordinate system type")

    call Grid_getDomainBoundBox(domain)
#if NDIM == 1
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), 0.0,    "Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), 0.0,    "Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0,    "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0,    "Incorrect high Z-coordinate")
#elif NDIM == 2
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), 0.0,    "Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), 0.0,    "Incorrect high Z-coordinate")
#elif NDIM == 3 
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
    call assertEqual(domain(LOW,  KAXIS), ZMIN_EX,"Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), ZMAX_EX,"Incorrect high Z-coordinate")
#endif

    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0**(ilev - 1)
        y_expected = YDELTA_EX / 2.0**(ilev - 1)
        z_expected = ZDELTA_EX / 2.0**(ilev - 1)
#if NDIM == 1
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),0.0,     "Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0,     "Incorrect high Z-coordinate")
#elif NDIM == 2
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),0.0,       "Incorrect high Z-coordinate")
#elif NDIM == 3 
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect high X-coordinate")
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect high Y-coordinate")
        call assertEqual(deltas(KAXIS),z_expected,"Incorrect high Z-coordinate")
#endif
    end do

    ! DEV: TODO Once unittest is refining mesh, check Grid_getMaxRefinement
    ! with mode that checks actual number of levels in use
    
    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0

    ! The tests using this iterator were designed specifically such 
    ! that tiling cannot be used.
    call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)

    call itor%currentTile(tileDesc)
    xBlkMin = tileDesc%limits(LOW,  IAXIS)
    xBlkMax = tileDesc%limits(HIGH, IAXIS)
    yBlkMin = tileDesc%limits(LOW,  JAXIS)
    yBlkMax = tileDesc%limits(HIGH, JAXIS)
    zBlkMin = tileDesc%limits(LOW,  KAXIS)
    zBlkMax = tileDesc%limits(HIGH, KAXIS)
    ! DEV: TODO Do better than this
    xMin = 1.0e10
    xMax = -xMin
    yMin = 1.0e10
    yMax = -yMin
    zMin = 1.0e10
    zMax = -zMin
    do while (itor%isValid())
        n_blocks = n_blocks + 1
        call itor%currentTile(tileDesc)

        call tileDesc%boundBox(boundBox)
        xMin = MIN(xMin, boundBox(LOW,  IAXIS))
        xMax = MAX(xMax, boundBox(HIGH, IAXIS))
        yMin = MIN(yMin, boundBox(LOW,  JAXIS))
        yMax = MAX(yMax, boundBox(HIGH, JAXIS))
        zMin = MIN(zMin, boundBox(LOW,  KAXIS))
        zMax = MAX(zMax, boundBox(HIGH, KAXIS))

        ! DEVNOTE: Should we leave this unittest with simple data
        ! that does not refine so that testing the block structure is easy?
        ! All blocks on coarsest level since no refining
!        call assertEqual(block%level, 1, "Incorrect block level")

        ! Check guard cells along all directions
        blkLimits   = tileDesc%limits
        blkLimitsGC = tileDesc%blkLimitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
#if NDIM == 1
        call assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), 0, "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), 0, "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")
#elif NDIM == 2
        call assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), 0, "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), 0, "Incorrect guard cell along Z-axis")
#elif NDIM == 3
        call assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(LOW,  KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
#endif

        ! Correct cells per block along each direction
        blkSize = blkLimits(HIGH, :) - blkLimits(LOW, :) + 1
#if NDIM == 1
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), 1, "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 2
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), 1, "Incorrect cells per block along Z-axis")
#elif NDIM == 3
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
        call assertEqual(blkSize(KAXIS), NZCELL_EX / NZBLK_EX, &
                         "Incorrect cells per block along Z-axis")
#endif

        xBlkMin = MIN(xBlkMin, blkLimits(LOW,  IAXIS))
        yBlkMin = MIN(yBlkMin, blkLimits(LOW,  JAXIS))
        zBlkMin = MIN(zBlkMin, blkLimits(LOW,  KAXIS))
        xBlkMax = MAX(xBlkMax, blkLimits(HIGH, IAXIS))
        yBlkMax = MAX(yBlkMax, blkLimits(HIGH, JAXIS))
        zBlkMax = MAX(zBlkMax, blkLimits(HIGH, KAXIS))

        call itor%next()
    end do

    call Grid_releaseTileIterator(itor)
    
    ! Confirm proper number of blocks and cells
    call assertEqual(xBlkMin, 1, "Incorrect origin X-coordinate")
    call assertEqual(yBlkMin, 1, "Incorrect origin Y-coordinate")
    call assertEqual(zBlkMin, 1, "Incorrect origin Z-coordinate")

    ! FIXME: This only works for one processor.  Need a reduction here.
#if NDIM == 1
    call assertEqual(n_blocks, NXBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, 1, "More than one cell along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, 1.0,     "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, 1.0,     "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, 1.0,     "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, 1.0,     "Incorrect maximum Z-coordinate found")
#elif NDIM == 2
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX, &
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax, 1, "More than one cell along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, 1.0,     "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, 1.0,     "Incorrect maximum Z-coordinate found")
#elif NDIM == 3
    call assertEqual(n_blocks, NXBLK_EX*NYBLK_EX*NZBLK_EX, &
                     "Incorrect total number of blocks")
    
    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
    call assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
    call assertEqual(zBlkMax, NZCELL_EX, &
                     "Incorrect total number of cells along Z-axis")
    
    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
    call assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
    call assertEqual(zMin, ZMIN_EX, "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, ZMAX_EX, "Incorrect maximum Z-coordinate found")
#endif

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")
    call assertEqual(domainBC(LOW,  KAXIS), ZL_BC_EX, "Incorrect Z-left BC")
    call assertEqual(domainBC(HIGH, KAXIS), ZH_BC_EX, "Incorrect Z-right BC")

    !!!!! CONFIRM REFINEMENT SETUP
    ! TODO: Get nrefs from AMReX

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        call itor%currentTile(tileDesc)

        associate(lo => tileDesc%limits(LOW, :), &
                  hi => tileDesc%limits(HIGH, :))
        call tileDesc%getDataPtr(solnData, CENTER)
        do           var = UNK_VARS_BEGIN, UNK_VARS_END 
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        call assertEqual(solnData(i, j, k, var), &
                                         1.1 * var, &
                                         "Incorrect initial condition in unk")
                    end do
                end do
            end do
        end do
        call tileDesc%releaseDataPtr(solnData, CENTER)
    !!!!! CONFIRM PROPER INITIAL CONDITIONS for FACEVARS
#if NFACE_VARS>0
        call tileDesc%getDataPtr(solnData, FACEX)
        do           var = 1, NFACE_VARS 
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)+1
                        call assertEqual(solnData(i, j, k, var), &
                                         (1.3*i*i) * var, &
                                         "Incorrect initial condition")
                    end do
                end do
            end do
        end do
        call tileDesc%releaseDataPtr(solnData, FACEX)
#if NDIM>1
        call tileDesc%getDataPtr(solnData, FACEY)
        do           var = 1, NFACE_VARS
            do         k = lo(KAXIS), hi(KAXIS)
                do     j = lo(JAXIS), hi(JAXIS)+1
                    do i = lo(IAXIS), hi(IAXIS)
                        call assertEqual(solnData(i, j, k, var), &
                                         (i*i+1.2*j*j) * var, &
                                         "Incorrect initial condition")
                    end do
                end do
            end do
        end do
        call tileDesc%releaseDataPtr(solnData, FACEY)
#endif
#if NDIM>2
        call tileDesc%getDataPtr(solnData, FACEZ)
        do           var = 1, NFACE_VARS
            do         k = lo(KAXIS), hi(KAXIS)+1
                do     j = lo(JAXIS), hi(JAXIS)
                    do i = lo(IAXIS), hi(IAXIS)
                        call assertEqual(solnData(i, j, k, var), &
                                         (i*i+j*j+1.1*k*k) * var, &
                                         "Incorrect initial condition")
                    end do
                end do
            end do
        end do
        call tileDesc%releaseDataPtr(solnData, FACEZ)
#endif
#endif
        end associate
        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    !!!!! CONFIRM CELL COORDINATE ACCESSORS
    ! TEST THAT COORDINATE FUNCTIONS ARE CORRECT
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       lo = tileDesc%limits(LOW,  :)
       hi = tileDesc%limits(HIGH, :)
       allocate(x_coords(lo(IAXIS):hi(IAXIS)))
       allocate(y_coords(lo(JAXIS):hi(JAXIS)))
       allocate(z_coords(lo(KAXIS):hi(KAXIS)))

       call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, LEFT_EDGE, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                  call assertEqual(x_coords(i), XMIN_EX + (i-1)*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), 0.0,                       "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                       "Bad Z-coordinate")
#elif NDIM == 2
                  call assertEqual(x_coords(i), XMIN_EX + (i-1)*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + (j-1)*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                       "Bad Z-coordinate")
#elif NDIM == 3
                  call assertEqual(x_coords(i), XMIN_EX + (i-1)*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + (j-1)*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), ZMIN_EX + (k-1)*ZDELTA_EX, "Bad Z-coordinate")
#endif
             end do
          end do
       end do

       call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                 call assertEqual(x_coords(i), XMIN_EX + (i-0.5)*XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), 0.0,                         "Bad Y-coordinate")
                 call assertEqual(z_coords(k), 0.0,                         "Bad Z-coordinate")
#elif NDIM == 2
                 call assertEqual(x_coords(i), XMIN_EX + (i-0.5)*XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), YMIN_EX + (j-0.5)*YDELTA_EX, "Bad Y-coordinate")
                 call assertEqual(z_coords(k), 0.0,                         "Bad Z-coordinate")
#elif NDIM == 3
                 call assertEqual(x_coords(i), XMIN_EX + (i-0.5)*XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), YMIN_EX + (j-0.5)*YDELTA_EX, "Bad Y-coordinate")
                 call assertEqual(z_coords(k), ZMIN_EX + (k-0.5)*ZDELTA_EX, "Bad Z-coordinate")
#endif
             end do
          end do
       end do

       call Grid_getCellCoords(IAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, RIGHT_EDGE, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                  call assertEqual(x_coords(i), XMIN_EX + i*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), 0.0,                   "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                   "Bad Z-coordinate")
#elif NDIM == 2
                  call assertEqual(x_coords(i), XMIN_EX + i*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + j*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), 0.0,                   "Bad Z-coordinate")
#elif NDIM == 3
                  call assertEqual(x_coords(i), XMIN_EX + i*XDELTA_EX, "Bad X-coordinate")
                  call assertEqual(y_coords(j), YMIN_EX + j*YDELTA_EX, "Bad Y-coordinate")
                  call assertEqual(z_coords(k), ZMIN_EX + k*ZDELTA_EX, "Bad Z-coordinate")
#endif
             end do
          end do
       end do

       deallocate(x_coords)
       deallocate(y_coords)
       deallocate(z_coords)

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    ! CONFIRM THAT WE CAN GET COORDINATES OF TILE FACES
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       lo = tileDesc%grownLimits(LOW,  :)
       hi = tileDesc%grownLimits(HIGH, :)
       allocate(x_coords(lo(IAXIS):hi(IAXIS)+1))
       allocate(y_coords(lo(JAXIS):hi(JAXIS)  ))
       allocate(z_coords(lo(KAXIS):hi(KAXIS)  ))
       call Grid_getCellCoords(IAXIS, FACES, tileDesc%level, &
                               lo, hi, x_coords)
       call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                               lo, hi, y_coords)
       call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                               lo, hi, z_coords)
       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)+1
#if   NDIM == 1
                 call assertEqual(x_coords(i), XMIN_EX + (i-1)*XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), 0.0,                       "Bad Y-coordinate")
                 call assertEqual(z_coords(k), 0.0,                       "Bad Z-coordinate")
#elif NDIM == 2
                 call assertEqual(x_coords(i), XMIN_EX + (i-1)  *XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), YMIN_EX + (j-0.5)*YDELTA_EX, "Bad Y-coordinate")
                 call assertEqual(z_coords(k), 0.0,                         "Bad Z-coordinate")
#elif NDIM == 3
                 call assertEqual(x_coords(i), XMIN_EX + (i-1)  *XDELTA_EX, "Bad X-coordinate")
                 call assertEqual(y_coords(j), YMIN_EX + (j-0.5)*YDELTA_EX, "Bad Y-coordinate")
                 call assertEqual(z_coords(k), ZMIN_EX + (k-0.5)*ZDELTA_EX, "Bad Z-coordinate")
#endif
             end do
          end do
       end do

       deallocate(x_coords)
       deallocate(y_coords)
       deallocate(z_coords)

       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

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

