#include "Simulation.h"
#include "constants.h"

!!REORDER(4): solnData
subroutine Grid_unitTest(fileUnit, perfect)
    use Grid_interface,        ONLY : Grid_getDomainBoundBox, &
                                      Grid_getCellCoords, &
                                      Grid_getCellFaceAreas, &
                                      Grid_getGeometry, &
                                      Grid_getDeltas, &
                                      Grid_getMaxRefinement, &
                                      Grid_getTileIterator, &
                                      Grid_releaseTileIterator
    use Grid_iterator,         ONLY : Grid_iterator_t
    use Grid_tile,             ONLY : Grid_tile_t
    use ut_testDriverMod

    implicit none

    integer, intent(in)    :: fileUnit
    logical, intent(inout) :: perfect

    !!!!! EXPECTED RESULTS BASED ON flash.par AND SETUP VALUES GIVEN ABOVE
    integer   :: NXCELL_EX   
    integer   :: NXBLK_EX    
    real      :: XMIN_EX     
    real      :: XMAX_EX     
    real      :: XDELTA_EX   
    integer   :: XL_BC_EX    
    integer   :: XH_BC_EX    
#if NDIM >= 2
    integer   :: NYCELL_EX   
    integer   :: NYBLK_EX    
    real      :: YMIN_EX     
    real      :: YMAX_EX      
    real      :: YDELTA_EX   
    integer   :: YL_BC_EX    
    integer   :: YH_BC_EX    
#endif
#if NDIM == 3
    integer   :: NZCELL_EX   
    integer   :: NZBLK_EX    
    real      :: ZMIN_EX     
    real      :: ZMAX_EX      
    real      :: ZDELTA_EX   
    integer   :: ZL_BC_EX     
    integer   :: ZH_BC_EX    
#endif
    integer,  parameter :: MAXLEVEL_EX =  1

    integer :: geometry
    real    :: domain(LOW:HIGH, MDIM)
    integer :: domainBC(LOW:HIGH, MDIM)
    real    :: deltas(1:MDIM)
    real    :: x_expected
#if NDIM >= 2
    real    :: y_expected
#endif
#if NDIM == 3
    real    :: z_expected
#endif
    integer :: max_level
    integer :: axis

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

    real, allocatable :: areas(:,:,:) 
    real, allocatable :: FaceAreasX(:,:,:)
    real, allocatable :: ComputedFaceAreasX(:,:,:)
#if NDIM >= 2
    real, allocatable :: FaceAreasY(:,:,:)
    real, allocatable :: ComputedFaceAreasY(:,:,:)
#endif
#if NDIM == 3
    real, allocatable :: FaceAreasZ(:,:,:)
    real, allocatable :: ComputedFaceAreasZ(:,:,:)
#endif
    real, allocatable :: volumes(:, :, :)
    
    integer :: offset(1:MDIM)
    integer :: lo(1:MDIM)
    integer :: hi(1:MDIM)

    integer :: rank
    integer :: ilev
    integer :: i, j, k, var

    interface
        subroutine get_coordinates(axis_min, daxis, cells, array)
                implicit none
                integer, intent(in) :: cells
                real, intent(in)    :: axis_min, daxis
                real, intent(out)   :: array(1:cells+1)
        end subroutine get_coordinates
    end interface

    nullify(solnData)

    call start_test_run

    ! Physical domain
    call Grid_getGeometry(geometry)
    if (geometry == CARTESIAN) then
            NXCELL_EX = 128
            NXBLK_EX  = 1
            XMIN_EX   = -64.0  
            XMAX_EX   = 64.0  
            XL_BC_EX  = OUTFLOW 
            XH_BC_EX  = OUTFLOW  
#if NDIM >= 2
            NYCELL_EX = 128
            NYBLK_EX  = 1
            YMIN_EX   = -64.0  
            YMAX_EX   = 64.0  
            YL_BC_EX  = OUTFLOW 
            YH_BC_EX  = OUTFLOW  
#endif
#if NDIM == 3
            NZCELL_EX = 128
            NZBLK_EX  = 1
            ZMIN_EX   = -64.0  
            ZMAX_EX   = 64.0  
            ZL_BC_EX  = OUTFLOW 
            ZH_BC_EX  = OUTFLOW  
#endif
    else if (geometry == CYLINDRICAL) then
            NXCELL_EX = 128
            NXBLK_EX  = 1
            XMIN_EX   = 0.0  
            XMAX_EX   = 128.0 
            XL_BC_EX  = REFLECTING 
            XH_BC_EX  = OUTFLOW  
#if NDIM >= 2
            NYCELL_EX = 128
            NYBLK_EX  = 1
            YMIN_EX   = -64.0 
            YMAX_EX   =  64.0 
            YL_BC_EX  = OUTFLOW 
            YH_BC_EX  = OUTFLOW  
#endif
#if NDIM == 3
            NZCELL_EX = 32
            NZBLK_EX  = 1
            ZMIN_EX   = 0.0 
            ZMAX_EX   = 2 * PI 
            ZL_BC_EX  = PERIODIC
            ZH_BC_EX  = PERIODIC
#endif
    else if (geometry == SPHERICAL) then
            NXCELL_EX = 128
            NXBLK_EX  = 1
            XMIN_EX   = 0.0  
            XMAX_EX   = 128.0  
            XL_BC_EX  = REFLECTING 
            XH_BC_EX  = OUTFLOW  
#if NDIM >= 2
            NYCELL_EX = 32
            NYBLK_EX  = 1
            YMIN_EX   = 0.0 
            YMAX_EX   = PI
            YL_BC_EX  = REFLECTING 
            YH_BC_EX  = REFLECTING  
#endif
#if NDIM == 3
            NZCELL_EX = 32
            NZBLK_EX  = 1
            ZMIN_EX   = 0.0 
            ZMAX_EX   = 2 * PI  
            ZL_BC_EX  = PERIODIC
            ZH_BC_EX  = PERIODIC  
#endif
    endif
           
            XDELTA_EX = (XMAX_EX-XMIN_EX)/NXCELL_EX
#if NDIM >= 2 
            YDELTA_EX = (YMAX_EX-YMIN_EX)/NYCELL_EX 
#endif
#if NDIM == 3
            ZDELTA_EX = (ZMAX_EX-ZMIN_EX)/NZCELL_EX  
#endif

    call Grid_getDomainBoundBox(domain)
    call assertEqual(domain(LOW,  IAXIS), XMIN_EX,"Incorrect low X-coordinate")
    call assertEqual(domain(HIGH, IAXIS), XMAX_EX,"Incorrect high X-coordinate")
#if NDIM>= 2
    call assertEqual(domain(LOW,  JAXIS), YMIN_EX,"Incorrect low Y-coordinate")
    call assertEqual(domain(HIGH, JAXIS), YMAX_EX,"Incorrect high Y-coordinate")
#endif
#if NDIM == 3
    call assertEqual(domain(LOW,  KAXIS), ZMIN_EX,"Incorrect low Z-coordinate")
    call assertEqual(domain(HIGH, KAXIS), ZMAX_EX,"Incorrect high Z-coordinate")
#endif
    !!!!! CONFIRM PROPER REFINEMENT
    call Grid_getMaxRefinement(max_level, mode=1)
    call assertEqual(max_level, MAXLEVEL_EX, "Incorrect maximum refine level")

    do ilev = 1, max_level
        call Grid_getDeltas(ilev, deltas)

        x_expected = XDELTA_EX / 2.0**(ilev - 1)
        call assertEqual(deltas(IAXIS),x_expected,"Incorrect delta of  X-coordinate")
#if NDIM >= 2
        y_expected = YDELTA_EX / 2.0**(ilev - 1)
        call assertEqual(deltas(JAXIS),y_expected,"Incorrect delta of  Y-coordinate")
#endif
#if NDIM == 3
        z_expected = ZDELTA_EX / 2.0**(ilev - 1)
        call assertEqual(deltas(KAXIS),z_expected,"Incorrect delta of Z-coordinate")
#endif
   
    end do

    !!!!! CONFIRM PROPER BLOCK/CELL STRUCTURE
    ! Walk across all blocks to test and collect info
    n_blocks = 0

    ! The tests that use this iterator are testing for block info
    ! Do not use tiling here
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

        call assertEqual(tileDesc%level, 1, "Incorrect block level")

        ! Check guard cells along all directions
        blkLimits   = tileDesc%limits
        blkLimitsGC = tileDesc%blkLimitsGC
        blkGC(LOW, :) = blkLimits(LOW, :) - blkLimitsGC(LOW, :)
        blkGC(HIGH, :) = blkLimitsGC(HIGH, :) - blkLimits(HIGH, :)
        call assertEqual(blkGC(LOW,  IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
        call assertEqual(blkGC(HIGH, IAXIS), NGUARD, &
                         "Incorrect guard cell along X-axis")
#if NDIM >= 2
        call assertEqual(blkGC(LOW,  JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
        call assertEqual(blkGC(HIGH, JAXIS), NGUARD, &
                         "Incorrect guard cell along Y-axis")
#endif
#if NDIM == 3
        call assertEqual(blkGC(LOW,  KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
        call assertEqual(blkGC(HIGH, KAXIS), NGUARD, &
                         "Incorrect guard cell along Z-axis")
#endif

        ! Correct cells per block along each direction
        blkSize = blkLimits(HIGH, :) - blkLimits(LOW, :) + 1
        call assertEqual(blkSize(IAXIS), NXCELL_EX / NXBLK_EX, &
                         "Incorrect cells per block along X-axis")
#if NDIM >= 2
        call assertEqual(blkSize(JAXIS), NYCELL_EX / NYBLK_EX, &
                         "Incorrect cells per block along Y-axis")
#endif
#if NDIM == 3
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
    call assertEqual(xBlkMin, 1, "Incorrect number of blocks along X-axis")
    call assertEqual(yBlkMin, 1, "Incorrect number of blocks along Y-axis")
    call assertEqual(zBlkMin, 1, "Incorrect number of blocks along Z-axis")

    call assertEqual(n_blocks, &
#if NDIM == 1
            NXBLK_EX * 1 * 1, &
#elif NDIM == 2
            NXBLK_EX * NYBLK_EX * 1, &
#else
            NXBLK_EX * NYBLK_EX * NZBLK_EX, &
#endif
                     "Incorrect total number of blocks")

    call assertEqual(xBlkMax, NXCELL_EX, &
                     "Incorrect total number of cells along X-axis")
#if NDIM >= 2
    call assertEqual(yBlkMax, NYCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
#endif
#if NDIM == 3
    call assertEqual(zBlkMax, NZCELL_EX, &
                     "Incorrect total number of cells along Y-axis")
#endif   

    call assertEqual(xMin, XMIN_EX, "Incorrect minimum X-coordinate found")
    call assertEqual(xMax, XMAX_EX, "Incorrect maximum X-coordinate found")
#if NDIM >= 2
    call assertEqual(yMin, YMIN_EX, "Incorrect minimum Y-coordinate found")
    call assertEqual(yMax, YMAX_EX, "Incorrect maximum Y-coordinate found")
#endif
#if NDIM == 3
    call assertEqual(zMin, ZMIN_EX, "Incorrect minimum Z-coordinate found")
    call assertEqual(zMax, ZMAX_EX, "Incorrect maximum Z-coordinate found")
#endif

    !!!!! CONFIRM PROPER BC
    call Grid_getDomainBC(domainBC)
    call assertEqual(domainBC(LOW,  IAXIS), XL_BC_EX, "Incorrect X-left BC")
    call assertEqual(domainBC(HIGH, IAXIS), XH_BC_EX, "Incorrect X-right BC")
#if NDIM >= 2
    call assertEqual(domainBC(LOW,  JAXIS), YL_BC_EX, "Incorrect Y-left BC")
    call assertEqual(domainBC(HIGH, JAXIS), YH_BC_EX, "Incorrect Y-right BC")
#endif
#if NDIM == 3
    call assertEqual(domainBC(LOW,  KAXIS), ZL_BC_EX, "Incorrect Z-left BC")
    call assertEqual(domainBC(HIGH, KAXIS), ZH_BC_EX, "Incorrect Z-right BC")
#endif

    !!!!! CONFIRM PROPER INITIAL CONDITIONS
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
        call itor%currentTile(tileDesc)
        call tileDesc%getDataPtr(solnData, CENTER)

        associate(lo => tileDesc%limits(LOW, :), &
                  hi => tileDesc%limits(HIGH, :))
            do           var = UNK_VARS_BEGIN, UNK_VARS_END 
                do         k = lo(KAXIS), hi(KAXIS)
                    do     j = lo(JAXIS), hi(JAXIS)
                        do i = lo(IAXIS), hi(IAXIS)
                            call assertEqual(solnData(var, i, j, k), &
                                             1.1 * var, &
                                             "Incorrect initial condition")
                        end do
                    end do
                end do
            end do
        end associate

        call tileDesc%releaseDataPtr(solnData, CENTER)

        call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    !!!!! CONFIRM CELL FACE AREAS ACCESSORS
    call Grid_getTileIterator(itor, LEAF, tiling=.TRUE.)
    do while (itor%isValid())
       call itor%currentTile(tileDesc)

       lo = tileDesc%limits(LOW,  :)
       hi = tileDesc%limits(HIGH, :)


       allocate(FaceAreasX(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(IAXIS, tileDesc%level, lo, hi, FaceAreasX)
       allocate(ComputedFaceAreasX(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       allocate(areas(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call cellfacearea_computed(IAXIS)

       ComputedFaceAreasX = areas
       deallocate(areas)
#if NDIM >=2
       allocate(FaceAreasY(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(JAXIS, tileDesc%level, lo, hi, FaceAreasY)
       allocate(ComputedFaceAreasY(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       allocate(areas(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call cellfacearea_computed(JAXIS)

       ComputedFaceAreasY = areas
       deallocate(areas)
#endif
#if NDIM == 3
       allocate(FaceAreasZ(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       allocate(ComputedFaceAreasZ(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call Grid_getCellFaceAreas(KAXIS, tileDesc%level, lo, hi, FaceAreasZ)

       allocate(areas(lo(IAXIS):hi(IAXIS), &
               lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)))
       call cellfacearea_computed(KAXIS)

       ComputedFaceAreasZ = areas
       deallocate(areas)
#endif

       do       k = lo(KAXIS), hi(KAXIS)
          do    j = lo(JAXIS), hi(JAXIS)
             do i = lo(IAXIS), hi(IAXIS)
                  print *, "(", i, j, k, ")"
                  call assertEqual(FaceAreasX(i,j,k), ComputedFaceAreasX(i,j,k), "Bad FaceArea along X-axis")
#if NDIM >= 2                  
                  call assertEqual(FaceAreasY(i,j,k), ComputedFaceAreasY(i,j,k), "Bad FaceArea along Y-axis")
                  !call assertAlmostEqual(FaceAreasY(i,j,k), ComputedFaceAreasY(i,j,k), 1.e-14, "Bad FaceArea along Y-axis")
#endif
#if NDIM == 3                 
                  call assertEqual(FaceAreasZ(i,j,k), ComputedFaceAreasZ(i,j,k), "Bad FaceArea along Z-axis")
#endif
                  end do
          end do
       end do

       deallocate(FaceAreasX)
       deallocate(ComputedFaceAreasX)
#if NDIM >= 2
       deallocate(FaceAreasY)
       deallocate(ComputedFaceAreasY)
#endif
#if NDIM == 3
       deallocate(FaceAreasZ)
       deallocate(ComputedFaceAreasZ)
#endif
       call itor%next()
    end do
    call Grid_releaseTileIterator(itor)

    perfect = finish_test_run()

contains

subroutine cellfacearea_computed(axis)
   use Driver_interface, ONLY : Driver_abort

   integer, intent(in) :: axis

   integer :: lo1(1:MDIM)
   integer :: hi1(1:MDIM) 
     
   real, allocatable :: rf(:)
#if NDIM >= 2
   real, allocatable :: thf(:)
#endif
   real    :: area, facebase
   integer :: i, j, k

   lo1(1:MDIM) = lo
   hi1(1:MDIM) = hi

   select case (geometry)
   case (CARTESIAN)
#if   NDIM == 1
         area = 1.0
#elif NDIM == 2
         if      (axis == IAXIS) then
            area = YDELTA_EX
         else if (axis == JAXIS) then
            area = XDELTA_EX
         else
            ! DEV: TODO Should this set area to 1?
            call Driver_abort("[Grid_getCellFaceAreas] Invalid axis for 2D")
         end if
#elif NDIM == 3
         if      (axis == IAXIS) then
            area = YDELTA_EX * ZDELTA_EX
         else if (axis == JAXIS) then
            area = XDELTA_EX * ZDELTA_EX
         else if (axis == KAXIS) then
            area = XDELTA_EX * YDELTA_EX
         end if
#endif

         do       k = lo1(KAXIS), hi1(KAXIS)
            do    j = lo1(JAXIS), hi1(JAXIS) 
               do i = lo1(IAXIS), hi1(IAXIS)
                  areas(i, j, k) = area
               end do
            end do
         end do
   case (CYLINDRICAL)
      if      (axis == IAXIS) then
         ! Get radial coordinate of face centers
         allocate(rf(1:NXCELL_EX+1))
         call get_coordinates(XMIN_EX, XDELTA_EX, NXCELL_EX, rf)
         
         do       k = lo1(KAXIS), hi1(KAXIS)
            do    j = lo1(JAXIS), hi1(JAXIS)
               do i = lo1(IAXIS), hi1(IAXIS)
#if   NDIM == 1
                     areas(i, j, k) = 2.0 * PI * ABS(rf(i))
#elif NDIM == 2
                     areas(i, j, k) = 2.0 * PI * ABS(rf(i)) * YDELTA_EX
#elif NDIM == 3
                     areas(i, j, k) = ABS(rf(i)) * YDELTA_EX * ZDELTA_EX
#endif
                end do
             end do
          end do
          deallocate(rf)
      else if (axis == JAXIS) then
         ! Get radial coordinate of face centers
         allocate(rf(1:NXCELL_EX+1))
         call get_coordinates(XMIN_EX, XDELTA_EX, NXCELL_EX, rf)

         do       k = lo1(KAXIS), hi1(KAXIS)
            do    j = lo1(JAXIS), hi1(JAXIS)
               do i = lo1(IAXIS), hi1(IAXIS)
#if   NDIM == 1
                     areas(i, j, k) = PI * (rf(i+1)**2 - rf(i)**2)
#elif NDIM == 2
                     areas(i, j, k) = PI * (rf(i+1)**2 - rf(i)**2)
#elif NDIM == 3
                     areas(i, j, k) = 0.5 * (rf(i+1)**2 - rf(i)**2) * ZDELTA_EX
#endif
               end do
            end do
         end do
         deallocate(rf)

      else if (axis == KAXIS) then
            do       k = lo1(KAXIS), hi1(KAXIS)
               do    j = lo1(JAXIS), hi1(JAXIS)
                  do i = lo1(IAXIS), hi1(IAXIS)
#if   NDIM == 1
                     areas(i, j, k) = XDELTA_EX
#elif NDIM == 2
                     areas(i, j, k) = XDELTA_EX * YDELTA_EX
#elif NDIM == 3
                     areas(i, j, k) = XDELTA_EX * YDELTA_EX
#endif
                  end do
               end do
            end do
      end if
   case (SPHERICAL)
      if      (axis == IAXIS) then
         ! Get coordinates of faces
         allocate(rf(1:NXCELL_EX+1))
         call get_coordinates(XMIN_EX, XDELTA_EX, NXCELL_EX, rf)

#if NDIM >= 2
         allocate(thf(1:NYCELL_EX+1))
         call get_coordinates(YMIN_EX, YDELTA_EX, NYCELL_EX, thf)
#endif

            do       k = lo1(KAXIS), hi1(KAXIS)
               do    j = lo1(JAXIS), hi1(JAXIS)
                  do i = lo1(IAXIS), hi1(IAXIS)
                     facebase = rf(i) * rf(i)
#if   NDIM == 1
                     areas(i, j, k) = facebase * 4.0 * PI
#elif NDIM == 2
                     areas(i, j, k) = facebase * ( cos(thf(j)) - cos(thf(j+1)) ) * 2.0 * PI
#elif NDIM == 3
                     areas(i, j, k) = facebase * ( cos(thf(j)) - cos(thf(j+1)) ) * ZDELTA_EX
#endif
                  end do
               end do
            end do
         deallocate(rf)
#if NDIM >= 2
         deallocate(thf)
#endif
      else if (axis == JAXIS) then
         ! Get coordinates of faces
         allocate(rf(NXCELL_EX+1))
         call get_coordinates(XMIN_EX, XDELTA_EX, NXCELL_EX, rf)

#if NDIM >= 2
         allocate(thf(NYCELL_EX+1))
         call get_coordinates(YMIN_EX, YDELTA_EX, NYCELL_EX, thf)
#endif

         do       k = lo1(KAXIS), hi1(KAXIS)
            do    j = lo1(JAXIS), hi1(JAXIS)
               do i = lo1(IAXIS), hi1(IAXIS)

                     facebase = (rf(i)+rf(i+1))*(rf(i+1)-rf(i))*0.5
#if   NDIM == 1
                     areas(i, j, k) = facebase * 2.0 * PI
#elif NDIM == 2
                     areas(i, j, k) = facebase * sin(thf(j)) * 2.0 * PI
#elif NDIM == 3
                     areas(i, j, k) = facebase * sin(thf(j)) * ZDELTA_EX
#endif
               end do
            end do
         end do
         deallocate(rf)
#if NDIM >= 2
         deallocate(thf)
#endif
      else if (axis == KAXIS) then
         ! Get coordinates of faces
         allocate(rf(NXCELL_EX+1))
         call get_coordinates(XMIN_EX, XDELTA_EX, NXCELL_EX, rf)

#if NDIM >= 2
         allocate(thf(NYCELL_EX+1))
         call get_coordinates(YMIN_EX, YDELTA_EX, NYCELL_EX, thf)
#endif

         do       k = lo1(KAXIS), hi1(KAXIS)
            do    j = lo1(JAXIS), hi1(JAXIS)
               do i = lo1(IAXIS), hi1(IAXIS)

                     facebase = 0.5 * (rf(i+1) + rf(i)) * (rf(i+1) - rf(i))
#if NDIM == 1
                     areas(i, j, k) = facebase * PI
#elif NDIM >= 2
                     areas(i, j, k) = facebase * YDELTA_EX
#endif
               end do
            end do
         end do
         deallocate(rf)
#if NDIM >= 2
         deallocate(thf)
#endif
      end if
   end select
end subroutine cellfacearea_computed
end subroutine Grid_unitTest

subroutine get_coordinates(axis_min, daxis, cells, array)
        implicit none
        integer, intent(in) :: cells
        real, intent(in)    :: axis_min, daxis
        real, intent(out)   :: array(1:cells+1)

        integer             :: i


        do i=1, cells+1
           array(i) = axis_min + (i - 1) * daxis
        end do
end subroutine get_coordinates
