!! Display 1
!! This is how we iterate over blocks in FLASH. 

  integer, dimension(MAXBLOCKS) :: blkList
  integer :: lb, count

  call Grid_getListOfBlocks(LEAF, blkList, count)
  do lb = 1, count

     ! Do something with the current block !

  end do


!! Display 2


  integer, dimension(MAXBLOCKS) :: blkList
  integer :: lb, count, blkID

  call Grid_getListOfBlocks(LEAF, blkList, count)
  do lb = 1, count
     blkID = blkList(lb)

     ! Do something with the block with this ID !

  end do


!! Display 3


  integer, dimension(MAXBLOCKS) :: blkList
  integer :: lb, count, blkID
  real, pointer:: solnData(:,:,:,:)

  call Grid_getListOfBlocks(LEAF, blkList, count)
  do lb = 1, count
     blkID = blkList(lb)
     call Grid_getBlkPtr(blkID, solnData)

     ! Do something with the data in solnData !

     call Grid_releaseBlkPtr(blkList(i), solnData)
  end do


!! Display 4
!! For example, this is how solution data are initialized:


  integer, dimension(MAXBLOCKS) :: blkList
  integer :: lb, count, blkID
  real, pointer:: solnData(:,:,:,:)

  call Grid_getListOfBlocks(LEAF, blkList, count)
  do lb = 1, count
     blkID = blkList(lb)
     call Grid_getBlkPtr(blkID, solnData)

     ! Call the user function to fill in initial data:
     call Simulation_initBlock(blockID)

     call Grid_releaseBlkPtr(blkList(i), solnData)
  end do


!! Display 5
!! We can get metainformation and operate on block data directly:


  integer, dimension(MAXBLOCKS) :: blkList
  integer :: lb, count, blkID
  real, pointer:: solnData(:,:,:,:)
  integer, dimension(2,MDIM) :: blkLimitsGC

  call Grid_getListOfBlocks(LEAF, blkList, count)
  do lb = 1, count
     blkID = blkList(lb)
     call Grid_getBlkIndexLimits(blkID, blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkID, solnData)

     do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
        do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
           do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
              solnData(DENS_VAR,i,j,k) = 0.001
           end do
        end do
     end do

     call Grid_releaseBlkPtr(blkList(i), solnData)
  end do

!! The solnData array always includes room for guard cells and interior cells.
!!
!! A local range for the indexes is used.
!! The data for each block are as if declared by (example: 2D case)
!!    real,dimension(NUNK_VARS,NXB+2*NGUARD,NYB+2*NGUARD,1) :: solnData


!!!! LEFTOVER !!!!
do k=1,NZB
   do j=1,NYB
      do i=1,NXB
         solnData(DENS_VAR,i,j,k) = 0.001
      end do
   end do
end do
