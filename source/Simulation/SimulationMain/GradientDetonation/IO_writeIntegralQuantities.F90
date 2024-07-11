!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: myPE, 
!!                                    integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   myPE - current processor number
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities (isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, &
                        io_globalMe
  use Grid_interface, ONLY : Grid_getSingleCellVol, &
     Grid_getTileIterator, Grid_releaseTileIterator
  use Simulation_data, ONLY : sim_tempAmbient
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

  !implicit none

#include "Flashx_mpi_implicitNone.fh"
#include "constants.h"
#include "Simulation.h"
  
  !integer, intent(in) :: myPE
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  !integer :: blockList(MAXBLOCKS)

  !integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  integer, parameter ::  nGlobalSum = 7          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k
  real :: dvol             !, del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

  real, allocatable :: xCoord(:)
  integer :: iSize
  real :: postshockmaxx, postshockmaxx_global
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tileDesc
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  !call Grid_getListOfBlocks(LEAF, blockList, count)
  call Grid_getTileIterator(itor, nodetype=LEAF)
  postshockmaxx = -HUGE(1.0e0)
  
  !do lb = 1, count
  do while(itor%isValid())
     !get the index limits of the block
     call itor%currentTile(tileDesc)
     !call Grid_getBlkIndexLimits(blockList(lb), blkLimits, blkLimitsGC)
     
     ! get x coordinate for locating shock front
     !iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     !allocate(xCoord(iSize))
     grownTileLimits = tileDesc%grownLimits
     allocate(xCoord(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS)))
     !call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, .true., xCoord, iSize)
     call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
               grownTileLimits(LOW,  :), grownTileLimits(HIGH,  :), xCoord)
     ! get a pointer to the current block of data
     !call Grid_getBlkPtr(blockList(lb), solnData)
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = tileDesc%limits(LOW,KAXIS),tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(point, tileDesc%level, dvol)
     
              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif !ifdef DENS_VAR

              ! find maximum x coordinate at which there is post-shock material
              ! post shock is temperature more than 0.1% above ambient temperature
              if ( solnData(TEMP_VAR,i,j,k) > 1.001*sim_tempAmbient ) then
                 if (xCoord(i) > postshockmaxx) postshockmaxx = xCoord(i)
              endif

           enddo
        enddo
     enddo
     !call Grid_releaseBlkPtr(blockList(lb), solnData)
     call tileDesc%getDataPtr(solnData, CENTER)

     deallocate(xCoord)
     call itor%next()
  enddo
  call Grid_releaseTileIterator(itor)

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, MPI_COMM_WORLD, error)

  call MPI_Reduce (postshockmaxx, postshockmaxx_global, 1, FLASH_REAL, MPI_MAX, &
       &                MASTER_PE, MPI_COMM_WORLD, error)
  

  if (io_globalMe == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     
     !No mater what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif
     
     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'Max x post-shock matter   '
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     
     write (funit, 12) simtime, gsum, postshockmaxx_global      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



