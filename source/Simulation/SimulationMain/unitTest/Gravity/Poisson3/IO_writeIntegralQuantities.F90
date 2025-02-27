!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/IO_writeIntegralQuantities
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
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities() 
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
!!   quantities other than default values in the flashx.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName
  use Grid_interface, ONLY : Grid_getCellVolumes, &
    Grid_getCellCoords, Grid_getTileIterator, Grid_releaseTileIterator
  use Simulation_data, ONLY: sim_xctr, sim_yctr, sim_zctr, sim_a1, sim_a3, &
       sim_a1inv, sim_a3inv, sim_initGeometry

  use IO_data, ONLY : io_globalMe
  use Grid_tile, ONLY: Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
  
  implicit none

#include "mpif.h"
#include "constants.h"
#include "Simulation.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: lo(MDIM), hi(MDIM)
  real, dimension(:), allocatable :: xCenter, yCenter, zCenter
  real, dimension(:,:,:), allocatable :: cellVolumes 
  real        :: radius2, xdist, ydist, zdist
  real        :: radiusInside2, a3x90, a3x90Inv

  integer, parameter ::  nGlobalSum = 11          ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  integer :: i, j, k, ii, jj, kk
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  real  :: dvol

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tileDesc

  NULLIFY(solnData)
  ! define a variable that is 90% of the polar radius
  a3x90 = 0.9*sim_a3
  a3x90Inv = 0.9*sim_a3inv

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while (itor%isValid())
     call itor%currentTile(tileDesc)
     lo(:) = tileDesc%limits(LOW,:)
     hi(:) = tileDesc%limits(HIGH,:)
     allocate(xCenter(lo(IAXIS):hi(IAXIS)))
     allocate(yCenter(lo(JAXIS):hi(JAXIS)))
     allocate(zCenter(lo(KAXIS):hi(KAXIS)))
     allocate(cellVolumes(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)))

     ! get a pointer to the current block of data
     call tileDesc%getDataPtr(solnData, CENTER)

     !! Get the cell coordinates
     call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level, lo, hi, xCenter)
     call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level, lo, hi, yCenter)
     call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level, lo, hi, zCenter)


     !! Get the cell volume for all cells
     call Grid_getCellVolumes(tileDesc%level, lo, hi, cellVolumes)

     ! Sum contributions from the indicated blkLimits of cells.
     do       k = lo(KAXIS), hi(KAXIS)
        zdist = (zCenter(k) - sim_zctr)
        do    j = lo(JAXIS), hi(JAXIS)
           ydist = (yCenter(j) - sim_yctr)
           do i = lo(IAXIS), hi(IAXIS)
              xdist = (xCenter(i) - sim_xctr)

              dvol = cellVolumes(i,j,k)
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
#endif ! ifdef DENS_VAR

!! Specialist information for the Maclaurin spheriod case

              lsum(8) = lsum(8) + solnData(ERRN_VAR,i,j,k)**2

              !! Divide inside and outside spheriod, as well as inside a sphere 
              !!   of radius 0.9*a3 (completely inside sphereoid)
              select case (sim_initGeometry)
              case (CYLINDRICAL)   ! 2d axisymmetric
                 radius2 = (xdist*sim_a1inv)**2 + (ydist*sim_a3inv)**2
                 radiusInside2 = (xdist*a3x90Inv)**2 + (ydist*a3x90Inv)**2
              case (CARTESIAN)       ! 3d cartesian
                 radius2 = (xdist*sim_a1inv)**2 + (ydist*sim_a1inv)**2 + (zdist*sim_a3inv)**2
                 radiusInside2 = (xdist*a3x90Inv)**2 + (ydist*a3x90Inv)**2 + (zdist*a3x90Inv)**2
              end select
              if (radius2 < 1.0) then     !! inside the spheroid
                 lsum(9) = lsum(9) + solnData(ERRN_VAR,i,j,k)**2
              else                        !! outside the spheriod
                 lsum(10) = lsum(10) + solnData(ERRN_VAR,i,j,k)**2
              end if
              if (radiusInside2 < 1.0) then   !! completely enclosed by spheriod
                 lsum(11) = lsum(11) + solnData(ERRN_VAR,i,j,k)**2
              end if
           enddo
        enddo
     enddo

     call tileDesc%releaseDataPtr(solnData, CENTER)
     deallocate(xCenter,yCenter,zCenter,cellVolumes)
     call itor%next()
  enddo
  call Grid_releaseTileIterator(itor)
 
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, MPI_Double_Precision, MPI_Sum, & 
       &                MASTER_PE, MPI_Comm_World, error)
  
  ! Take sqrts for L2 norms
  lsum(8) = sqrt(lsum(8))
  lsum(9) = sqrt(lsum(9))
  lsum(10) = sqrt(lsum(10))
  lsum(11) = sqrt(lsum(11))

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     if (isfirst == 0) then
        open (funit, file=trim(io_statsFileName), position='APPEND')
     else 
        if (.NOT. io_restart) then
           open (funit, file=trim(io_statsFileName)) 
           write (funit, 10)               &
                '#time                     ', &
                'mass                      ', &
                'x-momentum                ', &
                'y-momentum                ', & 
                'z-momentum                ', &
                'E_total                   ', &
                'E_kinetic                 ', &
                'E_internal                ', &
                'L2_total                  ', &
                'L2_insideSpheriod         ', &
                'L2_outsideSpheriod        ', &
                'L2_insidePoles        '
           
           
10         format (2x,50(a25, :, 1X))
           
        else
           open (funit, file=trim(io_statsFileName), position='APPEND')
           write (funit, 11) 
11         format('# simulation restarted')
        endif
     endif
     
     write (funit, 12) simtime, gsum      ! Write the global sums to the file.
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
    
  endif
  
  call MPI_Barrier (MPI_Comm_World, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities



