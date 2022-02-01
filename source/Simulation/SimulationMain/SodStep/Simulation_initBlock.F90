!!****if* source/Simulation/SimulationMain/SodStep/Simulation_initBlock
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
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  This version sets up the Sod shock-tube
!!  problem.
!!
!!  Reference: Sod, G. A., 1978, J. Comp. Phys., 27, 1
!!
!! 
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!
!!***

!!REORDER(4): solnData
subroutine Simulation_initBlock(solnData, tileDesc)

#include "constants.h"
#include "Simulation.h"

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getCellCoords

  use Simulation_data, ONLY: sim_posn, sim_xCos, sim_yCos, sim_zCos, &    
     &  sim_rhoLeft,  sim_pLeft, sim_uLeft, sim_rhoRight, sim_pRight, sim_uRight, &
     &  sim_smallX, sim_gamma, sim_smallP

  implicit none

  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: i, j, k, n
  integer :: iMax, jMax, kMax

  real :: xx, yy,  zz, xxL, xxR
  
  real :: lPosn0, lPosn
  

  real,allocatable, dimension(:) ::xCenter,xLeft,xRight,yCoord,zCoord

  real :: rhoZone, velxZone, velyZone, velzZone, presZone, & 
       eintZone, enerZone, ekinZone

  integer :: lo(1:MDIM)
  integer :: hi(1:MDIM)

  lo(:) = tileDesc%limits(LOW,  :)
  hi(:) = tileDesc%limits(HIGH, :)
  allocate(  xLeft(lo(IAXIS):hi(IAXIS)))
  allocate( xRight(lo(IAXIS):hi(IAXIS)))
  allocate(xCenter(lo(IAXIS):hi(IAXIS)))
  allocate( yCoord(lo(JAXIS):hi(JAXIS)))
  allocate( zCoord(lo(KAXIS):hi(KAXIS)))
  xLeft = 0.0
  xRight = 0.0
  xCenter = 0.0
  yCoord = 0.0
  zCoord = 0.0

#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          lo, hi, zCoord)
#endif
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          lo, hi, yCoord)
#endif

  call Grid_getCellCoords(IAXIS, LEFT_EDGE, tileDesc%level, &
                          lo, hi, xLeft)
  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          lo, hi, xCenter)
  call Grid_getCellCoords(IAXIS, RIGHT_EDGE, tileDesc%level, &
                          lo, hi, xRight)

#ifdef DEBUG_SIMULATION
98 format('initBlock:',A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format('initBlock:',A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print 99,"solnData" ,(lbound(solnData ,i),ubound(solnData ,i),i=1,4)
  print*,'tile limits:',tileDesc%limits
  print*,'grown tile limits:',tileDesc%limitsGC
#endif
!------------------------------------------------------------------------------

! Loop over cells in the block.  For each, compute the physical position of 
! its left and right edge and its center as well as its physical width.  
! Then decide which side of the initial discontinuity it is on and initialize 
! the hydro variables appropriately.

  do k = lo(KAXIS), hi(KAXIS)

     ! get the coordinates of the cell center in the z-direction
     zz = zCoord(k)
     
     ! Where along the x-axis the shock intersects the xz-plane at the current z.
     lPosn0 = sim_posn - zz*sim_zCos/sim_xCos
     
     do j = lo(JAXIS), hi(JAXIS)
        
        ! get the coordinates of the cell center in the y-direction
        yy = yCoord(j)

        ! The position of the shock in the current yz-row.
        lPosn = lPosn0 - yy*sim_yCos/sim_xCos
        
        do i = lo(IAXIS), hi(IAXIS)
           
           ! get the cell center, left, and right positions in x
           xx  = xCenter(i)
           
           xxL = xLeft(i)
           xxR = xRight(i)

           ! initialize cells to the left of the initial shock.
           if (xxR <= lPosn) then
              
              presZone = sim_pLeft
              rhoZone = sim_rhoLeft
              
              velxZone = sim_uLeft * sim_xCos
              velyZone = sim_uLeft * sim_yCos
              velzZone = sim_uLeft * sim_zCos 
              
              ! initialize cells which straddle the shock.  Treat them as though 1/2 of 
              ! the cell lay to the left and 1/2 lay to the right.
           elseif ((xxL < lPosn) .and. (xxR > lPosn)) then
              
              presZone = 0.5 * (sim_pLeft+sim_pRight)
              rhoZone = 0.5 * (sim_rhoLeft+sim_rhoRight)

              velxZone = 0.5 *(sim_uLeft+sim_uRight) * sim_xCos
              velyZone = 0.5 *(sim_uLeft+sim_uRight) * sim_yCos
              velzZone = 0.5 *(sim_uLeft+sim_uRight) * sim_zCos
              
              ! initialize cells to the right of the initial shock.
           else

              presZone = sim_pRight
              rhoZone = sim_rhoRight

              velxZone = sim_uRight * sim_xCos
              velyZone = sim_uRight * sim_yCos
              velzZone = sim_uRight * sim_zCos
              
           endif

#if NSPECIES > 0
           !put in value of default species
           solnData(SPECIES_BEGIN,i,j,k)=1.0-(NSPECIES-1)*sim_smallX

              !if there is only 1 species, this loop will not execute
           do n = SPECIES_BEGIN+1,SPECIES_END
              solnData(n,i,j,k)= sim_smallX
           enddo
#endif

           ! Compute the gas energy and set the gamma-values needed for the equation of 
           ! state.
           ekinZone = 0.5 * (velxZone**2 + & 
                velyZone**2 + & 
                velzZone**2)
           
           eintZone = presZone / (sim_gamma-1.)
           eintZone = eintZone / rhoZone
           enerZone = eintZone + ekinZone
           enerZone = max(enerZone, sim_smallP)
           
           ! Store the variables in the current zone.
           ! Data is stored one cell at a time through the solnData pointer.


           solnData(DENS_VAR, i,j,k) =  rhoZone
           solnData(PRES_VAR, i,j,k) =  presZone
           solnData(VELX_VAR, i,j,k) =  velxZone
           solnData(VELY_VAR, i,j,k) =  velyZone
           solnData(VELZ_VAR, i,j,k) =  velzZone

#ifdef ENER_VAR
           solnData(ENER_VAR, i,j,k) =  enerZone
#endif
#ifdef GAME_VAR          
           solnData(GAME_VAR, i,j,k) =  sim_gamma
#endif
#ifdef GAMC_VAR
           solnData(GAMC_VAR, i,j,k) =  sim_gamma
#endif


        enddo
     enddo
  enddo
 
  deallocate(xLeft)
  deallocate(xRight)
  deallocate(xCenter)
  deallocate(yCoord)
  deallocate(zCoord)

  return
end subroutine Simulation_initBlock










