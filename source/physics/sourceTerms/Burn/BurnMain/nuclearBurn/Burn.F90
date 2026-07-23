!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn
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
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( real, intent(IN) ::  dt  )
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   dt  --       passed to the internal bn_burner module
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useBurnTable -- Boolean, False.  Controls the generation of reaction rates.
!!                TRUE interpolates from a stored table; FALSE generates them
!!                analytically.
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  algebra -- Integer, 1, [1,2].  Controls choice of linear algebra package used
!!                for matrix solution.  1=Ma28 sparse package, 2=Gift hardwired package.
!!  odeStepper -- Integer, 1, [1,2].  Controls time integration routines.
!!                1=Bader-Deuflhard variable order, 2=Rosenbrock 4th order
!!  nuclearTempMin/Max -- Real, 1.1E+8/1.0E+12.  Minimum and maximum temperature
!!                ranges where burning can occur
!!  nuclearDensMin/Max -- Real, 1.0E-10/1.0E+14.  Minimum and maximum density range
!!                where burning can occur.
!!  nuclearNI56Max -- Real, 1.0.  Maximum mass fraction of nickel where burning
!!                can occur.
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.See Burn_computeDt for details.
!!
!! NOTES
!!
!!  The burning unit adds a new mesh variable ENUC_VAR which is the nuclear energy
!!             generation rate
!!
!!***

!!REORDER(4): solnData

subroutine Burn (  dt  )

  use Burn_data, ONLY : bn_useBurn, bn_useShockBurn, bn_gcMaskSD
  use bn_interface, ONLY : bn_mapNetworkToSpecies, bn_burner

  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getCellCoords, &
       Grid_getMaxRefinement, &
       Grid_getTileIterator, Grid_releaseTileIterator
  use Eos_interface, ONLY : Eos_multiDim
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Hydro_interface, ONLY : Hydro_shockStrength

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use Burn_interface, ONLY : Burn_burner, Burn_update

  implicit none

#define DEBUG_GRID_GCMASK
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

  !args
  real, intent(in) :: dt

  ! locals
  integer :: i, j, k, n, speciesMap

  real, dimension(NSPECIES) :: xIn, xOut
  real :: sdot, tmp, rho, ei, ek, enuc

  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel

  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  integer, dimension(1:MDIM) :: lo,hi,loGC,hiGC

  integer, parameter :: shock_mode = 1
  real, parameter :: shock_thresh = 0.33
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,parameter :: gcMaskLogged =.TRUE.
#endif

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits
  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc


  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

#ifndef SHOK_VAR
  call Driver_abort("[Burn.F90] SHOK_VAR is not defined")
#endif

  !---------------------------------------------------------------------------------
  nullify(solnData)
  ! start the timer ticking
  call Timers_start("burn")

  if (.NOT. bn_useShockBurn) then
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(bn_gcMaskSD, .FALSE., '[Burn]', 'gcWant[Detect]')
     end if
#endif
     call Grid_fillGuardCells(CENTER, ALLDIR, maskSize=NUNK_VARS, mask=bn_gcMaskSD,&
                              doLogMask=.NOT.gcMaskLogged)
#ifdef DEBUG_GRID_GCMASK
     gcMaskLogged = .TRUE.
#endif
  endif

  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE. )
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     ! get dimensions/limits and coordinates
     lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
     hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)
     loGC(1:MDIM) = tileDesc%grownLimits(LOW,1:MDIM)
     hiGC(1:MDIM) = tileDesc%grownLimits(HIGH,1:MDIM)

     blkLimits = tileDesc%limits(:,:)

     ! allocate space for dimensions
     allocate(xCoord(loGC(IAXIS):hiGC(IAXIS)))
     allocate(yCoord(loGC(JAXIS):hiGC(JAXIS)))
     allocate(zCoord(loGC(KAXIS):hiGC(KAXIS)))

     call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level,loGC,hiGC,xCoord)
     call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level,loGC,hiGC,yCoord)
     call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level,loGC,hiGC,zCoord)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Shock detector
     if (.NOT. bn_useShockBurn) then
        solnData(SHOK_VAR,:,:,:)=0.0
        call Hydro_shockStrength(solnData, lo,hi, loGC, hiGC, &
             xCoord,yCoord,zCoord,shock_thresh,shock_mode)
     endif

     call Burn_burner(solnData, loGC, blkLimits, dt)   ! NOTE: this may flag `burnedZone`

     call Burn_update(solnData, loGC, blkLimits, dt)

     ! ! we've altered the EI, let's equilabrate
     ! if (burnedZone) then
     !    lo(:) = tileDesc%blkLimitsGC(LOW,:)
     !    call Eos_multiDim(MODE_DENS_EI,tileDesc%limits,lo,solnData)
     ! end if

     call tileDesc%releaseDataPtr(solnData, CENTER)
     nullify(solnData)

     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)

     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)

  call Timers_stop("burn")

  return
end subroutine Burn
