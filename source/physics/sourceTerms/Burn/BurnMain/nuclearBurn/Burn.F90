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

  use Burn_data, ONLY : bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_useBurn, bn_gcMaskSD
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

  logical :: burnedZone
  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel

  real, allocatable, dimension(:) :: xCoord, yCoord, zCoord
  integer, dimension(1:MDIM) :: lo,hi,loHalo,hiHalo

  integer, parameter :: shock_mode = 1
  real, parameter :: shock_thresh = 0.33
  real, allocatable, dimension(:,:,:) :: shock
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,parameter :: gcMaskLogged =.TRUE.
#endif

  real, pointer, dimension(:,:,:,:) :: solnData

  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

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

     burnedZone = .FALSE.

     ! get dimensions/limits and coordinates
     lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
     hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)
     loHalo(1:MDIM) = tileDesc%grownLimits(LOW,1:MDIM)
     hiHalo(1:MDIM) = tileDesc%grownLimits(HIGH,1:MDIM)


     ! allocate space for dimensions
     allocate(xCoord(loHalo(IAXIS):hiHalo(IAXIS)))
     allocate(yCoord(loHalo(JAXIS):hiHalo(JAXIS)))
     allocate(zCoord(loHalo(KAXIS):hiHalo(KAXIS)))

     allocate(shock(loHalo(IAXIS):hiHalo(IAXIS),&
                    loHalo(JAXIS):hiHalo(JAXIS),&
                    loHalo(KAXIS):hiHalo(KAXIS)))

     call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level,loHalo,hiHalo,xCoord)
     call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level,loHalo,hiHalo,yCoord)
     call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level,loHalo,hiHalo,zCoord)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Shock detector
     if (.NOT. bn_useShockBurn) then
        call Hydro_shockStrength(solnData, shock, lo,hi, loHalo, hiHalo, &
             xCoord,yCoord,zCoord,shock_thresh,shock_mode)
     else
        shock(:,:,:) = 0.0
     endif

     ! AH: Aprox13 and Aprox19 are not currently thread-safe
     !!$omp parallel do &
     !!$omp collapse(3) &
     !!$omp default(shared) &
     !!$omp private(i,j,k,n,speciesMap,tmp,rho,sdot,xIn,xOut,ei,ek,enuc, &
     !!$omp         okBurnTemp,okBurnDens,okBurnShock,okBurnNickel)
     do k = lo(KAXIS), hi(KAXIS)
        do j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)

              tmp  = solnData(TEMP_VAR,i,j,k)
              rho  = solnData(DENS_VAR,i,j,k)
              sdot = 0.0e0

              okBurnTemp = (tmp >= bn_nuclearTempMin .AND. tmp <= bn_nuclearTempMax)
              okBurnDens = (rho >= bn_nuclearDensMin .AND. rho <= bn_nuclearDensMax)
              okBurnShock = (shock(i,j,k) <= 0.0 .OR. (shock(i,j,k) > 0.0 .AND. bn_useShockBurn))

              if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

                 if (NI56_SPEC /= NONEXISTENT) then
                    okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                 else    ! nickel is not even a species in this simulation, so we'll always burn
                    okBurnNickel = .TRUE.
                 endif

                 if (okBurnNickel) then

                    !$omp atomic write
                    burnedZone = .TRUE.
                    !$omp end atomic

                    ! Map the solution data into the order required by bn_burner
                    do n = 1, NSPECIES
                       call bn_mapNetworkToSpecies(n,speciesMap)
                       xIn(n) = solnData(speciesMap,i,j,k)
                    end do

                    ! Do the actual burn
                    call bn_burner(dt, tmp, rho, xIn, xOut, sdot)

                    !  Map it back NOTE someday make a nicer interface....
                    do n = 1, NSPECIES
                       call bn_mapNetworkToSpecies(n,speciesMap)
                       solnData(speciesMap,i,j,k) = xOut(n)
                    end do

                    ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                         solnData(VELY_VAR,i,j,k)**2 +  &
                         solnData(VELZ_VAR,i,j,k)**2)

                    ! internal energy, add on nuclear rate*timestep
                    enuc = dt*sdot
                    ei = solnData(ENER_VAR,i,j,k) + enuc - ek

#ifdef EINT_VAR
                    solnData(EINT_VAR,i,j,k) = ei
#endif
                    solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                    solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + dt*sdot
#endif
                    solnData(ENUC_VAR,i,j,k) = sdot

                 endif
              endif

           end do
        end do
     end do
     !!$omp end parallel do

     ! we've altered the EI, let's equilabrate
     if (burnedZone) then

        call Eos_multiDim(MODE_DENS_EI,tileDesc%limits,solnData)


     end if

     call tileDesc%releaseDataPtr(solnData, CENTER)
     nullify(solnData)

     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
     deallocate(shock)

     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)

  call Timers_stop("burn")

  return
end subroutine Burn
