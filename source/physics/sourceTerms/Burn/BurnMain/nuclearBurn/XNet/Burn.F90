!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/Burn
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

#include "Simulation.h"

subroutine Burn (  dt  )

  use bn_interface, ONLY : bn_burner, bn_azbar
  use bn_xnetData, ONLY : xnet_myid, xnet_nzbatchmx, xnet_inuc2unk
  use Burn_data, ONLY : bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &   bn_nuclearDensMax, bn_nuclearNI56Max, bn_useShockBurn, &
       &   bn_useBurn, bn_gcMaskSD, xmass
  use Burn_dataEOS, only: ytot1, bye
  use Driver_interface, ONLY : Driver_abort
  use Eos_interface, ONLY : Eos_multiDim
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_getCellCoords, &
       Grid_getTileIterator, Grid_releaseTileIterator
  use Logfile_interface, ONLY : Logfile_stampVarMask
  use Hydro_interface, ONLY : Hydro_shockStrength
  use Timers_interface, ONLY : Timers_start, Timers_stop

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

  use Grid_interface, ONLY : Grid_setWork

  !$ use omp_lib

  implicit none

#define DEBUG_GRID_GCMASK
#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

  !args
  real, intent(in) :: dt

  ! locals
  integer :: blockID, thisBlock, blockCount
  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:)   :: xCoord, yCoord, zCoord

  integer :: iSize, jSize, kSize, iSizeGC, jSizeGC, kSizeGC
  integer :: iSize_max, jSize_max, kSize_max
  integer,dimension(1:MDIM) :: lo, hi, loGC,hiGC, loHalo, hiHalo
  logical :: okBurnTemp, okBurnDens, okBurnShock, okBurnNickel
  logical, parameter :: getGuardCells = .true.

  real,    allocatable :: shock(:,:,:)
  real,    allocatable, target :: xIn(:,:,:,:,:), xOut(:,:,:,:,:)
  real,    allocatable, target :: sdot(:,:,:,:), tmp(:,:,:,:), rho(:,:,:,:)
  logical, allocatable, target :: burnedZone(:,:,:,:)
  integer, allocatable, target :: zone(:,:,:,:,:)

  real,    pointer :: xIn_batch(:,:,:), xOut_batch(:,:,:)
  real,    pointer :: sdot_batch(:,:), tmp_batch(:,:), rho_batch(:,:)
  logical, pointer :: burnedZone_batch(:,:)
  integer, pointer :: zone_batch(:,:,:)

  integer :: nzones, batchCount
  integer, dimension(:), allocatable :: sumBurn_TS_batch
  integer, dimension(:), allocatable :: batch_lo, batch_hi
  integer, dimension(:), allocatable :: sumBurn_TS

  integer, parameter :: shock_mode = 1
  real, parameter :: shock_thresh = 0.33
  real :: ei, ek, enuc
  integer :: i, j, k, m, n, ii, jj, kk, mm, nn
#ifdef DEBUG_GRID_GCMASK
  logical,save :: gcMaskLogged =.FALSE.
#else
  logical,parameter :: gcMaskLogged =.TRUE.
#endif

  type(Grid_iterator_t)  :: itor
  type(Grid_tile_t) :: tileDesc

  ! ----------------------- check if burning is requested in runtime parameters -------
  if (.not. bn_useBurn) return

  !---------------------------------------------------------------------------------
  nullify(solnData)
  ! start the timer ticking
  call Timers_start("burn")

  call Timers_start("burn_top")

  blockCount = 0
  iSize_max = 1
  jSize_max = 1
  kSize_max = 1
  call Grid_getTileIterator(itor, LEAF, tiling=.false.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     blockCount = blockCount+1
     lo(:) = tileDesc%limits(LOW,:)
     hi(:) = tileDesc%limits(HIGH,:)
     iSize = hi(IAXIS)-lo(IAXIS)+1
     jSize = hi(JAXIS)-lo(JAXIS)+1
     kSize = hi(KAXIS)-lo(KAXIS)+1
     iSize_max = max(iSize_max,iSize)
     jSize_max = max(jSize_max,jSize)
     kSize_max = max(kSize_max,kSize)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

  allocate(xIn(NSPECIES,iSize_max,jSize_max,kSize_max,blockCount))
  allocate(xOut(NSPECIES,iSize_max,jSize_max,kSize_max,blockCount))
  allocate(sdot(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(tmp(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(rho(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(burnedZone(iSize_max,jSize_max,kSize_max,blockCount))
  allocate(zone(4,iSize_max,jSize_max,kSize_max,blockCount))
  allocate(batch_lo(blockCount))
  allocate(batch_hi(blockCount))
  allocate(sumBurn_TS(blockCount))

  burnedZone = .FALSE.

  if (.NOT. bn_useShockBurn) then
#ifdef DEBUG_GRID_GCMASK
     if (.NOT.gcMaskLogged) then
        call Logfile_stampVarMask(bn_gcMaskSD, .FALSE., '[Burn]', 'gcWant[Shock]')
     end if
#endif
     call Grid_fillGuardCells(CENTER, ALLDIR, maskSize=NUNK_VARS, mask=bn_gcMaskSD,&
                              doLogMask=.NOT.gcMaskLogged)
#ifdef DEBUG_GRID_GCMASK
     gcMaskLogged = .TRUE.
#endif
  endif

  nzones = 0
  thisBlock = 0
  call Grid_getTileIterator(itor, LEAF, tiling=.false. )
  do while(itor%isValid())
     call itor%CurrentTile(tileDesc)
     lo(:)=tileDesc%limits(LOW,:)
     hi(:)=tileDesc%limits(HIGH,:)
     loHalo(:)=tileDesc%grownLimits(LOW,:)
     hihalo(:)=tileDesc%grownLimits(HIGH,:)

     thisBlock = thisBlock + 1

     ! get dimensions/limits and coordinates
     !! AD: Get Jared's help on this one

     iSize = hi(IAXIS)-lo(IAXIS)+1
     jSize = hi(JAXIS)-lo(JAXIS)+1
     kSize = hi(KAXIS)-lo(KAXIS)+1

     allocate(shock(loHalo(IAXIS):hiHalo(IAXIS),&
                    loHalo(JAXIS):hiHalo(JAXIS),&
                    loHalo(KAXIS):hiHalo(KAXIS)))

     ! identify the range of batches in each block (use floor/ceil in case of overlap)
     batch_lo(thisBlock) = nzones / xnet_nzbatchmx + 1
     nzones = nzones + iSize * jSize * kSize
     batch_hi(thisBlock) = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

     ! allocate space for dimensions
     allocate(xCoord(loHalo(IAXIS):hiHalo(IAXIS)))
     allocate(yCoord(loHalo(JAXIS):hiHalo(JAXIS)))
     allocate(zCoord(loHalo(KAXIS):hiHalo(KAXIS)))

     call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level,loHalo,hiHalo,xCoord)
     call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level,loHalo,hiHalo,yCoord)
     call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level,loHalo,hiHalo,zCoord)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Shock detector
     if (.NOT. bn_useShockBurn) then
        call Hydro_shockStrength(solnData, shock, lo, hi, loHalo, hiHalo, &
           xCoord,yCoord,zCoord,shock_thresh,shock_mode)
     else
        shock(:,:,:) = 0.0
     endif

     solnData(NMPI_VAR,:,:,:) = xnet_myid

     !$omp parallel do &
     !$omp collapse(3) &
     !$omp default(shared) &
     !$omp private(k,kk,j,jj,i,ii,okBurnTemp,okBurnDens,okBurnShock,okBurnNickel)
     do k = lo(KAXIS), hi(KAXIS)
        do j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)
              kk = k - lo(KAXIS) + 1
              jj = j - lo(JAXIS) + 1
              ii = i - lo(IAXIS) + 1

              tmp(ii,jj,kk,thisBlock)  = solnData(TEMP_VAR,i,j,k)
              rho(ii,jj,kk,thisBlock)  = solnData(DENS_VAR,i,j,k)
              sdot(ii,jj,kk,thisBlock) = 0.0e0

              zone(1,ii,jj,kk,thisBlock) = ii
              zone(2,ii,jj,kk,thisBlock) = jj
              zone(3,ii,jj,kk,thisBlock) = kk
              zone(4,ii,jj,kk,thisBlock) = thisBlock

              ! Map the solution data into the order required by bn_burner
              xIn(1:NSPECIES,ii,jj,kk,thisBlock) = solnData(xnet_inuc2unk,i,j,k)

              okBurnTemp = .FALSE.
              okBurnDens = .FALSE.
              okBurnShock = .FALSE.
              okBurnNickel = .FALSE.

              okBurnTemp = (tmp(ii,jj,kk,thisBlock) >= bn_nuclearTempMin .AND. tmp(ii,jj,kk,thisBlock) <= bn_nuclearTempMax)
              okBurnDens = (rho(ii,jj,kk,thisBlock) >= bn_nuclearDensMin .AND. rho(ii,jj,kk,thisBlock) <= bn_nuclearDensMax)
              okBurnShock = (shock(i,j,k) <= 0.0 .OR. (shock(i,j,k) > 0.0 .AND. bn_useShockBurn))

              if (okBurnTemp .AND. okBurnDens .AND. okBurnShock) then

                 if (NI56_SPEC /= NONEXISTENT) then
                    okBurnNickel = (solnData(NI56_SPEC,i,j,k) <  bn_nuclearNI56Max)
                 else    ! nickel is not even a species in this simulation, so we'll always burn
                    okBurnNickel = .TRUE.
                 endif

                 if (okBurnNickel) then
                    burnedZone(ii,jj,kk,thisBlock) = .TRUE.
                 endif

              endif

           end do
        end do
     end do
     !$omp end parallel do

     call tileDesc%releaseDataPtr(solnData, CENTER)
     nullify(solnData)
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
     deallocate(shock)

     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)

  call Timers_stop("burn_top")

  call Timers_start("burn_middle")

  ! get number of batches needed for all local zones (round up)
  batchCount = (nzones + xnet_nzbatchmx - 1) / xnet_nzbatchmx

  ! reshape all local zone data arrays into batches
  tmp_batch (1:xnet_nzbatchmx,1:batchCount) => tmp (:,:,:,:)
  rho_batch (1:xnet_nzbatchmx,1:batchCount) => rho (:,:,:,:)
  sdot_batch(1:xnet_nzbatchmx,1:batchCount) => sdot(:,:,:,:)
  xIn_batch (1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xIn (:,:,:,:,:)
  xOut_batch(1:NSPECIES,1:xnet_nzbatchmx,1:batchCount) => xOut(:,:,:,:,:)
  burnedZone_batch(1:xnet_nzbatchmx,1:batchCount) => burnedZone(:,:,:,:)
  zone_batch(1:4,1:xnet_nzbatchmx,1:batchCount) => zone(:,:,:,:,:)

  allocate(sumBurn_TS_batch(batchCount))

  !$omp parallel do &
  !$omp schedule(runtime) &
  !$omp default(shared)
  do m = 1, batchCount
     ! Do the actual burn
     call bn_burner(dt, tmp_batch(:,m), rho_batch(:,m), xIn_batch(:,:,m), &
          xOut_batch(:,:,m), sdot_batch(:,m), burnedZone_batch(:,m), zone_batch(:,:,m), sumBurn_TS_batch(m))
  end do
  !$omp end parallel do

  ! get average sumBurn_TS_batch per block
  sumBurn_TS = 0
  do thisBlock = 1, blockCount
     sumBurn_TS(thisBlock) = sum( sumBurn_TS_batch(batch_lo(thisBlock):batch_hi(thisBlock)) ) &
                           / ( batch_hi(thisBlock) - batch_lo(thisBlock) + 1 )
  end do

  deallocate(sumBurn_TS_batch)

  call Timers_stop("burn_middle")

  call Timers_start("burn_bottom")

  thisBlock = 0
  call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     thisBlock = thisBlock + 1

     ! get dimensions/limits and coordinates
     lo(:)=tileDesc%limits(LOW,:)
     hi(:)=tileDesc%limits(HIGH,:)

     ! Get a pointer to solution data
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Now put updated local data arrays back into unk through solnData pointer
     !$omp parallel do &
     !$omp collapse(3) &
     !$omp default(shared) &
     !$omp private(k,kk,j,jj,i,ii,ei,ek,enuc)
     do k = lo(KAXIS), hi(KAXIS)
        do j = lo(JAXIS), hi(JAXIS)
           do i = lo(IAXIS), hi(IAXIS)
              kk = k - lo(KAXIS) + 1
              jj = j - lo(JAXIS) + 1
              ii = i - lo(IAXIS) + 1

              ! Map the solution data into the order required by bn_burner
              solnData(xnet_inuc2unk,i,j,k) = xOut(1:NSPECIES,ii,jj,kk,thisBlock)

              !  NOTE should probably do something here with eintSwitch for consistency
              !  LBR will settle for simply using internal energy!
              ! kinetic energy
              ek = 0.5e0*(solnData(VELX_VAR,i,j,k)**2 +  &
                 solnData(VELY_VAR,i,j,k)**2 +  &
                 solnData(VELZ_VAR,i,j,k)**2)

              ! internal energy, add on nuclear rate*timestep
              enuc = dt*sdot(ii,jj,kk,thisBlock)
              ei = solnData(ENER_VAR,i,j,k) + enuc - ek

              solnData(ENUC_VAR,i,j,k) = sdot(ii,jj,kk,thisBlock)

              ! only update internal energy if the zone actually burned
              if (burnedZone(ii,jj,kk,thisBlock)) then
#ifdef EINT_VAR
                 solnData(EINT_VAR,i,j,k) = ei
#endif
                 solnData(ENER_VAR,i,j,k) = ei + ek
#ifdef EELE_VAR
                 solnData(EELE_VAR,i,j,k) = solnData(EELE_VAR,i,j,k) + enuc
#endif

                 xmass = xOut(1:NSPECIES,ii,jj,kk,thisBlock)
                 call bn_azbar

#ifdef YE_MSCALAR
                 solnData(YE_MSCALAR,i,j,k) = bye
#endif
#ifdef SUMY_MSCALAR
                 solnData(SUMY_MSCALAR,i,j,k) = ytot1
#endif
              end if
           end do
        end do
     end do
     !$omp end parallel do

     call Grid_setWork(tileDesc,real(sumBurn_TS(thisBlock)))
     solnData(MTSB_VAR,:,:,:) = sumBurn_TS(thisBlock)

     ! we've altered the EI, let's equilabrate
     if (any(burnedZone(:,:,:,thisBlock))) then

        call Eos_multiDim(MODE_DENS_EI,tileDesc%limits,solnData)

     end if

     call tileDesc%releaseDataPtr(solnData, CENTER)
     nullify(solnData)

     call itor%next()

  end do
  call Grid_releaseTileIterator(itor)

  call Timers_stop("burn_bottom")

  deallocate(xIn)
  deallocate(xOut)
  deallocate(sdot)
  deallocate(tmp)
  deallocate(rho)
  deallocate(burnedZone)
  deallocate(zone)
  deallocate(batch_lo)
  deallocate(batch_hi)
  deallocate(sumBurn_TS)

  call Timers_stop("burn")

  return
end subroutine Burn
