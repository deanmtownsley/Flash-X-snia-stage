!!****if* source/Simulation/SimulationMain/Cellular/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer, INTENT(in)::blockId,
!!                       integer, INTENT(in)::myPE)
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a planar detonation,
!!   perturbed with some noise.  This detonation should be unstable
!!   to becoming cellular at the detonation front.
!!
!! ARGUMENTS
!!
!!  blockId:        integer  the number of the block to initialize
!!  myPE   :        local processor number
!!
!!
!! PARAMETERS
!!
!!    xhe4               mass fraction of he4
!!    xc12               mass fraction of c12
!!    xn14               mass fraction of n14
!!    xo16               mass fraction of o16
!!    rhoAmbient         density of the cold upstream material 
!!    tempAmbient        temperature of the cold upstream material
!!    velxAmbient        x-velocity of the cold upstream material
!!    rhoPerturb         density of the post shock material
!!    tempPerturb        temperature of the post shock material
!!    velxPerturb        x-velocity of the post shock material
!!    radiusPerturb      distance below which the perturbation is applied
!!    xCenterPerturb     origin of the of the perturbation
!!    yCenterPerturb     origin of the of the perturbation
!!    zCenterPerturb     origin of the of the perturbation
!!    usePseudo1d        .true. for a 1d initial configuration, with the ??
!!                          copied along the y and z directions
!!                       .false. for a spherical configuration
!!    noiseAmplitude     amplitude of the white noise added to the perturbation
!!    noiseDistance      distances above and below radiusPerturb get noise added
!!    xmax               boundary of domain
!!    xmin               boundary of domain
!!    ymax               boundary of domain
!!    ymin               boundary of domain
!!    zmax               boundary of domain
!!    zmin               boundary of domain
!!
!!    smallx             smallest allowed abundance
!!    smlrho             smallest allowed density
!!
!!  NOTES
!!    Species used in this simulation are HE4 (helium-4), C12 (carbon-12), O16 (oxygen-16)
!!
!! NOTES
!!   See paper: Timmes, FX; Zingale, M; Olson, K; Fryxell, B; The Astrophysical
!!               Journal, Nov. 10, 2000 : 543: 938-954
!!  
!!***

subroutine Simulation_initBlock(solnData,tileDesc)

  use Simulation_data, ONLY: sim_smallRho, sim_smallx, sim_radiusPerturb, sim_usePseudo1d, &
     sim_xhe4, sim_xc12, sim_xn14, sim_xo16, &
     sim_rhoAmbient, sim_tempAmbient, sim_velxAmbient, &
     sim_rhoPerturb, sim_tempPerturb, sim_velxPerturb, &
     sim_noiseAmplitude, sim_noiseDistance, &
     sim_xCenterPerturb, sim_yCenterPerturb, sim_zCenterPerturb, &
     sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax, &
     sim_ign_keep_pres
  use Driver_interface, ONLY : Driver_abort
  use sim_local_interface, ONLY : sim_find_ign_state_keeping_pres
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
    Multispecies_getSumFrac
  use Grid_interface, ONLY : Grid_getCellCoords
  use Eos_interface, ONLY : Eos
  use Grid_tile, ONLY : Grid_tile_t
  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  !integer,INTENT(in) ::  blockId, myPE

!  Local variables
  real,              pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc
  
  integer,dimension(LOW:HIGH,MDIM) :: tileLimits
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits
  !real :: deltas(1:MDIM)

  real :: abar, zbar                 ! something to do with sum of mass fractions
  real :: xx, yy, zz, dist
  logical, parameter :: useGuardCell = .TRUE.

  !integer, dimension(2,MDIM), save               :: blockRange,blockExtent

  real, dimension(SPECIES_BEGIN:SPECIES_END) ::  massFraction
  real, dimension(1,NSPECIES) :: massFraction_buffer

  real,allocatable,dimension(:) :: xCoordsCell,yCoordsCell,zCoordsCell
  !integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  !integer :: sizeX,sizeY,sizeZ

  real, dimension(MDIM) :: deltas



  integer :: i, j, k, n
  integer, dimension(MDIM) :: iPosition   !for putting data with Grid_putData

  integer :: icount
  integer, parameter :: ifail = -1
  real, allocatable  :: rvec(:)                   ! for the random number generator
  integer            :: rvecSize=0             ! number of random numbers needed,
                                                     ! calculated below
  integer, parameter :: iseed = -867690
  integer            :: iseedUse

  ! variables needed for the eos call
  real :: temp_zone, rho_zone, vel_zone
  real :: ptot, eint, etot, gamma
  real, dimension(EOS_NUM)  :: eosData
  real, dimension(1,EOS_VARS) :: eosData_buffer


  real :: left, right, tleft, tright

  iseedUse = iseed
! ----------------------------------------------------------------------------------------------

 ! determine size of blocks
  !call Grid_getBlkIndexLimits(blockId,blockRange,blockExtent)
  tileLimits = tileDesc%limits
  grownTileLimits = tileDesc%grownLimits
    
  call tileDesc%deltas(deltas)

  ! Get the indices of the blocks
  !call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  !sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoordsCell(grownTileLimits(LOW, IAXIS):grownTileLimits(HIGH, IAXIS)))
  !sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoordsCell(grownTileLimits(LOW, JAXIS):grownTileLimits(HIGH, JAXIS)))
  !sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoordsCell(grownTileLimits(LOW, KAXIS):grownTileLimits(HIGH, KAXIS)))

  !call Grid_getDeltas(blockId,deltas)

!  if (NDIM == 3)  &
   !call Grid_getCellCoords(KAXIS, blockId, CENTER, useGuardCell, zCoordsCell, sizeZ)
   call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, & 
                  grownTileLimits(LOW,:), grownTileLimits(HIGH,  :), zCoordsCell)

!  if (NDIM >= 2)  &
   call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                  grownTileLimits(LOW,:), grownTileLimits(HIGH,  :), yCoordsCell)

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                  grownTileLimits(LOW,:), grownTileLimits(HIGH,  :), xCoordsCell)

   !call Grid_getCellCoords(JAXIS, blockId, CENTER, useGuardCell, yCoordsCell, sizeY)
   !call Grid_getCellCoords(IAXIS, blockId, CENTER, useGuardCell, xCoordsCell, sizeX)

  !! NOTE the incorrect syntax below --  this causes a crash for KAXIS when NDIM < 3, 
  !!    works OK if you substitute 1 at the MAXCELLS location
!  call Grid_getCellCoords(KAXIS, blockID,CENTER,useGuardCell,zCoordsCell,1)
!  call Grid_getCellCoords(JAXIS, blockID,CENTER,useGuardCell,yCoordsCell,MAXCELLS)
!  call Grid_getCellCoords(IAXIS, blockID,CENTER,useGuardCell,xCoordsCell,MAXCELLS)

  ! the initial composition
  massFraction(:)    = sim_smallx 
  if (HE4_SPEC > 0) massFraction(HE4_SPEC) = max(sim_xhe4,sim_smallx)
  if (C12_SPEC > 0) massFraction(C12_SPEC) = max(sim_xc12,sim_smallx)
#ifdef N14_SPEC
  if (N14_SPEC > 0) massFraction(N14_SPEC) = max(sim_xn14,sim_smallx)
#else
  if (sim_xn14 > 0.0) then
    call Driver_abort("Nonzero N14 fraction used in sim without N14")
  endif
#endif
  if (O16_SPEC > 0) massFraction(O16_SPEC) = max(sim_xo16,sim_smallx)

  call Multispecies_getSumInv(A,abar,massFraction)
  abar = 1.0 / abar
  call Multispecies_getSumFrac(Z,zbar,massFraction)
  zbar = abar * zbar


  !..get a blocks worth of random numbers between 0.0 and 1.0
  rvecSize = grownTileLimits(HIGH,KAXIS)*grownTileLimits(HIGH,JAXIS)*grownTileLimits(HIGH,IAXIS)
  allocate(rvec(rvecSize))
  call sim_ranmar(iseedUse, rvec, rvecSize)

  icount = 0

  ! now fill the master arrays

  do k = grownTileLimits(LOW,KAXIS), grownTileLimits(HIGH,KAXIS)
     if (NDIM == 3) then
        iPosition(3) = k
        zz = zCoordsCell(k)
     endif

     do j = grownTileLimits(LOW,JAXIS), grownTileLimits(HIGH,JAXIS)
        if (NDIM >= 2) then
           iPosition(2) = j
           yy = yCoordsCell(j)
        endif

        do i = grownTileLimits(LOW,IAXIS), grownTileLimits(HIGH,IAXIS)
           iPosition(1) = i
           xx = xCoordsCell(i)
           icount = icount + 1


           ! compute the distance from the center
           if (NDIM .EQ. 1) then
              dist = xx - sim_xCenterPerturb
           else if (NDIM .EQ. 2) then
              if (sim_usePseudo1d) then
                 dist = xx - sim_xCenterPerturb
              else
                 dist = sqrt((xx - sim_xCenterPerturb)**2 + & 
                      &                 (yy - sim_yCenterPerturb)**2)
              endif
           elseif (NDIM .EQ. 3) then
              if (sim_usePseudo1d) then
                 dist = xx - sim_xCenterPerturb
              else
                 dist = sqrt((xx - sim_xCenterPerturb)**2 + & 
                      &                 (yy - sim_yCenterPerturb)**2 + & 
                      &                 (zz - sim_zCenterPerturb)**2)
              endif
           endif


           rho_zone = sim_rhoAmbient
           vel_zone = sim_velxAmbient
           ! calculate temperature by averaging the gradient profile over the cell
           ! cell edges
           left = max(0.0,dist-0.5*deltas(IAXIS))
           right = dist+0.5*deltas(IAXIS)
           if (right <= sim_radiusPerturb) then
              ! cell fully in gradient region
              tleft = sim_tempPerturb - &
                           (sim_tempPerturb-sim_tempAmbient)*left/sim_radiusPerturb
              tright = sim_tempPerturb - &
                           (sim_tempPerturb-sim_tempAmbient)*right/sim_radiusPerturb
              temp_zone = 0.5*(tleft+tright)
              if ( sim_ign_keep_pres ) then
                 eosData(EOS_TEMP) = sim_tempAmbient
                 eosData(EOS_DENS) = sim_rhoAmbient
                 
                 call sim_find_ign_state_keeping_pres(temp_zone, eosData,massFraction)
                 
                 temp_zone = eosData(EOS_TEMP)
                 rho_zone = eosData(EOS_DENS)
              endif
           else if (left <= sim_radiusPerturb) then
              ! cell straddles edge of gradient region
              tleft = sim_tempPerturb - &
                           (sim_tempPerturb-sim_tempAmbient)*left/sim_radiusPerturb
              temp_zone = 0.5*(tleft+sim_tempAmbient)*(sim_radiusPerturb-left)/(right-left) &
                          + sim_tempAmbient*(right-sim_radiusPerturb)/(right-left)
              if ( sim_ign_keep_pres ) then
                 eosData(EOS_TEMP) = sim_tempAmbient
                 eosData(EOS_DENS) = sim_rhoAmbient

                 call sim_find_ign_state_keeping_pres(temp_zone,eosData,massFraction)
                 temp_zone = eosData(EOS_TEMP)
                 rho_zone = eosData(EOS_DENS)
              endif

           else
              temp_zone=sim_tempAmbient
           endif
              

           !..seed the initial conditions with some white noise
           if ( abs(dist - sim_radiusPerturb) .le. sim_noiseDistance) then
            !!  print *, 'sim_noiseAmplitude = ',sim_noiseAmplitude
            !!  print *,'rvec, icount', icount, rvec(icount)
              rho_zone = rho_zone *  & 
                   &               (1.0 + sim_noiseAmplitude * (1.0 - 2.0 * rvec(icount)))
           end if


           !  Need input of density and temperature
           eosData(EOS_TEMP) = temp_zone
           eosData(EOS_DENS) = rho_zone

           eosData_buffer(1,:) = eosData(1:EOS_VARS)
           massFraction_buffer(1,:) = massFraction(SPECIES_BEGIN: SPECIES_END)
           call Eos_vector(MODE_DENS_TEMP,1,eosData_buffer,massFraction_buffer)
           massFraction(SPECIES_BEGIN: SPECIES_END) = massFraction_buffer(1,:)
           eosData(1:EOS_VARS) = eosData_buffer(1,:)

           temp_zone = eosData(EOS_TEMP)
           rho_zone = eosData(EOS_DENS)
           ptot = eosData(EOS_PRES)
           eint = eosData(EOS_EINT)
           gamma = eosData(EOS_GAMC)

           ! calculate kinetic energy and total energy
           !! this was NOT done in flash2
           etot = eint + 0.5*vel_zone**2

           ! store the values
           ! fill the flash arrays

           !call Grid_putPointData(blockId,CENTER,TEMP_VAR,EXTERIOR,iPosition,temp_zone)
           solnData(DENS_VAR, i, j, k) = rho_zone
           solnData(TEMP_VAR, i, j, k) = temp_zone
           solnData(PRES_VAR, i, j, k) = ptot

           solnData(PRES_VAR, i, j, k) = ptot
           solnData(EINT_VAR, i, j, k) = eint
           solnData(PRES_VAR, i, j, k) = etot
           solnData(PRES_VAR, i, j, k) = gamma
           solnData(PRES_VAR, i, j, k) = (ptot/(etot*sim_rhoAmbient) + 1.0)
           solnData(PRES_VAR, i, j, k) = vel_zone
           solnData(PRES_VAR, i, j, k) = 0.0e0
           solnData(PRES_VAR, i, j, k) = 0.0e0
           !call Grid_putPointData(blockId,CENTER,DENS_VAR,EXTERIOR,iPosition,rho_zone)
           !call Grid_putPointData(blockId,CENTER,PRES_VAR,EXTERIOR,iPosition,ptot)
           !call Grid_putPointData(blockId,CENTER,EINT_VAR,EXTERIOR,iPosition,eint)
           !call Grid_putPointData(blockId,CENTER,ENER_VAR,EXTERIOR,iPosition,etot)
           !call Grid_putPointData(blockId,CENTER,GAMC_VAR,EXTERIOR,iPosition,gamma)
           !call Grid_putPointData(blockId,CENTER,GAME_VAR,EXTERIOR,iPosition,(ptot/(etot*sim_rhoAmbient) + 1.0))
           !call Grid_putPointData(blockId,CENTER,VELX_VAR,EXTERIOR,iPosition,vel_zone)
           !! No need to do dimensional check here, as the VELZ and VELY are always defined in the config file
           !call Grid_putPointData(blockId,CENTER,VELY_VAR,EXTERIOR,iPosition,0.0e0)
           !call Grid_putPointData(blockId,CENTER,VELZ_VAR,EXTERIOR,iPosition,0.0e0)

           do n = SPECIES_BEGIN,SPECIES_END
              !call Grid_putPointData(blockID,CENTER,n,EXTERIOR,iPosition,massFraction(n))
              solnData(n,i,j,k) = massFraction(n)
           enddo

           !..end of 3d loops
        enddo  ! end of k loop
     enddo     ! end of j loop
  enddo        ! end of i loop

  ! cleanup
  deallocate(rvec)
  deallocate(xCoordsCell)
  deallocate(yCoordsCell)
  deallocate(zCoordsCell)



  return
end subroutine Simulation_initBlock


