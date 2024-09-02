!!****f* source/physics/sourceTerms/Burn/Burn
!!
!! NAME
!!
!!  Burn
!!
!!
!! SYNOPSIS
!!
!!   call Burn ( integer, intent(IN)    :: blockCount, 
!!               integer(:), intent(IN) :: blockList, 
!!               real, intent(IN)       ::  dt  )    
!!
!! DESCRIPTION
!!
!!  Apply burner to all blocks in specified list.
!!
!! ARGUMENTS
!!
!!   blockCount -- dimension of blockList
!!   blockList -- array of blocks which should receive burning
!!   dt  --       passed to the internal bn_burner module  
!!
!! PARAMETERS
!!
!!  useBurn -- Boolean, True.  Turns on burning module
!!  useShockBurn -- Boolean, FALSE.  Controls whether burning is allowed inside
!!                a regime experiencing shocks
!!  enucDtFactor -- Real, 1.0E+30.  Timestep limiter.  See Burn_computeDt for details.               
!!
!! NOTES
!!
!!  The burning unit adds a new mesh variable ENUC_VAR which is the nuclear energy 
!!             generation rate
!!
!!***

subroutine Burn (dt)    
  use omp_lib
  use Burn_data, ONLY: bn_useBurn, bn_useShockBurn, bn_map_fi_to_mi, bn_limiter_max_dlnT
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Grid_interface, ONLY  : Grid_fillGuardCells, &
       Grid_getBlkIndexLimits, Grid_getCellCoords, &
       Grid_getSingleCellVol, Grid_getMinCellSize, &
       Grid_getTileIterator, Grid_releaseTileIterator
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t
  use Eos_interface, ONLY   : Eos_multiDim, Eos
  use Hydro_interface, ONLY : Hydro_detectShock
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Multispecies_interface, ONLY: Multispecies_getSumFrac

  use mesa_flash_nucreact, only: burn_const_dens_temp

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
#include "Multispecies.h"

  !args
  !integer, INTENT(in)                        :: blockCount
  !integer, INTENT(in), DIMENSION(blockCount)  :: blockList
  real,    INTENT(in)                        :: dt

  logical :: gcMask(NUNK_VARS)
  !integer :: blockiteri, blockID
  integer :: i,j,k, fi

  real, pointer, dimension(:,:,:,:)                    :: solnData

  !integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC

  real, allocatable, dimension(:)         :: xCoord, yCoord, zCoord
  integer,dimension(LOW:HIGH,MDIM) :: grownTileLimits
  real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC) :: shock

  real, dimension(NSPECIES):: fabund, mabund
  real, dimension(1,NSPECIES) :: fabund_buff
  real :: qbar_ini, qbar_fin, Ye_ini, Ye_fin, de
  real :: cgsMeVperGram = 9.6485e17
  real :: pe_n_mdiff = 938.2723+0.5110-939.5656
  real :: temp,rho, t_sc, dx_min, dt_nuc, de_limit
  real :: eosData(EOS_NUM)
  real, dimension(1,EOS_VARS) :: eosData_buff
  real, dimension(1,EOS_VARS+1,EOS_NUM) :: derivs
  logical :: eosMask(EOS_VARS+1:EOS_NUM)
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tileDesc
  


  ! --------------------------------------------------------
  ! 0. do nothing if burn is turned off and initalize some things
  if (.not. bn_useBurn) return

  ! maybe should live outside of this routine, but this is current convention
  call Timers_start("burn")

  ! initialize some things that are unchanged in various loops
  !eosMask(:) = .false.
  !eosMask(EOS_CV) = .true.
  !eosMask(EOS_DET) = .true.
  call Grid_getMinCellSize(dx_min)
  ! ---------------------------------------------------
  ! 1. bring guardcells up to date
  ! ---------------------------------------------------
  gcMask = .FALSE.
  ! If burning in shocks is disabled (the normal case)
  ! we will need the pressure and velocity in order to do shock detection
  if (.not. bn_useShockBurn) then
     gcMask(PRES_VAR) = .TRUE.
     gcMask(VELX_VAR) = .TRUE.
#if NDIM > 1
     gcMask(VELY_VAR) = .TRUE.
#endif
#if NDIM > 2
     gcMask(VELZ_VAR) = .TRUE.
#endif
     ! fill guardcells
     ! we don't need to do the EOS because if there is a shock at the
     ! edge of the block, then the next block should be fully refined also and
     ! thus at the same refinement as this one.  The EOS call is only
     ! really needed when the neighboring block is at a different
     ! refinement level.
     call Grid_fillGuardCells(CENTER, ALLDIR, &
       doEos=.false., maskSize=NUNK_VARS,  mask = gcMask, makemaskConsistent=.true.)
  endif

  ! -------------------------------------------------------
  ! 2.  iterate over blocks, apply burning evolution and update quantities
  ! -------------------------------------------------------
  call Grid_getTileIterator(itor, nodetype=LEAF)
  !do blockiteri = 1, blockCount
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     !blockID = blockList(blockiteri)
     call tileDesc%getDataPtr(solnData, CENTER)
     
     ! -------------------------------
     ! 2.1 initialize quantities for this block
     !     including shock and proximity, which require unchanged neighbor cells
     ! -------------------------------
     !call Grid_getBlkPtr(blockID,solnData)
     !call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)

     ! shock detect if burning is turned off in shocks
     if (.NOT. bn_useShockBurn) then
        ! get coordinate positions, used for shock detection
        grownTileLimits = tileDesc%grownLimits
        allocate(xCoord(grownTileLimits(LOW,IAXIS):grownTileLimits(HIGH,IAXIS)))
        allocate(yCoord(grownTileLimits(LOW,JAXIS):grownTileLimits(HIGH,JAXIS)))
        allocate(zCoord(grownTileLimits(LOW,KAXIS):grownTileLimits(HIGH,KAXIS)))
        call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          xCoord)
        call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          yCoord)
        call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          grownTileLimits(LOW,  :), &
                          grownTileLimits(HIGH, :), &
                          zCoord)

        call Hydro_detectShock(solnData, shock, tileDesc%limits, grownTileLimits, (/0,0,0/), &
             xCoord,yCoord,zCoord)

        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
     else
        shock(:,:,:) = 0
     endif

#ifdef SHK_VAR
     solnData(SHK_VAR,tileDesc%limits(LOW,IAXIS):tileDesc%limits(HIGH,IAXIS),   &
                      tileDesc%limits(LOW,JAXIS):tileDesc%limits(HIGH,JAXIS),   &
                      tileDesc%limits(LOW,KAXIS):,tileDesc%limits(HIGH,KAXIS)) = &
                shock(tileDesc%limits(LOW,IAXIS):tileDesc%limits(HIGH,IAXIS),   &
                      tileDesc%limits(LOW,JAXIS):tileDesc%limits(HIGH,JAXIS),   &
                      tileDesc%limits(LOW,KAXIS):tileDesc%limits(HIGH,KAXIS))
#endif

     ! --------------------------------
     ! 2.2 loop over all interior zones and apply evolution
     ! --------------------------------
#if NDIM==3
     !$omp parallel do private(i,j,k,temp,rho,fi,fabund,mabund,de,&
     !$omp qbar_ini,Ye_ini,qbar_fin,Ye_fin, t_sc, de_limit, eosData, dt_nuc)
#endif
     do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
#if NDIM==2
     !$omp parallel do private(i,j,temp,rho,fi,fabund,mabund,de,&
     !$omp qbar_ini,Ye_ini,qbar_fin,Ye_fin, t_sc, de_limit, eosData, dt_nuc)
#endif
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
#if NDIM==1
     !$omp parallel do private(i,temp,rho,fi,fabund,mabund,de &
     !$omp ,qbar_ini,Ye_ini,qbar_fin,Ye_fin, t_sc, de_limit, eosData, dt_nuc)
#endif
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
              ! skip this cell if shock burning is turned off (normally) and it is in a shock
              if ( (.not. bn_useShockBurn) .and. ( shock(i,j,k) == 1.0 ) ) then
                 solnData(ENUC_VAR,i,j,k) = 0.0
                 cycle
              endif
              temp = solnData(TEMP_VAR,i,j,k)
              rho  = solnData(DENS_VAR,i,j,k)
              fabund(1:NSPECIES) = solnData(SPECIES_BEGIN:SPECIES_END,i,j,k)

              ! need c_V to compute de_limit
              eosData(EOS_TEMP) = temp
              eosData(EOS_DENS) = rho
              fabund_buff(1,:) = fabund(:)
              eosData_buff(1,:) = eosData(1:EOS_VARS) 
              call Eos_vector(MODE_DENS_TEMP, 1, eosData_buff, fabund_buff, derivs=derivs)
              eosData(1:EOS_VARS) = eosData_buff(1.:)
              eosData(EOS_VARS+1:EOS_NUM) = derivs(:)

              t_sc = dx_min/sqrt(solnData(GAMC_VAR,i,j,k)*solnData(PRES_VAR,i,j,k)/rho)
              de_limit = bn_limiter_max_dlnT*temp*eosData(EOS_CV)/t_sc * dt

              call Multispecies_getSumFrac(EB, qbar_ini, fabund)
              call Multispecies_getSumFrac(Z, Ye_ini, fabund)

              ! reorder abundances to mesa order
              do fi=1,NSPECIES
                 mabund( bn_map_fi_to_mi(fi) ) = fabund(fi)
              end do
              dt_nuc = dt
              ! TODO should get neutrino losses
              ! call the mesa wrapper to integrate the network
              call burn_const_dens_temp( dt_nuc, de_limit, solnData(DENS_VAR,i,j,k), solnData(TEMP_VAR,i,j,k), mabund )
              ! restuff
              do fi=1,NSPECIES
                 fabund(fi) = mabund( bn_map_fi_to_mi(fi) )
              end do

              call Multispecies_getSumFrac(EB, qbar_fin, fabund)
              call Multispecies_getSumFrac(Z, Ye_fin, fabund)

              solnData(SPECIES_BEGIN:SPECIES_END,i,j,k) = fabund(1:NSPECIES)

              de = ( qbar_fin-qbar_ini + (Ye_fin-Ye_ini)*(pe_n_mdiff) )*cgsMeVperGram

              solnData(ENER_VAR,i,j,k) = solnData(ENER_VAR,i,j,k) + de
#ifdef EINT_VAR
              solnData(EINT_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + de
#endif
              solnData(ENUC_VAR,i,j,k) = de/dt

           enddo
#if NDIM==1
    !$omp end parallel do
#endif
        enddo
#if NDIM==2
    !$omp end parallel do
#endif
     enddo
#if NDIM==3
    !$omp end parallel do
#endif

     ! --------------------------------
     ! 2.3 Finish up:
     !     Update interior EOS quantities for this block and release stuff
     ! --------------------------------
     call Eos_multiDim(MODE_DENS_EI,tileDesc%limits,solnData)

     !call Grid_releaseBlkPtr(blockID,solnData)
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()

  ! end of iteration over blockslist passed in
  end do
  call Grid_releaseTileIterator(itor)
  call Timers_stop("burn")

  return

end subroutine Burn
