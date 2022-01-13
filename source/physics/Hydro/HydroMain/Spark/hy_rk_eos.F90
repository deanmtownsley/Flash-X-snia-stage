!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_eos
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!  NAME
!!
!!  hy_rk_eos
!!
!!  SYNOPSIS
!!
!!  call hy_rk_eos ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_eos(limits)

  use Hydro_data, ONLY : hy_starState, hy_threadWithinBlock,hy_eosData, hy_mfrac
  use Eos_interface, ONLY : Eos_putData, Eos_getData, Eos

  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
  integer :: i,j,k,vecLen
  integer,dimension(MDIM) :: pos, datasize
  integer, dimension(LOW:HIGH,MDIM)  :: range
  !Note this allocation works around the inability of FLASH5 to produce MAXCELLS constant in 
  !Simulation.h using preprocessors alone.  This may want to be changed for future.
  !real, dimension(NSPECIES*MAXCELLS) :: hy_mfrac
  !real, dimension(EOS_NUM*MAXCELLS) :: hy_eosData
  real, pointer :: tempData(:,:,:,:)

  datasize = limits(HIGH,:)-limits(LOW,:)+1 


  tempData   => hy_starState
  !$omp parallel if (hy_threadWithinBlock .AND. NDIM > 1) &
  !$omp default(none) &
  !$omp firstprivate(vecLen) &
  !$omp private(i,j,k,pos,hy_eosData,hy_mfrac,range)&
  !$omp shared(limits,hy_starState,tempData)

  pos(IAXIS) = limits(LOW,IAXIS)
  vecLen     = limits(HIGH,IAXIS)-limits(LOW,IAXIS)+1

  !  Begin loop over zones
  !$omp do schedule(static) collapse(2)
  do k = limits(LOW,KAXIS), limits(HIGH,KAXIS)
     do j = limits(LOW,JAXIS), limits(HIGH,JAXIS)

       !pos(JAXIS) = j
       !pos(KAXIS) = k

       !Here we are getting/putting data by vectors along the X direction
       !X Dim
       range(LOW,1) = limits(LOW,IAXIS)
       range(HIGH,1) = limits(HIGH,IAXIS)
       !Y dim
       range(LOW,2) = j
       range(HIGH,2) = j
       !Z dim
       range(LOW,3) = k
       range(HIGH,3) = k
!!!*** hand this part over to the composer
       call Eos_getData(range,vecLen,tempData,CENTER,hy_eosData,hy_mfrac)
       call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
       call Eos_putData(range,vecLen,tempData,CENTER,hy_eosData)
     end do
  end do
  !$omp end do nowait
  !$omp end parallel

  nullify(tempData)
end subroutine hy_rk_eos





!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_eos_offloaded
!!
!!  NAME
!!
!!  hy_rk_eos_offloaded
!!
!!  SYNOPSIS
!!
!!  call hy_rk_eos_offloaded ( integer(IN) :: blockID )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine hy_rk_eos_offloaded(limits)

  use Hydro_data, ONLY : hy_threadWithinBlock
  ! use Eos_interface, ONLY : Eos_putData, Eos_getData, Eos
  use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_mapLookup, eos_meshMe
  use Driver_interface, ONLY : Driver_abortFlash
#ifdef OMP_OL
  use eos_idealGammaData, ONLY: eos_gammam1
#endif
  use Eos_data, ONLY : eos_gasConstant, eos_gamma, eos_singleSpeciesA, eos_singleSpeciesZ
  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
  integer :: i,j,k,vecLen

  if(CENTER == SCRATCH) then
    call Driver_abortFlash("Eos_getData : the use of SCRATCH is deprecated")
  end if

  vecLen     = limits(HIGH,IAXIS)-limits(LOW,IAXIS)+1

!!$!   ! !  Begin loop over zones
!!$!   ! !$omp do schedule(static) collapse(2)
!!$!  ! !$omp target map(to:limits,i,k,range,vecLen,eos_eintSwitch,eos_smalle,eos_mapLookup,eos_gammam1,eos_gasConstant,eos_gamma,eos_singleSpeciesA,eos_singleSpeciesZ,eos_meshMe,hy_eosData_size) map(alloc:hy_eosData)
!!$!   !$omp parallel if (hy_threadWithinBlock .AND. NDIM > 1) &
!!$!   !$omp default(none) &
!!$!   !$omp private(i,j,k)&
!!$!   !$omp shared(limits,hy_starState, eos_gammam1, eos_gasConstant, eos_gamma, eos_singleSpeciesA, eos_singleSpeciesZ, eos_mapLookup,eos_smalle,eos_eintSwitch) &
!!$!   !$omp do collapse(3)
#ifdef OMP_OL
  !$omp target teams distribute parallel do &
  !$omp shared(limits,eos_gammam1, eos_gasConstant, eos_gamma, eos_singleSpeciesA, &
  !$omp eos_singleSpeciesZ, eos_mapLookup,eos_smalle,eos_eintSwitch,vecLen) &
  !$omp map(to:limits,eos_gammam1, eos_gasConstant, eos_gamma, eos_singleSpeciesA, &
  !$omp eos_singleSpeciesZ, eos_mapLookup,eos_smalle,eos_eintSwitch,vecLen) &
  !$omp private(i,j,k) collapse(3) default(none)
  do k = limits(LOW,KAXIS), limits(HIGH,KAXIS)
     do j = limits(LOW,JAXIS), limits(HIGH,JAXIS)
        do i = limits(LOW,IAXIS), limits(HIGH,IAXIS)
          call eos_idealGamma_offloaded(vecLen,i,j,k, eos_gammam1, &
               eos_gasConstant, eos_gamma, eos_singleSpeciesA, eos_singleSpeciesZ, &
               eos_mapLookup,eos_smalle,eos_eintSwitch)
        end do
     end do
  end do
#endif

end subroutine hy_rk_eos_offloaded



subroutine eos_idealGamma_offloaded(vecLen,i, j, k, eos_gammam1, &
     eos_gasConstant, eos_gamma, eos_singleSpeciesA, eos_singleSpeciesZ, &
     eos_mapLookup,eos_smalle,eos_eintSwitch)

  !==============================================================================
    ! use Eos_data, ONLY : eos_gasConstant, eos_gamma, &
    !      eos_singleSpeciesA, eos_singleSpeciesZ
    ! use eos_idealGammaData, ONLY: eos_gammam1
    ! use Driver_interface, ONLY : Driver_abortFlash
  use Hydro_data, only : hy_starState

  implicit none
  
#include "Eos.h"
#include "constants.h"
#include "Simulation.h"
#include "Eos_map.h"

  !     Arguments
  real, intent(IN) :: eos_gammam1, eos_gasConstant, eos_gamma, eos_singleSpeciesA, eos_singleSpeciesZ
  integer, dimension(1:EOSMAP_NUM_ROLES, 1:2, 1:5), intent(in) :: eos_mapLookup
  integer, intent(in) :: i,j,k,vecLen
  real, intent(in) :: eos_smalle
  real, intent(in) :: eos_eintSwitch
  real ::  ggprod, ggprodinv, gam1inv
  integer :: dens, temp, pres, eint, abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, c_v, c_p, gamc, pel, ne, eta
  integer :: ilo,ihi, m, n
  integer :: pres_map,entr_map,gamc_map,temp_map
  integer :: eint_map,ener_map,game_map, dens_map
  integer :: velx_map, vely_map, velz_map, sumy_map, ye_map
  real :: kineticEnergy, internalEnergy
  real :: massFrac(SPECIES_END-SPECIES_BEGIN+1)
  
  !$omp declare target

  ggprod = eos_gammam1 * eos_gasConstant
  !============================================================================

  ! load massFrac
  m = 1
  do n = SPECIES_BEGIN,SPECIES_END
    massFrac(m) = hy_starState(n,i,j,k)
    m=m+1
  end do

  pres_map = eos_mapLookup(EOSMAP_PRES,EOS_IN,CENTER)
  dens_map = eos_mapLookup(EOSMAP_DENS,EOS_IN,CENTER)
  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_IN,CENTER)
  gamc_map = eos_mapLookup(EOSMAP_GAMC,EOS_IN,CENTER)
  game_map = eos_mapLookup(EOSMAP_GAME,EOS_IN,CENTER)
  eint_map = eos_mapLookup(EOSMAP_EINT,EOS_IN,CENTER)
  ener_map = eos_mapLookup(EOSMAP_ENER,EOS_IN,CENTER)
  velx_map = eos_mapLookup(EOSMAP_VELX,EOS_IN,CENTER)
  vely_map = eos_mapLookup(EOSMAP_VELY,EOS_IN,CENTER)
  velz_map = eos_mapLookup(EOSMAP_VELZ,EOS_IN,CENTER)
  sumy_map = eos_mapLookup(EOSMAP_SUMY,EOS_IN,CENTER)
  ye_map   = eos_mapLookup(EOSMAP_YE,  EOS_IN,CENTER)
  entr_map = eos_mapLookup(EOSMAP_ENTR,EOS_IN,CENTER)

  ! Get interal energy
  if (velx_map > 0 .AND. vely_map > 0 .AND. velz_map > 0) then
    kineticEnergy  = 0.5*(hy_starState(velx_map,i,j,k)**2 + &
                          hy_starState(vely_map,i,j,k)**2 + &
                          hy_starState(velz_map,i,j,k)**2)
  else
    kineticEnergy = 0.0
  end if

  if(eint_map /= NONEXISTENT) then
    internalEnergy  = hy_starState(eint_map,i,j,k)
    if(ener_map /= NONEXISTENT) then
        if ( hy_starState(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
          internalEnergy = hy_starState(ener_map,i,j,k) - kineticEnergy
        end if
    end if
  else if(game_map /= NONEXISTENT) then ! This case should be usable for R(elativistic)HD - KW
    internalEnergy  = hy_starState(pres_map,i,j,k) / hy_starState(dens_map,i,j,k) / &
                          (hy_starState(game_map,i,j,k) - 1.0)
    if(ener_map /= NONEXISTENT) then
        if ( hy_starState(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
          internalEnergy = hy_starState(ener_map,i,j,k) - kineticEnergy
        end if
    end if
  else if(ener_map /= NONEXISTENT) then
    internalEnergy = hy_starState(ener_map,i,j,k)-kineticEnergy
  else
    internalEnergy = eos_smalle
  endif
  
  ! save internal and kinetic energy
  internalEnergy = max(internalEnergy, eos_smalle)

  
  pres_map = eos_mapLookup(EOSMAP_PRES,EOS_OUT,CENTER)
  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_OUT,CENTER)
  gamc_map = eos_mapLookup(EOSMAP_GAMC,EOS_OUT,CENTER)
  game_map = eos_mapLookup(EOSMAP_GAME,EOS_OUT,CENTER)
  eint_map = eos_mapLookup(EOSMAP_EINT,EOS_OUT,CENTER)
  ener_map = eos_mapLookup(EOSMAP_ENER,EOS_OUT,CENTER)
  entr_map = eos_mapLookup(EOSMAP_ENTR,EOS_OUT,CENTER)

  ! density, internal energy taken as input
  ggprodinv = 1. / ggprod
  gam1inv   = 1. / eos_gammam1
  hy_starState(pres_map,i,j,k) = hy_starState(dens_map,i,j,k) * &
                                internalEnergy * gam1inv
  hy_starState(temp_map,i,j,k) = internalEnergy * ggprodinv * &
                                eos_singleSpeciesA
  if(entr_map /= NONEXISTENT) then
    hy_starState(entr_map,i,j,k) = (hy_starState(pres_map,i,j,k)/hy_starState(dens_map,i,j,k) +  &
          &  internalEnergy)/hy_starState(temp_map,i,j,k)
  end if
  
  if(ener_map /= NONEXISTENT)hy_starState(ener_map,i,j,k) = internalEnergy+kineticEnergy


  hy_starState(game_map,i,j,k) = hy_starState(pres_map,i,j,k)/&
    (internalEnergy *hy_starState(dens_map,i,j,k)) +1
  return
end subroutine eos_idealGamma_offloaded
  
