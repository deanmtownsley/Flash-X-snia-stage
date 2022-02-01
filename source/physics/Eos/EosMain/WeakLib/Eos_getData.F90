!!****if* source/physics/Eos/EosMain/WeakLib/Eos_getData
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
!!  Eos_getData
!! 
!!  For more details see the documentation of the NULL implementation
!!
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_getData &
  (range,vecLen,solnData,gridDataStruct,eosData,massFrac, eosMask)

  use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_mapLookup
  use Driver_interface, ONLY : Driver_abort
  implicit none
  
#include "Eos.h"
#include "Eos_map.h"
#include "constants.h"
#include "Simulation.h"
 
  integer, intent(in) :: vecLen, gridDataStruct
  integer, dimension(LOW:HIGH,MDIM), intent(in) :: range
  real, dimension(EOS_NUM*vecLen),intent(INOUT) :: eosData
  real,dimension(:),optional,intent(OUT) :: massFrac
  logical, optional, INTENT(INOUT),dimension(EOS_VARS+1:) :: eosMask
  real, pointer:: solnData(:,:,:,:)


  integer :: i,j,k,n,m,pres,dens,gamc,temp,abar,zbar,eint,ekin,entr,elef
  integer :: pres_map,dens_map,gamc_map,game_map,temp_map,entr_map
  integer :: eint_map,ener_map, velx_map, vely_map, velz_map, sumy_map, ye_map
  integer :: ib,ie,jb,je,kb,ke
  real :: kineticEnergy, internalEnergy

!! ---------------------------------------------------------------------------------
  ! Test calling arguments

  ! Initializations:   grab the solution data from UNK and determine
  !   the length of the data being operated upon

  ib=range(LOW,IAXIS)
  jb=range(LOW,JAXIS)
  kb=range(LOW,KAXIS)
  ie=range(HIGH,IAXIS)
  je=range(HIGH,JAXIS)
  ke=range(HIGH,KAXIS)
  
!!$  select case(axis)
!!$  case(IAXIS)
!!$     ie=ie+vecLen-1
!!$  case(JAXIS)
!!$     je=je+vecLen-1
!!$  case(KAXIS)
!!$     ke=ke+vecLen-1
!!$  end select

  ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  ekin = (EOS_EKIN-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen
  elef = (EOS_YE-1)*vecLen

  pres_map = eos_mapLookup(EOSMAP_PRES,EOS_IN,gridDataStruct)
  dens_map = eos_mapLookup(EOSMAP_DENS,EOS_IN,gridDataStruct)
  temp_map = eos_mapLookup(EOSMAP_TEMP,EOS_IN,gridDataStruct)
  gamc_map = eos_mapLookup(EOSMAP_GAMC,EOS_IN,gridDataStruct)
  game_map = eos_mapLookup(EOSMAP_GAME,EOS_IN,gridDataStruct)
  eint_map = eos_mapLookup(EOSMAP_EINT,EOS_IN,gridDataStruct)
  ener_map = eos_mapLookup(EOSMAP_ENER,EOS_IN,gridDataStruct)
  velx_map = eos_mapLookup(EOSMAP_VELX,EOS_IN,gridDataStruct)
  vely_map = eos_mapLookup(EOSMAP_VELY,EOS_IN,gridDataStruct)
  velz_map = eos_mapLookup(EOSMAP_VELZ,EOS_IN,gridDataStruct)
  sumy_map = eos_mapLookup(EOSMAP_SUMY,EOS_IN,gridDataStruct)
  ye_map   = eos_mapLookup(EOSMAP_YE,  EOS_IN,gridDataStruct)
  entr_map = eos_mapLookup(EOSMAP_ENTR,EOS_IN,gridDataStruct)

  if(gridDataStruct == SCRATCH) then
     call Driver_abort("Eos_getData : the use of SCRATCH is deprecated")
  end if

  if(present(massFrac)) then
     m=1
     do k = kb,ke
        do j = jb,je
           do i = ib,ie
              do n = SPECIES_BEGIN,SPECIES_END
                 massFrac(m) = solnData(n,i,j,k)
                 m=m+1
              end do
           end do
        end do
     end do
  end if

  n = 0
  do k = kb,ke
     do j = jb,je
        do i = ib,ie
           if (velx_map > 0 .AND. vely_map > 0 .AND. velz_map > 0) then
              kineticEnergy  = 0.5*(solnData(velx_map,i,j,k)**2 + &
                                    solnData(vely_map,i,j,k)**2 + &
                                    solnData(velz_map,i,j,k)**2)
           else
              kineticEnergy = 0.0
           end if
           
           n=n+1
           eosData(ekin+n) = kineticEnergy
           !! kineticEnergy holds velocity vector information -- 1/2 * Vmag**2
           !! internalEnergy holds eint (directly)  or energyTotal - ekinetic (calculated),
           !!          depending upon eintSwitch
           if(eint_map /= NONEXISTENT) then
              internalEnergy  = solnData(eint_map,i,j,k)
              if(ener_map /= NONEXISTENT) then
                 if ( solnData(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
                    internalEnergy = solnData(ener_map,i,j,k) - kineticEnergy
                 end if
              end if
           else if(game_map /= NONEXISTENT) then ! This case should be usable for R(elativistic)HD - KW
              internalEnergy  = solnData(pres_map,i,j,k) / solnData(dens_map,i,j,k) / &
                                   (solnData(game_map,i,j,k) - 1.0)
              if(ener_map /= NONEXISTENT) then
                 if ( solnData(ener_map,i,j,k) - kineticEnergy > max(eos_smalle, eos_eintSwitch*kineticEnergy)) then
                    internalEnergy = solnData(ener_map,i,j,k) - kineticEnergy
                 end if
              end if
           else if(ener_map /= NONEXISTENT) then
              internalEnergy = solnData(ener_map,i,j,k)-kineticEnergy
           else
              internalEnergy = eos_smalle
           endif
           
           internalEnergy = max(internalEnergy, eos_smalle)
           eosData(eint+n) = internalEnergy
           
           eosData(pres+n) = solnData(pres_map,i,j,k)
           eosData(dens+n) = solnData(dens_map,i,j,k)
           eosData(temp+n) = solnData(temp_map,i,j,k)
           eosData(gamc+n) = solnData(gamc_map,i,j,k)

           if (ye_map /= NONEXISTENT) then
!!!!             eosData(abar+n) =  1.0d0
!!!!             eosData(zbar+n) = solnData(ye_map,i,j,k)
             eosData(elef+n) = solnData(ye_map,i,j,k)
           endif

           if((ye_map /= NONEXISTENT).and.(sumy_map /= NONEXISTENT)) then
              !! cal says abar=1/sumy
              !! cal says zbar=ye / sumy and he claims sumy are never zero
              eosData(abar+n) =  1.0 /  solnData(sumy_map,i,j,k)
              eosData(zbar+n) = solnData(ye_map,i,j,k) /  solnData(sumy_map,i,j,k)
              !PRINT*,'sumy_map'
           endif
           if(entr_map /= NONEXISTENT) eosData(entr+n) = solnData(entr_map,i,j,k)
        end do
     end do
  end do
  
  return
end subroutine Eos_getData 
