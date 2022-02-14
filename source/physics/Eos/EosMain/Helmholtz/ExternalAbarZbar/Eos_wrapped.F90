!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/Eos_wrapped
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
!!  Eos_wrapped
!! 
!!  For more details see the documentation of the NULL implementation
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_wrapped(mode,range,blockID,gridDataStruct)

  use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_meshMe
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Logfile_interface, ONLY: Logfile_stampMessage 
  use Eos_interface, ONLY : Eos
  use eos_localInterface, ONLY : eos_externalComputeAbarZbar

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  integer,intent(in) :: blockID
  integer, optional, intent(IN) :: gridDataStruct

  real, pointer:: solnData(:,:,:,:)

#ifndef FIXEDBLOCKSIZE
  real, allocatable, dimension(:):: energyKinetic,energyInternal
  real, allocatable :: eosData(:),massFraction(:)
  integer, allocatable,dimension(:) :: iFlag
#else
  real, dimension(MAXCELLS):: energyKinetic,energyInternal
  real, dimension(NSPECIES*MAXCELLS) :: massFraction
  real, dimension(EOS_NUM*MAXCELLS) :: eosData
  integer, dimension(MAXCELLS) :: iFlag
#endif


  integer :: ierr, dataStruct
  integer :: i,j,k, vecLen,pres,dens,gamc,temp,abar,zbar,eint,entr

!! ---------------------------------------------------------------------------------
  ! Test calling arguments
#ifdef DEBUG
  ierr = 1
  select case (mode)
  case (MODE_DENS_PRES)
     ierr = 0
  case (MODE_DENS_TEMP)
     ierr = 0
  case (MODE_DENS_EI)
     ierr = 0
  end select

  if(ierr /= 0) then
     call Driver_abort("[Eos_wrapped] invalid mode: must be MODE_DENS_PRES, MODE_DENS_TEMP, or MODE_DENSE_EI")
  end if

! Sanity check
  if (present(gridDataStruct)) then
     if (gridDataStruct .NE. CENTER) then
        call Driver_abort("Eos_wrapped: Support for gridDataStruct other than CENTER not implemented in this version!")
     end if
  end if
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

  ! Initializations:   grab the solution data from UNK and determine
  !   the length of the data being operated upon

  if(present(gridDataStruct))then
     dataStruct=gridDataStruct
  else
     dataStruct=CENTER
  end if

  call Grid_getBlkPtr(blockID,solnData,dataStruct)
  vecLen = range(HIGH,IAXIS)-range(LOW,IAXIS)+1

  ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  gamc = (EOS_GAMC-1)*vecLen
  eint = (EOS_EINT-1)*vecLen
  abar = (EOS_ABAR-1)*vecLen
  zbar = (EOS_ZBAR-1)*vecLen
  entr = (EOS_ENTR-1)*vecLen

#ifndef FIXEDBLOCKSIZE
  allocate(energyInternal(vecLen))
  allocate(energyKinetic(vecLen))
  allocate(massFraction(NSPECIES*vecLen))
  allocate(eosData(EOS_NUM*vecLen))
  allocate(iFlag(vecLen))
#endif  

  do k = range(LOW,KAXIS), range(HIGH,KAXIS)
     do j = range(LOW,JAXIS), range(HIGH,JAXIS)

        !! Fill up two scratch arrays. 
        !! energyKinetic holds velocity vector information -- 1/2 * Vmag**2
        !! energyInternal holds eint (directly)  or energyTotal - ekinetic (calculated),
        !!          depending upon eintSwitch

        energyKinetic(1:vecLen) = 0.5*(solnData(VELX_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
    &                                  solnData(VELY_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
    &                                  solnData(VELZ_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2)


#ifdef EINT_VAR
        energyInternal(1:vecLen) = solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)

        do i = 1,vecLen
           if (solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) > &
                (1.+ eos_eintSwitch)*energyKinetic(i)) then
              energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - energyKinetic(i)
           end if
           energyInternal(i) = max(energyInternal(i), eos_smalle)
           massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
                solnData(SPECIES_BEGIN:SPECIES_END,range(LOW,IAXIS)+i-1,j,k)
        end do
#else
        do i = 1,vecLen
           energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - energyKinetic(i)
           energyInternal(i) = max(energyInternal(i), eos_smalle)
           massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
                solnData(SPECIES_BEGIN:SPECIES_END,range(LOW,IAXIS)+i-1,j,k)
        end do
#endif

        eosData(pres+1:pres+vecLen) = &
             solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
        eosData(dens+1:dens+vecLen) = &
             solnData(DENS_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
        eosData(temp+1:temp+vecLen) = &
             solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
        eosData(gamc+1:gamc+vecLen) = &
             solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
        eosData(eint+1:eint+vecLen) = energyInternal(1:vecLen)

        ! Call out for external computation of abar and zbar based on
        ! other scalars.  This calls the callback for whoever
        ! requested the ExternalAbarZbar implementation.
        call eos_externalComputeAbarZbar(solnData(SPECIES_BEGIN:UNK_VARS_END,&
             range(LOW,IAXIS):range(HIGH,IAXIS),j,k), &
             eosData(abar+1:abar+vecLen), eosData(zbar+1:zbar+vecLen) )

        call Eos(mode,vecLen,eosData,massFraction)
        
        solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(pres+1:pres+vecLen)
        solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(temp+1:temp+vecLen)
        solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(gamc+1:gamc+vecLen)
#ifdef EINT_VAR
        solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(eint+1:eint+veclen)
#endif
        solnData(ENER_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(eint+1:eint+veclen) + energyKinetic(1:vecLen)
#ifdef ENTR_VAR
        solnData(ENTR_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(entr+1:entr+veclen)
#endif

        ! check for zero values before calculating gamma
        iFlag = 0
        where ( (eosData(eint+1:eint+vecLen) .eq. 0.) .or. (eosData(dens+1:dens+vecLen) .eq. 0.))
           iFlag(1:vecLen) = 1
        end where

        !maybe there was a wrong flag set
        if (maxval(iFlag) .gt. 0) then
           if (eos_meshMe .EQ. MASTER_PE) then
              write(*,*) "ERROR After calling Eos, eosData(EOS_EINT) or eosData(EOS_DENS) are zero"
              write(*,*) "  Perhaps the initialization routine is wrong..... or"
              write(*,*) "  perhaps the runtime parameter eosMode is wrong."
              write(*,*) "  This routine Eos_wrapped was called with mode= ", mode
              write(*,*) "     Check constants.h to determine value of MODE_DENS_??"
           endif
           call Logfile_stampMessage('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
           call Driver_abort('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
        end if

        solnData(GAME_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
             eosData(pres+1:pres+veclen)/&
             (eosData(eint+1:eint+veclen) *eosData(dens+1:dens+veclen)) + 1.0

     end do
  end do
  call Grid_releaseBlkPtr(blockID,solnData,dataStruct)

#ifndef FIXEDBLOCKSIZE
  deallocate(energyKinetic)
  deallocate(energyInternal)
  deallocate(eosData)
  deallocate(iFlag)
  deallocate(massFraction)
#endif
  return
end subroutine Eos_wrapped



