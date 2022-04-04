!!****if* source/physics/Eos/EosMain/Helmholtz/Ye/Eos_wrapped
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
!! This version works in EOS_YE mode.  Here, abar is calculated from
!! the Ye mass scalar and zbar from the Ye and SumY mass scalars,
!! rather than leaving these to be determined (at a lower level,
!! i.e. probably in the Eos() implementation) from multiple species
!! information (or from single species runtime parameters).
!! Code can check for EOS_YE mode by testing whether the preprocessor
!! symbol USE_EOS_YE is defined.
!!
!!  For more details see the documentation of the NULL implementation
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


subroutine Eos_wrapped(mode,range,solnData,gridDataStruct)

  use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_meshMe, &
       eos_threadWithinBlock
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Driver_interface, ONLY : Driver_abort

  !$ use omp_lib    
  use Eos_interface, ONLY : Eos

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Simulation.h"

  integer, intent(in) :: mode
  integer, dimension(2,MDIM), intent(in) :: range
  real, pointer:: solnData(:,:,:,:)

  integer,optional,intent(in) :: gridDataStruct


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


  integer :: ierr
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
#endif

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *

! Sanity check. Remove when gridDataStruct .NE. CENTER gets implemented, or move into the 'ifdef DEBUG'
! section above if performance is a concern.
  if (present(gridDataStruct)) then
     if (gridDataStruct .NE. CENTER) then
        call Driver_abort("Eos_wrapped: This gridDataStruct is not implemented in this version!")
     end if
  end if

  ! Initializations:   grab the solution data from UNK and determine
  !   the length of the data being operated upon

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

  !$omp parallel if (eos_threadWithinBlock .and. NDIM > 1) &
  !$omp default(none) &
  !$omp shared(mode,range,blockID,solnData,eos_eintSwitch,eos_smalle,&
  !$omp eos_meshMe,eos_threadWithinBlock) &
  !$omp firstprivate(vecLen,pres,dens,temp,gamc,eint,abar,zbar,entr) &
  !$omp private(i,j,k,energyInternal,energyKinetic,massFraction,eosData,iFlag)

#ifndef FIXEDBLOCKSIZE
  allocate(energyInternal(vecLen))
  allocate(energyKinetic(vecLen))
  allocate(massFraction(NSPECIES*vecLen))
  allocate(eosData(EOS_NUM*vecLen))
  allocate(iFlag(vecLen))
#endif  

#if NDIM==3
  !$omp do schedule(static)
#endif
  do k = range(LOW,KAXIS), range(HIGH,KAXIS)
#if defined(_OPENMP) && defined(DEBUG_THREADING) && NDIM==3
     if (eos_threadWithinBlock) then
        write (6,'(a,i3,a,i3,a,i3)') 'Thread', omp_get_thread_num(), &
             " of", omp_get_num_threads(), " assigned k loop iteration", k
     end if
#endif

#if NDIM==2
     !$omp do schedule(static)
#endif
     do j = range(LOW,JAXIS), range(HIGH,JAXIS)
#if defined(_OPENMP) && defined(DEBUG_THREADING) && NDIM==2
     if (eos_threadWithinBlock) then
        write (6,'(a,i3,a,i3,a,i3)') 'Thread', omp_get_thread_num(), &
             " of", omp_get_num_threads(), " assigned j loop iteration", j
     end if
#endif

     !! Fill up two scratch arrays. 
     !! energyKinetic holds velocity vector information -- 1/2 * Vmag**2
     !! energyInternal holds eint (directly)  or energyTotal - ekinetic (calculated),
     !!          depending upon eintSwitch
     ! NOTE this internal energy stuff is quite irrelevant for anything other than mode DENS_EI

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

        !!  These values may be overwritten within the internal Eos routines.
        !!  They are added here for cases where the mass fractions will not be used
        !!  to calculate abar and zbar.
#ifdef YE_MSCALAR
        !! cal says abar=1/sumy
        !! cal says zbar=ye / sumy and he claims sumy are never zero
        eosData(abar+1:abar+vecLen) =  1.0 / &
             solnData(SUMY_MSCALAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
        eosData(zbar+1:zbar+vecLen) = &
             solnData(YE_MSCALAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) / & 
             solnData(SUMY_MSCALAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
#else
        !! yeah, yeah this should be done in initialization, but trying to avoid too much code duplication
        call Driver_abort("[Eos_wrapped] This routine was called in YE mode, but no YE_MSCALAR is defined")
#endif

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

        ! the EOS_YE variables don't ever need to written back into the solution vector
        ! check for zero values
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
             (eosData(eint+1:eint+veclen) *eosData(dens+1:dens+veclen)) +1

     end do
#if NDIM==2
     !$omp end do nowait
#endif

  end do
#if NDIM==3
  !$omp end do nowait
#endif

#ifndef FIXEDBLOCKSIZE
  deallocate(energyKinetic)
  deallocate(energyInternal)
  deallocate(eosData)
  deallocate(iFlag)
  deallocate(massFraction)
#endif

  !$omp end parallel

  return
end subroutine Eos_wrapped



