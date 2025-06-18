!!****if* source/physics/Eos/EosMain/Helmholtz/ExternalAbarZbar/Eos_wrapped
!! NAME
!!
!!  Eos_wrapped
!! 
!! SYNOPSIS
!!
!!  call Eos_wrapped(  integer(IN) :: mode,
!!                     integer(IN) :: range(HIGH, MDIM),
!!                     integer(IN) :: blockID,
!!            optional,integer(IN) :: gridDataStruct )
!!
!! DESCRIPTION
!!
!! Josh Martin 2025
!!
!! This function is provided for the user's convenience and acts as a simple
!! wrapper to the Eos interface. The Eos interface uses a single, flexible data
!! structure "eosData" to pass the thermodynamic quantities in and out of the
!! funtion (see Eos). The wrapper hides formation and use of eosData
!! from the users.
!!
!! While Eos does not know anything about blocks, Eos_wrapped takes its
!! input thermodynamic state variables from a given block's storage area.
!! It works by taking a selected section of a block described by array
!! "range" and translating it to eosData before calling the Eos routine.
!! Upon return from Eos, Eos_wrapper updates certain state variables in
!! the same section of the block's storage area. Which variables are taken
!! as input, and which are updated, depends on the "mode" argument.
!!
!! If you want to return the derived quantities defined from EOS_VAR+1:EOS_NUM
!! in Eos.h, then you must use the direct interface Eos().
!!
!!
!!  ARGUMENTS 
!!
!!   
!!   mode : determines which variables are used as Eos input.
!!          The valid values are MODE_DENS_EI (where density and internal
!!          energy are inputs), MODE_DENS_PRES (density and pressure as inputs)
!!          MODE_DENS_TEMP (density and temperature are inputs).
!!          These quantities are defined in constants.h, the argument is 
!!          forwarded unchanged to the Eos function call.
!!          Note that internal energy is grid variable EINT_VAR, not ENER_VAR.
!!
!! 
!!   range: an array that holds the lower and upper indices of the section
!!          of block on which Eos is to be applies. The example shows how
!!          the array describes the block section.
!!
!!   blockID: current block number
!!
!!   gridDataStruct : the grid data structure on whose data Eos is to be applied
!!
!!
!!  EXAMPLE 
!!      if range(LOW,IAXIS)=1,range(HIGH,IAXIS)=iguard,
!!         range(LOW,JAXIS)=1,range(HIGH,JAXIS)=jguard,
!!         range(LOW,KAXIS)=1,range(HIGH,KAXIS)=kguard,
!!      then Eos is applied to the lower left hand corner of the guard
!!      cells in the block. 
!!
!!      However, if the value were
!!         range(LOW,IAXIS)=iguard+1,range(HIGH,IAXIS)=iguard+nxb,
!!         range(LOW,JAXIS)=jguard+1,range(HIGH,JAXIS)=jguard+nyb,
!!         range(LOW,KAXIS)=kguard+1,range(HIGH,KAXIS)=kguard+nzb,
!!      then Eos is applied to all the interior cells in the block.
!!
!!  NOTES
!!      This interface is defined in Fortran Module 
!!      Eos_interface. All functions calling this routine should include
!!      a statement like
!!      use Eos_interface, ONLY : Eos_wrapped
!!
!!      This routine cannot use "INTERIOR" mode of indexing the range.  In the
!!      second example given above, although only the interior cells are being
!!      calculated with EOS, the range indices still must include the guard cells.
!!      See, for example, IsentropicVortex/Simulation_initBlock where the data is
!!      generated on INTERIOR cells with Grid_putRowData, but the same indices can't
!!      be used for the EOS call.
!!
!!  SEE ALSO
!!
!!     Eos
!!     Eos.h
!!
!!***

! solnData depends on the ordering on unk
!!REORDER(4): solnData


!JM subroutine Eos_multiDim(mode,range,blockID,gridDataStruct)
subroutine Eos_multiDim(mode,range,solnData)

   use Eos_data, ONLY: eos_eintSwitch, eos_smalle, eos_meshMe
   use Driver_interface, ONLY : Driver_abort
   !JMJN use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
   use Logfile_interface, ONLY: Logfile_stampMessage 
   !JM use Eos_interface, ONLY : Eos
   use Eos_interface, ONLY : Eos_vector
   use eos_localInterface, ONLY : eos_externalComputeAbarZbar
 
   implicit none
 
 #include "Eos.h"
 #include "constants.h"
 #include "Simulation.h"
 
   integer, intent(in) :: mode
   integer, dimension(2,MDIM), intent(in) :: range
  !JMJN integer,intent(in) :: blockID
  !JM integer, optional, intent(IN) :: gridDataStruct
 
   real, pointer:: solnData(:,:,:,:)
 
 !JM #ifndef FIXEDBLOCKSIZE
   real, allocatable, dimension(:):: energyKinetic,energyInternal
   !JM real, allocatable :: eosData(:),massFraction(:)
   real, allocatable :: eosData(:,:) !JM
   integer, allocatable,dimension(:) :: iFlag
 !JM #else
 !JM  real, dimension(MAXCELLS):: energyKinetic,energyInternal
   !JM real, dimension(NSPECIES*MAXCELLS) :: massFraction
   !JM real, dimension(EOS_NUM*MAXCELLS) :: eosData
 !JM   real, dimension(MAXCELLS, EOS_NUM) :: eosData
 !JM   integer, dimension(MAXCELLS) :: iFlag
 !JM #endif
 
 
   integer :: ierr, dataStruct
   integer :: i,j,k, vecLen,pres,dens,gamc,temp,abar,zbar,eint,entr
 
   integer :: index, index_start, index_end !JMJN
 
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
 !JM   if (present(gridDataStruct)) then
 !JM      if (gridDataStruct .NE. CENTER) then
 !JM         call Driver_abort("Eos_wrapped: Support for gridDataStruct other than CENTER not implemented in this version!")
 !JM      end if
 !JM   end if
 #endif
 
   if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *
 
   ! Initializations:   grab the solution data from UNK and determine
   !   the length of the data being operated upon
 
 !JM   if(present(gridDataStruct))then
 !JM      dataStruct=gridDataStruct
 !JM   else
 !JM      dataStruct=CENTER
 !JM   end if
 
   !JM call Grid_getBlkPtr(blockID,solnData,dataStruct)
   !JM vecLen = range(HIGH,IAXIS)-range(LOW,IAXIS)+1 (I don't know why this would only be 1D)
 
    !JM I think VecLen is the number of rows in eosData
    vecLen = (range(HIGH,IAXIS)-range(LOW,IAXIS)+1)*& !JM
             (range(HIGH,JAXIS)-range(LOW,JAXIS)+1)*& !JM
             (range(HIGH,KAXIS)-range(LOW,KAXIS)+1) !JM
 
    if (vecLen==0) return !JM * Return immediately for empty IAXIS range! (for efficiency and avoiding index range errors)
 
   
   ! These integers are indexes into the location in eosData just before the storage area for the appropriate variable.
   !JM pres = (EOS_PRES-1) !JM *vecLen
   !JM dens = (EOS_DENS-1) !JM *vecLen
   !JM temp = (EOS_TEMP-1) !JM *vecLen
   !JM gamc = (EOS_GAMC-1) !JM *vecLen
   !JM eint = (EOS_EINT-1) !JM *vecLen
   !JM abar = (EOS_ABAR-1) !JM *vecLen
   !JM zbar = (EOS_ZBAR-1) !JM *vecLen
   !JM entr = (EOS_ENTR-1) !JM *vecLen
 
 !#ifndef FIXEDBLOCKSIZE
 allocate(energyInternal(vecLen))
 allocate(energyKinetic(vecLen))
 !JM allocate(massFraction(NSPECIES*vecLen))
 !JM allocate(eosData(EOS_NUM*vecLen))
 allocate(eosData(vecLen, EOS_NUM)) !DQ What is the difference between EOS_NUM and EOS_VARS?
 allocate(iFlag(vecLen)) !JM I still need to change this?
 !#endif  
 
   !JM do k = range(LOW,KAXIS), range(HIGH,KAXIS)
      !JM do j = range(LOW,JAXIS), range(HIGH,JAXIS)
 
         !! Fill up two scratch arrays. 
         !! energyKinetic holds velocity vector information -- 1/2 * Vmag**2
         !! energyInternal holds eint (directly)  or energyTotal - ekinetic (calculated),
         !!          depending upon eintSwitch
 
         !JM energyKinetic(1:vecLen) = 0.5*(solnData(VELX_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
     !JM &                                  solnData(VELY_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2 + &
     !JM &                                  solnData(VELZ_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)**2)
 
 
 !JM My version of energyKinetic (per gram) and energyInternal (per gram)
 
 !JM Creating an index to iterate through the for loops. index is actually local cell number.
 index = 1
 
 do k = range(LOW,KAXIS), range(HIGH,KAXIS)
    do j = range(LOW,JAXIS), range(HIGH,JAXIS)
       index_start = index
 
       do i = range(LOW,IAXIS), range(HIGH,IAXIS)
                   
          energyKinetic(index) = 0.5*(solnData(VELX_VAR,i,j,k)**2 + &
          &                              solnData(VELY_VAR,i,j,k)**2 + &
          &                              solnData(VELZ_VAR,i,j,k)**2)
 #ifdef EINT_VAR
             !JM initial values of EINT that may be overwritten
             !DQ Does having to access and then reaccess the array like this slow the program down?
             !JM energyInternal(index) = solnData(EINT_VAR,i,j,k)
             !JM Instead of doing potential overwriting, I am just doing a classic if/then/else
 
             !JM Here comes the potential overwriting. EInt should be
             !JM ETotal - Ekinetic, but this can only be calculated
             !JM to good precision if ETotal and EKinetic are different
             !JM enough. (Have Dean check this understanding). !DQ
             if (solnData(ENER_VAR,i,j,k) > &
                (1.+ eos_eintSwitch)*energyKinetic(index)) then
                energyInternal(index) = solnData(ENER_VAR,i,j,k) - energyKinetic(index)
             else 
                energyInternal(index) = solnData(EINT_VAR,i,j,k)
             end if
             !JM bringing Eint up to an Efloor if necessary
             energyInternal(index) = max(energyInternal(index), eos_smalle)
 #else
             !JM So now if we don't have the EINT_VAR stored, we just
             !JM simply do the E_T - E_kin and then check that E_int is
             !JM above the energy floor of our choosing.
             energyInternal(index) = solnData(ENER_VAR,i,j,k) - energyKinetic(index)
             energyInternal(index) = max(energyInternal(index), eos_smalle)
 #endif
 
          !JM Now grabbing all the grid data and populating eosData
          eosData(index, EOS_PRES) = solnData(PRES_VAR,i,j,k)
          eosData(index, EOS_DENS) = solnData(DENS_VAR,i,j,k)
          eosData(index, EOS_TEMP) = solnData(TEMP_VAR,i,j,k)
          eosData(index, EOS_GAMC) = solnData(GAMC_VAR,i,j,k)
          eosData(index, EOS_EINT) = energyInternal(index)
 
          if (index > vecLen) then !DEBUG
             call Driver_abort("Index exceeds vecLen in Eos_multiDim!") !DEBUG
          end if !DEBUG
          
          index = index + 1
 
       end do
 
       index_end = index - 1
 
       ! Call out for external computation of abar and zbar based on
       ! other scalars.  This calls the callback for whoever
       ! requested the ExternalAbarZbar implementation.
 
       call eos_externalComputeAbarZbar(solnData(SPECIES_BEGIN:UNK_VARS_END,&
       range(LOW,IAXIS):range(HIGH,IAXIS),j,k), &
       eosData(index_start:index_end, EOS_ABAR), eosData(index_start:index_end, EOS_ZBAR) )
 
    end do
 end do
 
 !JM Now that eosData has obtained all necessary data from grid, we
 !JM call Eos_vector
 
 call Eos_vector(mode,vecLen,eosData)
 
 !JM Now that Eos_vector has been called, we need to write
 !JM the completed Eos variable calculations (eosData) back
 !JM to the grid (which is solnData).
 
 index = 1
 iFlag = 0 !DQ is this like iFlag = np.zeros(vecLen)?
 
 do k = range(LOW,KAXIS), range(HIGH,KAXIS)
    do j = range(LOW,JAXIS), range(HIGH,JAXIS)
       do i = range(LOW,IAXIS), range(HIGH,IAXIS)
          
          solnData(PRES_VAR,i,j,k) = eosData(index, EOS_PRES)
          solnData(TEMP_VAR,i,j,k) = eosData(index, EOS_TEMP)
          solnData(GAMC_VAR,i,j,k) = eosData(index, EOS_GAMC)
 
 #ifdef EINT_VAR
          solnData(EINT_VAR,i,j,k) = eosData(index, EOS_EINT)
 #endif
          solnData(ENER_VAR,i,j,k) = eosData(index, EOS_EINT) + energyKinetic(index)
 #ifdef ENTR_VAR
          solnData(ENTR_VAR,i,j,k) = eosData(index, EOS_ENTR)
 #endif
 
          ! check for zero values before calculating gamma
          !JM I changed this where statement to an if/then statement
          if ( (eosData(index, EOS_EINT) .eq. 0.) .or. (eosData(index, EOS_DENS) .eq. 0.)) then
             iFlag(index) = 1
          end if
          
          if (index > vecLen) then !DEBUG
             call Driver_abort("Index exceeds vecLen in Eos_multiDim!") !DEBUG
          end if !DEBUG
 
          index = index + 1
          
       end do
    end do
 end do
 
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
 
 
 index = 1
 
 do k = range(LOW,KAXIS), range(HIGH,KAXIS)
    do j = range(LOW,JAXIS), range(HIGH,JAXIS)
       do i = range(LOW,IAXIS), range(HIGH,IAXIS)
 
          !JM calculating gamma_e !DQ What is gamma_e?
          solnData(GAME_VAR,i,j,k) = eosData(index, EOS_PRES)/&
          (eosData(index, EOS_EINT) *eosData(index, EOS_DENS)) + 1.0
 
          if (index > vecLen) then !DEBUG
             call Driver_abort("Index exceeds vecLen in Eos_multiDim!") !DEBUG
          end if !DEBUG
 
          index = index + 1
 
       end do
    end do
 end do
                
 
 
 !JM #ifdef EINT_VAR
 !JM          energyInternal(1:vecLen) = solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
 !JM 
 !JM          do i = 1,vecLen
 !JM             if (solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) > &
 !JM                (1.+ eos_eintSwitch)*energyKinetic(i)) then
 !JM                energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - energyKinetic(i)
 !JM             end if
 !JM             energyInternal(i) = max(energyInternal(i), eos_smalle)
 !JM             massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
 !JM                solnData(SPECIES_BEGIN:SPECIES_END,range(LOW,IAXIS)+i-1,j,k)
 !JM          end do
 !JM #else
 !JM         do i = 1,vecLen
 !JM            energyInternal(i) = solnData(ENER_VAR,range(LOW,IAXIS)+i-1,j,k) - energyKinetic(i)
 !JM            energyInternal(i) = max(energyInternal(i), eos_smalle)
 !JM            massFraction((i-1)*NSPECIES+1:i*NSPECIES) = &
 !JM                 solnData(SPECIES_BEGIN:SPECIES_END,range(LOW,IAXIS)+i-1,j,k)
 !JM         end do
 !JM #endif
 
 !JM         eosData(pres+1:pres+vecLen) = &
 !JM              solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
 !JM         eosData(dens+1:dens+vecLen) = &
 !JM              solnData(DENS_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
 !JM         eosData(temp+1:temp+vecLen) = &
 !JM              solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
 !JM         eosData(gamc+1:gamc+vecLen) = &
 !JM              solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k)
 !JM         eosData(eint+1:eint+vecLen) = energyInternal(1:vecLen)
 
         ! Call out for external computation of abar and zbar based on
         ! other scalars.  This calls the callback for whoever
         ! requested the ExternalAbarZbar implementation.
 !JM         call eos_externalComputeAbarZbar(solnData(SPECIES_BEGIN:UNK_VARS_END,&
 !JM              range(LOW,IAXIS):range(HIGH,IAXIS),j,k), &
 !JM              eosData(abar+1:abar+vecLen), eosData(zbar+1:zbar+vecLen) )
 
 !JM         call Eos(mode,vecLen,eosData,massFraction)
         
 !JM         solnData(PRES_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(pres+1:pres+vecLen)
 !JM         solnData(TEMP_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(temp+1:temp+vecLen)
 !JM         solnData(GAMC_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(gamc+1:gamc+vecLen)
 !JM #ifdef EINT_VAR
 !JM         solnData(EINT_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(eint+1:eint+veclen)
 !JM #endif
 !JM         solnData(ENER_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(eint+1:eint+veclen) + energyKinetic(1:vecLen)
 !JM #ifdef ENTR_VAR
 !JM         solnData(ENTR_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(entr+1:entr+veclen)
 !JM #endif
 
         ! check for zero values before calculating gamma
 !JM         iFlag = 0
 !JM         where ( (eosData(eint+1:eint+vecLen) .eq. 0.) .or. (eosData(dens+1:dens+vecLen) .eq. 0.))
 !JM            iFlag(1:vecLen) = 1
 !JM         end where
 
         !maybe there was a wrong flag set
 !JM         if (maxval(iFlag) .gt. 0) then
 !JM            if (eos_meshMe .EQ. MASTER_PE) then
 !JM               write(*,*) "ERROR After calling Eos, eosData(EOS_EINT) or eosData(EOS_DENS) are zero"
 !JM               write(*,*) "  Perhaps the initialization routine is wrong..... or"
 !JM               write(*,*) "  perhaps the runtime parameter eosMode is wrong."
 !JM               write(*,*) "  This routine Eos_wrapped was called with mode= ", mode
 !JM               write(*,*) "     Check constants.h to determine value of MODE_DENS_??"
 !JM            endif
 !JM            call Logfile_stampMessage('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
 !JM            call Driver_abort('[Eos_wrapped] ERROR Density or Internal Energy are zero after a call to EOS!')
 !JM         end if
 
 !JM         solnData(GAME_VAR,range(LOW,IAXIS):range(HIGH,IAXIS),j,k) = &
 !JM              eosData(pres+1:pres+veclen)/&
 !JM              (eosData(eint+1:eint+veclen) *eosData(dens+1:dens+veclen)) + 1.0
 
 !JM      end do
 !JM   end do
 
 !JM call Grid_releaseBlkPtr(blockID,solnData,dataStruct)
 
 !JM #ifndef FIXEDBLOCKSIZE
   deallocate(energyKinetic)
   deallocate(energyInternal)
   deallocate(eosData)
   deallocate(iFlag)
 !JM  deallocate(massFraction)
 !JM #endif
 !JM  return
 end subroutine Eos_multiDim