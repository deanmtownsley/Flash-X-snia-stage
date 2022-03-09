 !!****if* source/physics/Eos/EosMain/Eos
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
!!  Eos
!!
!!  For more details see the documentation of the NULL implementation
!!
!!***

subroutine Eos(mode, vecLen, eosData, massFrac, mask, vecBegin,vecEnd, diagFlag)

!==============================================================================
  use Driver_interface, ONLY : Driver_abort
  use Eos_data, ONLY : eos_meshMe, eos_type
  use eos_localInterface, ONLY : eos_idealGamma, eos_mgamma, eos_helmholtz,&
      eos_weaklib
  implicit none
#include "constants.h"
#include "Eos.h"
#include "Simulation.h"
  integer, INTENT(in) :: mode, vecLen
  real,INTENT(inout), dimension(EOS_NUM*vecLen) :: eosData 
  logical, optional, INTENT(in),target,dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, optional, INTENT(in),dimension(NSPECIES*vecLen)    :: massFrac
  integer,optional,INTENT(in) :: vecBegin,vecEnd
  integer, optional, INTENT(out)    :: diagFlag
  logical :: pMassFrac_and_mask, pMassFrac, pMask

  if (present(diagFlag)) diagFlag = 0

  if (mode==MODE_EOS_NOP) return ! * Return immediately for MODE_EOS_NOP! *
  if (mode==MODE_EOS_WRAPPERONLY) return ! * Return immediately for MODE_EOS_WRAPPERONLY! *

  pMassFrac = present(massFrac)
  pMask = present(mask)
  pMassFrac_and_mask = pMassFrac.and.pMask
  
  if(pMassFrac_and_mask) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, mask=mask, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, mask=mask)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, massFrac=massFrac, mask=mask)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, massFrac, mask)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  elseif (pMassFrac) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, massFrac=massFrac)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, massFrac=massFrac)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, massFrac)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  elseif (pMask) then
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, mask=mask, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd, mask=mask)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData, mask=mask)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData, mask=mask)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  else
     select case(eos_type)
     case(EOS_GAM)
        call eos_idealGamma(mode, vecLen, eosData, vecBegin,vecEnd, diagFlag=diagFlag)
     case(EOS_MGAM)
        call eos_mgamma(mode, vecLen, eosData, vecBegin,vecEnd)
     case(EOS_HLM)
        call eos_helmholtz(mode, vecLen, eosData)
     case(EOS_WL)
        call eos_weaklib(mode, vecLen, eosData)
     case default
        if (eos_meshMe==MASTER_PE) print*,'[Eos] unrecognized eos_type',eos_type
        call Driver_abort('[Eos] unrecognized eos_type.')
     end select
  end if
  return
end subroutine Eos
