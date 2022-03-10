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
  use eos_localInterface, ONLY : eos_weaklib
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
     call eos_weaklib(mode, vecLen, eosData, massFrac, mask)
  elseif (pMassFrac) then
     call eos_weaklib(mode, vecLen, eosData, massFrac)
  elseif (pMask) then
     call eos_weaklib(mode, vecLen, eosData, mask=mask)
  else
     call eos_weaklib(mode, vecLen, eosData)
  end if
  return
end subroutine Eos
