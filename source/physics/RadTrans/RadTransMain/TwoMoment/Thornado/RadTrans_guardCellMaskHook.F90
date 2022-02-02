!!****if* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_guardCellMaskHook
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
!!  RadTrans_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call RadTrans_guardCellMaskHook(logical,intent(INOUT)  :: ccMask,
!!                                  logical,intent(IN)  :: needEos)
!!
!! DESCRIPTION
!! GC mask hook
!!
!! ARGUMENTS
!!
!!   ccmask : cc mask
!!
!!   needeos : check for eos
!!
!!
!!
!!***

subroutine RadTrans_guardCellMaskHook(ccMask, needEos)

  use Driver_interface, ONLY : Driver_abort
  use RadTrans_data, ONLY : rt_eosModeGc

  implicit none

#include "Simulation.h"

  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

#ifndef YE_MSCALAR
  call Driver_abort("RadTrans_guardCellMaskHook: YE_MSCALAR not defined.")
#else
  ccMask(YE_MSCALAR) = .TRUE.
#ifdef FLASH_EOS_WEAKLIB
  if(needEos) then
     ccMask(DENS_VAR) = .TRUE.
     ccMask(TEMP_VAR) = .TRUE.
  end if
#endif
#endif
#if NTHORNADO > 0
  ccMask(THORNADO_BEGIN:THORNADO_END) = .TRUE.
#endif
  !needEos = (.not.(rt_eosModeGc==MODE_DENS_EI)).and.(ccMask(ENER_VAR).or.ccMask(EINT_VAR))

end subroutine RadTrans_guardCellMaskHook

