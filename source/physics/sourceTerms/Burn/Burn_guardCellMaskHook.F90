!!****f* source/physics/sourceTerms/Burn/Burn_guardCellMaskHook
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
!!  Burn_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call Burn_guardCellMaskHook(logical,intent(INOUT)  :: ccMask,
!!                              logical,intent(IN)  :: needEos)
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

subroutine Burn_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

end subroutine Burn_guardCellMaskHook

