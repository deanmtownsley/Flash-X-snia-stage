!!****f* source/Simulation/Simulation_guardCellMaskHook
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
!!  Simulation_guardCellMaskHook
!!
!! SYNOPSIS
!!
!!  call Simulation_guardCellMaskHook(logical(INOUT)  :: ccmask,
!!                                    logical(IN)  :: needeos)
!!
!! DESCRIPTION
!!
!!  A hook that lets a simulation modify the mask to use for guard cell filling.
!!
!!  Indirectly called from gr_makeMaskConsistent, which may get called from
!!  Grid_fillGuardCells (depending on the arguments with which Grid_fillGuardCells
!!  is called).
!!
!! ARGUMENTS
!!
!!   ccmask : the mask
!!
!!   needeos : switch for the need of Eos
!!
!!
!!***


subroutine Simulation_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  ! Stub does nothing.
end subroutine Simulation_guardCellMaskHook

