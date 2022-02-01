!!****if* source/Simulation/SimulationMain/HydroStatic/Simulation_guardCellMaskHook
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
!!  call Simulation_guardCellMaskHook(logical(INOUT),dimension(*) :: ccmask,
!!                                    logical(IN)                 :: needeos)
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

#include "Simulation.h"

subroutine Simulation_guardCellMaskHook(ccMask, needEos)
  implicit none
  logical,intent(INOUT) :: ccMask(*)
  logical,intent(IN)    :: needEos

  !!  Additional logic necessary due to the hydrostatic boundary conditions.
  ! For constant isothermal BC we need temperature, as well as material info
  ! (such as what is stored in mass fractions or mass scalars used in some
  ! configurations), in the interior.
  ! This temperature setting does not trigger the "needEos" flag.
  if ( ccMask(DENS_VAR) .or. ccMask(EINT_VAR) .or. ccMask(TEMP_VAR) &
       .or. ccMask(PRES_VAR) .or. ccMask(ENER_VAR) &
       .or. ccMask(GAMC_VAR) .or. ccMask(GAME_VAR) ) then
     ccMask(DENS_VAR) = .true.
     ccMask(TEMP_VAR) = .true.
     ccMask(SPECIES_BEGIN:NUNK_VARS) = .true.
  endif

end subroutine Simulation_guardCellMaskHook

