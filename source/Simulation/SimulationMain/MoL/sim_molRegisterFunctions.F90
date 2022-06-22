!!****if* source/Simulation/SimulationMain/MoL/sim_molRegisterFunctions
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
!!  NAME 
!!
!!      sim_molRegisterFunctions
!!
!!  SYNOPSIS
!!
!!      call sim_molRegisterFunctions
!!
!!  DESCRIPTION 
!!
!!      Registers various MoL functions
!!
!!***
subroutine sim_molRegisterFunctions()
    use sim_molInterface, only: sim_molExplicitRHS,     &
                                sim_molImplicitRHS,     &
                                sim_molImplicitRHS,     &
                                sim_molImplicitUpdate,  &
                                sim_molPostUpdate,      &
                                sim_molPostFastUpdate

    use MoL_interface, only: MoL_registerFunction

#include "MoL.h"

    implicit none

    call MoL_registerFunction(MOL_RHS_EXPLICIT,     sim_molExplicitRHS)
    call MoL_registerFunction(MOL_RHS_IMPLICIT,     sim_molImplicitRHS)
    call MoL_registerFunction(MOL_RHS_FAST,         sim_molFastRHS)
    call MoL_registerFunction(MOL_IMPLICIT_UPDATE,  sim_molImplicitUpdate)
    call MoL_registerFunction(MOL_POST_UPDATE,      sim_molPostUpdate)
    call MoL_registerFunction(MOL_POST_UPDATE_FAST, sim_molPostFastUpdate)
end subroutine sim_molRegisterFunctions
