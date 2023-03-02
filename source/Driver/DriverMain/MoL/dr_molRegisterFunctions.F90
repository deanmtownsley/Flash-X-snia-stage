!!****if* source/Driver/DriverMain/MoL/dr_molRegisterFunctions
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!      dr_molRegisterFunctions
!!
!!  SYNOPSIS
!!
!!      call dr_molRegisterFunctions
!!
!!  DESCRIPTION
!!
!!      Registers various MoL functions
!!
!!***
subroutine dr_molRegisterFunctions()
   use dr_molInterface, only: dr_molExplicitRHS, &
                              dr_molImplicitRHS, &
                              dr_molFastRHS, &
                              dr_molImplicitUpdate, &
                              dr_molPostUpdate, &
                              dr_molPostFastUpdate

   use MoL_interface, only: MoL_registerRHS, MoL_registerUpdate, MoL_registerPostUpdate

#include "MoL.h"

   implicit none

   call MoL_registerRHS(MOL_RHS_EXPLICIT, dr_molExplicitRHS)
   call MoL_registerRHS(MOL_RHS_IMPLICIT, dr_molImplicitRHS)
   call MoL_registerRHS(MOL_RHS_FAST, dr_molFastRHS)

   call MoL_registerUpdate(MOL_IMPLICIT_UPDATE, dr_molImplicitUpdate)

   call MoL_registerPostUpdate(MOL_POST_UPDATE, dr_molPostUpdate)
   call MoL_registerPostUpdate(MOL_POST_UPDATE_FAST, dr_molPostFastUpdate)
end subroutine dr_molRegisterFunctions
