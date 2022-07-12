!!****if* source/numericalTools/MoL/MoLMain/MoL_releaseFunctions
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
!!      MoL_releaseFunctions
!!
!!  SYNOPSIS
!!
!!      call MoL_releaseFunctions()
!!
!!  DESCRIPTION
!!
!!      Release all registered MoL functions
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine MoL_releaseFunctions()
    use MoL_functions

    implicit none

    if (associated(MoL_rhsE)) nullify(MoL_rhsE)
    if (associated(MoL_rhsI)) nullify(MoL_rhsI)
    if (associated(MoL_rhsF)) nullify(MoL_rhsF)

    if (associated(MoL_implicitUpdate)) nullify(MoL_implicitUpdate)
    
    if (associated(MoL_postUpdate)) nullify(MoL_postUpdate)
    if (associated(MoL_postUpdateFast)) nullify(MoL_postUpdateFast)

    MoL_rhsE => MoL_rhsE_default
    MoL_rhsI => MoL_rhsI_default
    MoL_rhsF => MoL_rhsF_default

    MoL_implicitUpdate => MoL_implicitUpdate_default

    MoL_postUpdate => MoL_postUpdate_default
    MoL_postUpdateFast => MoL_postUpdateFast_default
end subroutine MoL_releaseFunctions
