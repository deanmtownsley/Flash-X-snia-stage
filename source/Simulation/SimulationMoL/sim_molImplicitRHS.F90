!!****f* source/Simulation/SimulationMoL/sim_molImplicitRHS
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
!!      sim_molImplicitRHS
!!
!!  SYNOPSIS
!!
!!      call sim_molImplicitRHS(Grid_tile_t, intent(in) :: tileDesc
!!                              real, pointer           :: rhs(:,:,:,:)
!!                              real, pointer           :: U(:,:,:,:)
!!                              real, intent(in)        :: t)
!!
!!  DESCRIPTION 
!!
!!      Calculate implicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      tileDesc : Current tile descriptor
!!      rhs      : Pointer to the RHS storage to fill
!!      U        : Pointer to the current value of the evolved variables
!!      t        : Current time
!!
!!***
subroutine sim_molImplicitRHS(tileDesc, rhs, U, t)
    use Grid_tile, only: Grid_tile_t

    implicit none

    class(Grid_tile_t), intent(in) :: tileDesc
    real, dimension(:,:,:,:), pointer :: rhs, U
    real, intent(in) :: t

    return
end subroutine sim_molImplicitRHS
