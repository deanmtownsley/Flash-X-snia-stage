!!****h* source/Simulation/SimulationMain/MoL/sim_molInterface
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
!!      sim_molInterface
!!
!!  SYNOPSIS
!!
!!      use sim_molInterface
!!
!!  DESCRIPTION
!!
!!      MoL-specific simulation features
!!
!!***
module sim_molInterface

    use Grid_tile, only: Grid_tile_t

    implicit none


    interface
        subroutine sim_molExplicitRHS(tileDesc, rhs, U, t)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(in) :: tileDesc
            real, dimension(:,:,:,:), pointer :: rhs, U
            real, intent(in) :: t
        end subroutine sim_molExplicitRHS
    end interface

    interface
        subroutine sim_molImplicitRHS(tileDesc, rhs, U, t)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(in) :: tileDesc
            real, dimension(:,:,:,:), pointer :: rhs, U
            real, intent(in) :: t
        end subroutine sim_molImplicitRHS
    end interface

    interface
        subroutine sim_molFastRHS(tileDesc, rhs, U, t)
            import :: Grid_tile_t
            class(Grid_tile_t), intent(in) :: tileDesc
            real, dimension(:,:,:,:), pointer :: rhs, U
            real, intent(in) :: t
        end subroutine sim_molFastRHS
    end interface

    interface
        subroutine sim_molImplicitUpdate(t, dt)
            real, intent(in) :: t, dt
        end subroutine sim_molImplicitUpdate
    end interface

    interface
        subroutine sim_molPostUpdate(t)
            real, intent(in) :: t
        end subroutine sim_molPostUpdate
    end interface

    interface
        subroutine sim_molPostFastUpdate(t)
            real, intent(in) :: t
        end subroutine sim_molPostFastUpdate
    end interface

    interface
        subroutine sim_molPreEvolve(t)
            real, intent(in) :: t
        end subroutine sim_molPreEvolve
    end interface

    interface
        subroutine sim_molRegisterFunctions
        end subroutine sim_molRegisterFunctions
    end interface

end module sim_molInterface
