!!****if* source/Simulation/SimulationMain/AdvectDiffuseReact/sim_molFastRHS
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
!!      sim_molFastRHS
!!
!!  SYNOPSIS
!!
!!      call sim_molFastRHS(Grid_tile_t, intent(in) :: tileDesc
!!                          real, pointer           :: rhs(:,:,:,:)
!!                          real, pointer           :: U(:,:,:,:)
!!                          real, intent(in)        :: t)
!!
!!  DESCRIPTION 
!!
!!      Calculate fast RHS terms
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
subroutine sim_molFastRHS(tileDesc, rhs, vars, t)
    use Simulation_data, only: U_RHS, V_RHS, W_RHS, a => sim_a, b => sim_b, &
                               eps => sim_epsilon

    use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

    implicit none

    class(Grid_tile_t), intent(in) :: tileDesc
    real, dimension(:,:,:,:), pointer :: rhs, vars
    real, intent(in) :: t

    integer, dimension(LOW:HIGH,MDIM) :: lim, bcs
    integer :: i, j, k
    real :: u,v,w

    call tileDesc%faceBCs(bcs)

    lim = tileDesc%limits

    if(bcs(LOW,IAXIS) .ne. NOT_BOUNDARY) lim(LOW,IAXIS) = lim(LOW,IAXIS) + 1
    if(bcs(HIGH,IAXIS) .ne. NOT_BOUNDARY) lim(HIGH,IAXIS) = lim(HIGH,IAXIS) - 1

    do k = lim(LOW,KAXIS), lim(HIGH,KAXIS)
        do j = lim(LOW,JAXIS), lim(HIGH,JAXIS)
            do i = lim(LOW,IAXIS), lim(HIGH,IAXIS)
                u = vars(U_VAR,i,j,k)
                v = vars(V_VAR,i,j,k)
                w = vars(W_VAR,i,j,k)

                rhs(U_RHS,i,j,k) = rhs(U_RHS,i,j,k) + a         - (w+1d0)*u + v*u**2
                rhs(V_RHS,i,j,k) = rhs(V_RHS,i,j,k)             +       w*u - v*u**2
                rhs(W_RHS,i,j,k) = rhs(W_RHS,i,j,k) + (b-w)/eps -       w*u
            end do ! i
        end do ! j
    end do ! k
end subroutine sim_molFastRHS
