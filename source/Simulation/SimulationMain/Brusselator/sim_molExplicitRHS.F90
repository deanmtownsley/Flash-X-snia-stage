!!****if* source/Simulation/SimulationMain/Brusselator/sim_molExplicitRHS
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
!!      sim_molExplicitRHS
!!
!!  SYNOPSIS
!!
!!      call sim_molExplicitRHS(Grid_tile_t, intent(in) :: tileDesc
!!                              real, pointer           :: rhs(:,:,:,:)
!!                              real, pointer           :: vars(:,:,:,:)
!!                              real, intent(in)        :: t)
!!
!!  DESCRIPTION 
!!
!!      Calculate explicit RHS terms
!!
!!
!!  ARGUMENTS
!!
!!      tileDesc : Current tile descriptor
!!      rhs      : Pointer to the RHS storage to fill
!!      vars     : Pointer to the current value of the evolved variables
!!      t        : Current time
!!
!!***
subroutine sim_molExplicitRHS(tileDesc, rhs, vars, t)
    use Simulation_data

    use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

    implicit none

    class(Grid_tile_t), intent(in) :: tileDesc
    real, dimension(:,:,:,:), pointer :: rhs, vars
    real, intent(in) :: t

    integer :: i, j, k, ip, im

    integer, dimension(LOW:HIGH,MDIM) :: lim, bcs
    real :: del(MDIM), idx

    if (sim_rho .lt. 0d0) then
        ip = 0
        im = -1
    else
        ip = 1
        im = 0
    end if

    call tileDesc%faceBCs(bcs)

    lim = tileDesc%limits

    if(bcs(LOW,IAXIS) .ne. NOT_BOUNDARY) lim(LOW,IAXIS) = lim(LOW,IAXIS) + 1
    if(bcs(HIGH,IAXIS) .ne. NOT_BOUNDARY) lim(HIGH,IAXIS) = lim(HIGH,IAXIS) - 1

    call tileDesc%deltas(del)
    idx = 1d0/del(IAXIS)

    do k = lim(LOW,KAXIS), lim(HIGH,KAXIS)
        do j = lim(LOW,JAXIS), lim(HIGH,JAXIS)
            do i = lim(LOW,IAXIS), lim(HIGH,IAXIS)
                rhs(U_RHS,i,j,k) = rhs(U_RHS,i,j,k) &
                    + sim_rho*(vars(U_VAR,i+ip,j,k) - vars(U_VAR,i+im,j,k))*idx

                rhs(V_RHS,i,j,k) = rhs(V_RHS,i,j,k) &
                    + sim_rho*(vars(V_VAR,i+ip,j,k) - vars(V_VAR,i+im,j,k))*idx

                rhs(W_RHS,i,j,k) = rhs(W_RHS,i,j,k) &
                    + sim_rho*(vars(W_VAR,i+ip,j,k) - vars(W_VAR,i+im,j,k))*idx
            end do ! i
        end do ! j
    end do ! k
end subroutine sim_molExplicitRHS
