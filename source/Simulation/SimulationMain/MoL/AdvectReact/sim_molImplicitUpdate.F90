!!****if* source/Simulation/SimulationMain/MoL/AdvectReact/sim_molImplicitUpdate
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
!!      sim_molImplicitUpdate
!!
!!  SYNOPSIS
!!
!!      call sim_molImplicitUpdate(real, intent(in) :: t
!!                                 real, intent(in) :: dt)
!!
!!  DESCRIPTION 
!!
!!      Implicitly update evolved variables from t to t+dt
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!      dt : Size of the time step to take
!!
!!***
subroutine sim_molImplicitUpdate(t, dt)
    use Simulation_data

    use Grid_interface, only : Grid_getTileIterator, Grid_releaseTileIterator
    use Grid_iterator,  only: Grid_iterator_t
    use Grid_tile,      only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

    implicit none

    real, intent(in) :: t, dt

    integer, dimension(LOW:HIGH,MDIM) :: lim
    integer :: i, j, k

    real, dimension(:,:,:,:), pointer :: U
    real :: fac

    type(Grid_tile_t) :: tileDesc
    type(Grid_iterator_t) :: itor

    ! No need to do this if decay term is zero
    if (sim_beta .eq. 0d0) return

    ! For the implicit update
    !   U^n+1 = U^n - dt*beta*U^n+1
    ! We can directly solve this (in-place) as
    !   U^n+1 = U^n / (1 + dt*beta)
    fac = 1d0/(1d0 + sim_beta*dt)

    call Grid_getTileIterator(itor, LEAF)

    do ! not using while since this is technically deprecated
        if (.not. itor%isValid()) exit

        call itor%currentTile(tileDesc)

        lim = tileDesc%limits

        call tileDesc%getDataPtr(U, CENTER)

        do k = lim(LOW,KAXIS), lim(HIGH,KAXIS)
            do j = lim(LOW,JAXIS), lim(HIGH,JAXIS)
                do i = lim(LOW,IAXIS), lim(HIGH,IAXIS)
                    U(PHI0_VAR,i,j,k) = U(PHI0_VAR,i,j,k)*fac
                    U(PHI1_VAR,i,j,k) = U(PHI1_VAR,i,j,k)*fac
                end do ! i
            end do ! j
        end do ! k

        call tileDesc%releaseDataPtr(U, CENTER)

        call itor%next()
    end do ! itor

    call Grid_releaseTileIterator(itor)
end subroutine sim_molImplicitUpdate
