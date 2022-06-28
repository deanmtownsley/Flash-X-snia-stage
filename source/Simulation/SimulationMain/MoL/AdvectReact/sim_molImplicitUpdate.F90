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
    real :: Uold(2), Unew(2)

    type(Grid_tile_t) :: tileDesc
    type(Grid_iterator_t) :: itor

    nullify(U)

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
                    Uold(1) = U(PHI0_VAR,i,j,k)
                    Uold(2) = U(PHI1_VAR,i,j,k)

                    ! U(PHI0_VAR,i,j,k) = U(PHI0_VAR,i,j,k)*fac
                    ! U(PHI1_VAR,i,j,k) = U(PHI1_VAR,i,j,k)*fac

                    call newton(root_func, jac_func, Uold, Unew, 1d-8, 1d-8, 100)

                    U(PHI0_VAR,i,j,k) = Unew(1)
                    U(PHI1_VAR,i,j,k) = Unew(2)
                end do ! i
            end do ! j
        end do ! k

        call tileDesc%releaseDataPtr(U, CENTER)

        call itor%next()
    end do ! itor

    call Grid_releaseTileIterator(itor)

contains

function root_func(x) result(y)
    implicit none
    
    real :: y(2)
    real, intent(in) :: x(2)

    real :: phi0, phi1

    phi0 = x(1); phi1 = x(2)

    y(1) = phi0 - Uold(1) + dt*sim_beta*phi0
    y(2) = phi1 - Uold(2) + dt*sim_beta*phi1
end function root_func

function jac_func(x) result(dy)
    implicit none
    
    real :: dy(2,2)
    real, intent(in) :: x(2)

    real :: phi0, phi1

    phi0 = x(1); phi1 = x(2)

    dy(1,1) = 1d0 + dt*sim_beta
    dy(1,2) = 0d0

    dy(2,1) = 0d0
    dy(2,2) = 1d0 + dt*sim_beta
end function jac_func

subroutine newton(func, jac, x0, root, rtol, atol, maxiter)
    implicit none

    interface
        function func(x) result(f)
            real :: f(2)
            real, intent(in) :: x(2)
        end function

        function jac(x) result(df)
            real :: df(2,2)
            real, intent(in) :: x(2)
        end function
    end interface

    real,    intent(in)  :: x0(2)
    real,    intent(out) :: root(2)
    real,    intent(in)  :: rtol
    real,    intent(in)  :: atol
    integer, intent(in)  :: maxiter

    real :: y(2), dy(2,2), idy(2,2), idetdy
    real :: x1(2), x2(2), dx(2), pivot(2)
    integer :: iter, info

    x1 = x0

    do iter = 1, maxiter
        y = func(x1)

        ! Check if we hit a zero exactly
        if (norm2(y) == 0d0) then
            root = x1
            return
        end if

        dy = jac(x1)

        ! Take a newton step
        ! idetdy = 1d0/(dy(1,1)*dy(2,2)-dy(1,2)*dy(2,1))
        ! idy(1,1) =  dy(2,2)*idetdy
        ! idy(1,2) = -dy(2,1)*idetdy
        ! idy(2,1) = -dy(1,2)*idetdy
        ! idy(2,2) =  dy(1,1)*idetdy

        ! dx(1) = idy(1,1)*y(1) + idy(1,2)*y(2)
        ! dx(2) = idy(2,1)*y(1) + idy(2,2)*y(2)

        dx = y
        call dgesv(2, 1, dy, 2, pivot, dx, 2, info)

        x2 = x1 - dx

        ! Check if we converged
        if (norm2(x2-x1) <= (atol + rtol*norm2(x1))) then
            root = x2
            return
        end if

        x1 = x2
    end do ! iter

    root = x1

end subroutine newton

end subroutine sim_molImplicitUpdate
