!!****if* source/Simulation/SimulationMain/Brusselator/Simulation_molImplicitUpdate
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
!!      Simulation_molImplicitUpdate
!!
!!  SYNOPSIS
!!
!!      call Simulation_molImplicitUpdate(real, intent(in) :: t
!!                                        real, intent(in) :: dt)
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

!!REORDER(4): vars

subroutine Simulation_molImplicitUpdate(t, dt)
   use Simulation_data, only: sim_alpha

   use Grid_interface, only: Grid_getTileIterator, Grid_releaseTileIterator
   use Grid_iterator, only: Grid_iterator_t
   use Grid_tile, only: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none

   real, intent(in) :: t, dt

   real :: del(MDIM), dx2, fac
   integer, dimension(LOW:HIGH, MDIM) :: lim, bcs
   integer :: i, j, k, ilo, ihi

   real, dimension(:, :, :, :), pointer :: vars
   real :: u0, v0, w0

   real, dimension(:), allocatable :: DL, D, DU
   real, dimension(:, :), allocatable :: B
   integer :: N, info

   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

   ! LAPACK tri-diagonal matrix equation solver
   ! https://netlib.org/lapack/explore-html/d4/d62/group__double_g_tsolve_ga2bf93f2ddefa5e671866eb2191dc19d4.html
   external :: dgtsv

   nullify (vars)

   call Grid_getTileIterator(itor, LEAF)

   do ! not using while since this is technically deprecated
      if (.not. itor%isValid()) exit

      call itor%currentTile(tileDesc)

      call tileDesc%faceBCs(bcs)

      lim = tileDesc%limits
      ilo = lim(LOW, IAXIS)
      ihi = lim(HIGH, IAXIS)

      if (bcs(LOW, IAXIS) .ne. NOT_BOUNDARY) ilo = ilo + 1
      if (bcs(HIGH, IAXIS) .ne. NOT_BOUNDARY) ihi = ihi - 1

      call tileDesc%deltas(del)
      dx2 = del(IAXIS)**2
      fac = (sim_alpha*dt)/dx2

      N = ihi - ilo + 1
      allocate (D(N)); allocate (B(N, 3))
      allocate (DL(N - 1)); allocate (DU(N - 1))

      ! This is equivalent to `call MoL_getDataPtr(tileDesc, vars, MOL_EVOLVED)`
      ! that is utilized in the RHS procedures of this simulation
      call tileDesc%getDataPtr(vars, CENTER)

      do k = lim(LOW, KAXIS), lim(HIGH, KAXIS)
         do j = lim(LOW, JAXIS), lim(HIGH, JAXIS)
            ! Solve for U,V,W
            ! Note: DGTSV overwrites D,DL,DU so these need to be set everytime
            DU = -fac
            DL = -fac

            D = 1.0 + 2.0*fac

            B(:, 1) = vars(U_VAR, ilo:ihi, j, k)
            B(:, 2) = vars(V_VAR, ilo:ihi, j, k)
            B(:, 3) = vars(W_VAR, ilo:ihi, j, k)

            B(1, 1) = B(1, 1) + fac*vars(U_VAR, ilo - 1, j, k)
            B(1, 2) = B(1, 2) + fac*vars(V_VAR, ilo - 1, j, k)
            B(1, 3) = B(1, 3) + fac*vars(W_VAR, ilo - 1, j, k)

            B(N, 1) = B(N, 1) + fac*vars(U_VAR, ihi + 1, j, k)
            B(N, 2) = B(N, 2) + fac*vars(V_VAR, ihi + 1, j, k)
            B(N, 3) = B(N, 3) + fac*vars(W_VAR, ihi + 1, j, k)

            call dgtsv(N, 3, DL, D, DU, B, N, info)

            vars(U_VAR, ilo:ihi, j, k) = B(:, 1)
            vars(V_VAR, ilo:ihi, j, k) = B(:, 2)
            vars(W_VAR, ilo:ihi, j, k) = B(:, 3)
         end do ! j
      end do ! k

      call tileDesc%releaseDataPtr(vars, CENTER)

      deallocate (D); deallocate (DL); deallocate (DU); deallocate (B)

      call itor%next()
   end do ! itor

   call Grid_releaseTileIterator(itor)

   call Simulation_molPostUpdate(t)

end subroutine Simulation_molImplicitUpdate
