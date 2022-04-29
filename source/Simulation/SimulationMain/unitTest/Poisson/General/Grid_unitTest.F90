!!****if* source/Simulation/SimulationMain/unitTest/Poisson/General/Grid_unitTest
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
!! NAME
!!
!!  Grid_unitTest
!!
!! SYNOPSIS
!!
!!  call Grid_unitTest(integer(in):: fileUnit,
!!                     logical(inout)::perfect  )
!!
!! DESCRIPTION
!!
!!  This unit test is to test the multigrid solver of AMReX library.
!!  It uses known analytical function as well as harmonic
!!
!! ARGUMENTS
!!
!!  fileUnit - open f90 write unit
!!  perfect - returns a true if the test passed, false otherwise
!!
!! NOTES
!!
!!***

!!REORDER(4): solnData

subroutine Grid_unitTest(fileUnit, perfect)

   use Grid_interface, ONLY: GRID_PDE_BND_PERIODIC, GRID_PDE_BND_DIRICHLET, GRID_PDE_BND_NEUMANN, &
                             Grid_solvePoisson, &
                             Grid_getDeltas, Grid_fillGuardCells, &
                             Grid_getTileIterator, Grid_releaseTileIterator

   use gr_interface, ONLY: gr_findMean
   use Grid_data, ONLY: gr_meshMe, gr_meshComm, gr_domainBC
   use Grid_iterator, ONLY: Grid_iterator_t
   use Grid_tile, ONLY: Grid_tile_t

#include "Simulation.h"
#include "constants.h"

   implicit none
#include "Flashx_mpi.h"

   integer, intent(in)           :: fileUnit ! Output to file
   logical, intent(inout)        :: perfect  ! Flag to indicate errors

   real, dimension(:, :, :, :), pointer :: solnData
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc

   integer, dimension(LOW:HIGH, MDIM) :: tileLimits
   integer, dimension(2*MDIM) :: bcTypes
   real, dimension(2, 2*MDIM) :: bcValues

   logical :: gcMask(NUNK_VARS), evaluateData

   real :: poisfact
   integer :: blkCount = 0, lb, i, j, k
   real:: del(MDIM)
   real meanASOL, meanPFFT
   integer nx, ny, nz
   real, allocatable, dimension(:, :, :) :: fsrc, usol
   integer blkpoints, blkpointsaux, blkCountaux
   real L2_err, L2_erraux, Linf_err, Linf_erraux, Tvol, Tvolaux, vcell
   integer :: refinelevel
   real, parameter :: tol = 1.e-10
   real, parameter :: tolInf = 4.1e-3
   real, parameter :: tol2 = 2.24e-3

   integer TA(2), count_rate, ierr
   real :: ET

   ! -------------------------------------------------------------------
   nullify (solnData)

   if (gr_domainBC(LOW, IAXIS) == DIRICHLET) then
      bcTypes(:) = GRID_PDE_BND_DIRICHLET
   else
      bcTypes(:) = GRID_PDE_BND_NEUMANN
   end if

   bcValues(:, :) = 0.

   call mpi_barrier(gr_meshComm, ierr)
   if (gr_meshMe .eq. 0) CALL SYSTEM_CLOCK(TA(1), count_rate)
   poisfact = 1.
   call Grid_solvePoisson(NSOL_VAR, RHDS_VAR, bcTypes, bcValues, poisfact)
   call mpi_barrier(gr_meshComm, ierr)
   if (gr_meshMe .eq. 0) then
      CALL SYSTEM_CLOCK(TA(2), count_rate)
      ET = REAL(TA(2) - TA(1))/count_rate
      write (*, *) ' '
      write (*, *) '3 PERIODIC Poisson Solver time = ', ET, ' sec.'
   end if

   gcMask = .false.
   gcMask(NSOL_VAR) = .true.

   call Grid_fillGuardCells(CENTER, ALLDIR, masksize=NUNK_VARS, mask=gcMask)

   ! Check error in the solution:
   L2_err = 0.
   blkpoints = 0.
   Linf_err = 0.
   Tvol = 0.
   ! Get Tile iterator
   call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
   do while (itor%isValid())
      call itor%currentTile(tileDesc)
      !get the index limits of the tile
      tileLimits = tileDesc%limits

      ! get a pointer to the current tile of data
      call tileDesc%getDataPtr(solnData, CENTER)
      call tileDesc%deltas(del)

      select case (NDIM)
      case (1)
         vcell = del(IAXIS)
      case (2)
         vcell = del(IAXIS)*del(JAXIS)
      case (3)
         vcell = del(IAXIS)*del(JAXIS)*del(KAXIS)
      end select

      blkCount = blkCount + 1
      blkpoints = blkpoints + &
                  (tileLimits(HIGH, IAXIS) - tileLimits(LOW, IAXIS) + 1)* &
                  (tileLimits(HIGH, JAXIS) - tileLimits(LOW, JAXIS) + 1)* &
                  (tileLimits(HIGH, KAXIS) - tileLimits(LOW, KAXIS) + 1)

      Tvol = Tvol + vcell*real((tileLimits(HIGH, IAXIS) - tileLimits(LOW, IAXIS) + 1)* &
                               (tileLimits(HIGH, JAXIS) - tileLimits(LOW, JAXIS) + 1)* &
                               (tileLimits(HIGH, KAXIS) - tileLimits(LOW, KAXIS) + 1))

      solnData(DIFF_VAR, tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS), &
               tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS), &
               tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS)) = &
         abs(solnData(NSOL_VAR, tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS), &
                      tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS), &
                      tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS)) - &
             solnData(ASOL_VAR, tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS), &
                      tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS), &
                      tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS)))

      ! L2 norm of error:
      L2_err = L2_err + sum(vcell*solnData(DIFF_VAR, tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS), &
                                           tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS), &
                                           tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS))**2.)

      ! Linf norm of error:
      Linf_err = max(Linf_err, maxval(solnData(DIFF_VAR, tileLimits(LOW, IAXIS):tileLimits(HIGH, IAXIS), &
                                               tileLimits(LOW, JAXIS):tileLimits(HIGH, JAXIS), &
                                               tileLimits(LOW, KAXIS):tileLimits(HIGH, KAXIS))))

      call tileDesc%releaseDataPtr(solnData, CENTER)
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)
! #if defined(__GFORTRAN__) && (__GNUC__ <= 4)
!   call destroy_iterator(itor)
! #endif

   ! Sum processors Volumes
   Tvolaux = Tvol
   call MPI_Allreduce(Tvolaux, Tvol, 1, FLASH_REAL, &
                      MPI_SUM, gr_meshComm, ierr)

   ! Sum processors points
   blkpointsaux = blkpoints
   call MPI_Allreduce(blkpointsaux, blkpoints, 1, FLASH_INTEGER, &
                      MPI_SUM, gr_meshComm, ierr)

   ! Sum processors L2 norm of error squared
   L2_erraux = L2_err
   call MPI_Allreduce(L2_erraux, L2_err, 1, FLASH_REAL, &
                      MPI_SUM, gr_meshComm, ierr)

   ! Compute L2 norm for whole domain
   L2_err = sqrt(L2_err/Tvol)

   ! Sum processors Linf norm of error
   Linf_erraux = Linf_err
   call MPI_Allreduce(Linf_erraux, Linf_err, 1, FLASH_REAL, &
                      MPI_MAX, gr_meshComm, ierr)

   ! Fill GuardCells:
   gcMask = .TRUE.
   call Grid_fillGuardCells(CENTER, ALLDIR, &
                            maskSize=NUNK_VARS, mask=gcMask)

   ! Export to Tecplot:
   !call outtotecplot(gr_meshMe,0.0,1.,1,0,0.0,blkList,blkCount,0)

   call gr_findMean(ASOL_VAR, 2, .false., meanASOL)
   call gr_findMean(NSOL_VAR, 2, .false., meanPFFT)

   !Unit test gives a mean analytical solution of zero.  Ensure the absolute
   !value of the mean numerical solution is between 0 and the tolerance value.
   perfect = abs(meanPFFT) < MAX(TINY(1.), tol)

   if (perfect) perfect = (abs(Linf_err) .LE. tolInf)
   if (perfect) perfect = (L2_err .LE. tol2)

   ! Sum processors points
   call MPI_Allreduce(blkCount, blkCountaux, 1, FLASH_INTEGER, &
                      MPI_SUM, gr_meshComm, ierr)

   if (gr_meshMe .eq. 0) then
      write (*, *) ' '
      write (*, '(A,2g16.8)') ' Mean Analytical, Numerical Sol=', meanASOL, meanPFFT
      write (*, '(A,1g16.8)') " ||Phi - PhiAnalytical||inf =", Linf_err
      write (*, '(A,1g16.8)') " ||Phi - PhiAnalytical||2   =", L2_err
      write (*, *) " Total Volume =", Tvol
      write (*, *) " Total Number of Tiles on Leaf Blocks=", blkCountaux
      if (perfect) print *, "All tests PASSED"
      if (.not. perfect) print *, "Test FAILED"
      write (*, *) ' '
   end if
   return

end subroutine Grid_unitTest
