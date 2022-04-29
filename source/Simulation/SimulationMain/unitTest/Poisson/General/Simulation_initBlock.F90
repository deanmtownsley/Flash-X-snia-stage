!!****if* source/Simulation/SimulationMain/unitTest/Poisson/General/Simulation_initBlock
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
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(in) :: blockID)
!!
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified tile.
!!
!!  Reference:
!!
!!
!! ARGUMENTS
!!
!!  tile -          the tile to update
!!
!!
!!
!!
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(solnData, tileDesc)

!  use Simulation_data
   use Simulation_data, ONLY: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax
   use Grid_interface, ONLY: Grid_getCellCoords, Grid_getDomainBC
   use Grid_tile, ONLY: Grid_tile_t

   implicit none

#include "constants.h"
#include "Simulation.h"

   !Arguments-----------------------
   real, dimension(:, :, :, :), pointer :: solnData
   type(Grid_tile_t), intent(in) :: tileDesc
   integer :: tileDescID

   integer :: i, j, k
   integer, dimension(MDIM) :: lo, hi
   real, allocatable, dimension(:) ::xCenter, yCenter, zCenter
   real :: Lx, Ly, Lz, xi, yi, zi, Phi_ijk, F_ijk
   logical :: gcell = .true.
   integer :: domainBC(LOW:HIGH, 1:MDIM)

   !----------------------------------------------------------------------
   call Grid_getDomainBC(domainBC)

   lo = tileDesc%blkLimitsGC(LOW, :)
   hi = tileDesc%blkLimitsGC(HIGH, :)

   allocate (xCenter(lo(IAXIS):hi(IAXIS)))
   allocate (yCenter(lo(JAXIS):hi(JAXIS)))
   allocate (zCenter(lo(KAXIS):hi(KAXIS)))

   xCenter = 0.0
   yCenter = 0.0
   zCenter = 0.0

   call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, lo, hi, xCenter)
   if (NDIM >= 2) call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, lo, hi, yCenter)
   if (NDIM == 3) call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, lo, hi, zCenter)

#ifdef DEBUG_SIMULATION
98 format('initTile:', A4, '(', I3, ':   ,', I3, ':   ,', I3, ':   ,', I3, ':   )')
99 format('initTile:', A4, '(', I3, ':', I3, ',', I3, ':', I3, ',', I3, ':', I3, ',', I3, ':', I3, ')')
   print 99, "solnData", (lbound(solnData, i), ubound(solnData, i), i=1, 4)
   print *, 'blkLim  :', lo, hi
#endif
!------------------------------------------------------------------------------

   Lx = sim_xMax - sim_xMin
   Ly = sim_yMax - sim_yMin
   Lz = sim_zMax - sim_zMin

   do k = lo(KAXIS), hi(KAXIS)
      do j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS)
            xi = xCenter(i)
            yi = yCenter(j)
            zi = zCenter(k)

#if NDIM == 2
            if (domainBC(LOW, IAXIS) == DIRICHLET) then
               Phi_ijk = sin(2.*PI*xi)*sin(2.*PI*yi)
            else
               Phi_ijk = cos(2.*PI*xi)*cos(2.*PI*yi)
            end if

            F_ijk = -8*PI*PI*Phi_ijk

#else
            if (domainBC(LOW, IAXIS) == DIRICHLET) then
               Phi_ijk = sin(2.*PI*xi)*sin(2.*PI*yi)*sin(2.*PI*zi)
            else
               Phi_ijk = cos(2.*PI*xi)*cos(2.*PI*yi)*cos(2.*PI*zi)
            end if

            F_ijk = -12*PI*PI*Phi_ijk
#endif

            solnData(ASOL_VAR, i, j, k) = Phi_ijk
            solnData(RHDS_VAR, i, j, k) = F_ijk
            solnData(RFIN_VAR, i, j, k) = 1 - abs(Phi_ijk)

         end do
      end do
   end do

   ! set values for other variables
   solnData(DIFF_VAR, lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) = 0.0
   solnData(NSOL_VAR, lo(IAXIS):hi(IAXIS), lo(JAXIS):hi(JAXIS), lo(KAXIS):hi(KAXIS)) = 0.0

   deallocate (xCenter, yCenter, zCenter)

   return

end subroutine Simulation_initBlock
