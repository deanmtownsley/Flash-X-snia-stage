!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_initBlock
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
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: blockDesc  )
!!
!!
!! DESCRIPTION
!!
!!   Initialize solution data in one block for a streaming sine wave
!!
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  blockDesc -        describes the block to initialize
!!
!! PARAMETERS
!!
!!  
!!***

!!REORDER(4): solnData


subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getCellCoords, Grid_getDeltas, &
       Grid_getGeometry

  use ProgramHeaderModule, ONLY : nE, nDOF, nDOFE, nNodesX, nNodesE
  use RadiationFieldsModule, ONLY : nSpecies, nCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  use ReferenceElementModuleX, ONLY: NodesX_q
  use UnitsModule, ONLY : Centimeter, Gram, Second

  implicit none

#include "constants.h"
#include "Simulation.h"
  
  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  logical, parameter :: useGuardCell = .TRUE.

  real, allocatable, dimension(:) :: xLeft, xCenter, xRight
  real, allocatable, dimension(:) :: yLeft, yCenter, yRight
  real, allocatable, dimension(:) :: zLeft, zCenter, zRight
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real, dimension(MDIM) :: delta
  real :: dx_coarse, dy_coarse, dz_coarse
  real :: xCenter_coarse, yCenter_coarse, zCenter_coarse

  integer :: ilo, ihi
  integer :: level

  integer :: meshGeom

  integer, dimension(1:MDIM) :: lo, hi, u_lo, u_hi
  integer :: i, j, k, n, ii, jj, kk, ic, jc, kc
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff, ivar
  integer :: nX(3)
  real :: xnode, ynode, znode, ss, ye

  real, parameter :: conv_x = Centimeter
  real, parameter :: conv_J = Gram/Second**2/Centimeter
  real, parameter :: conv_H = Gram/Second**3

  level     = tileDesc%level

  ! get dimensions/limits and coordinates
  lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
  hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)

  !! allocate all needed space
  allocate(xLeft  (lo(IAXIS):hi(IAXIS)) )
  allocate(xCenter(lo(IAXIS):hi(IAXIS)) )
  allocate(xRight (lo(IAXIS):hi(IAXIS)) )
  allocate(yLeft  (lo(JAXIS):hi(JAXIS)) )
  allocate(yCenter(lo(JAXIS):hi(JAXIS)) )
  allocate(yRight (lo(JAXIS):hi(JAXIS)) )
  allocate(zLeft  (lo(KAXIS):hi(KAXIS)) )
  allocate(zCenter(lo(KAXIS):hi(KAXIS)) )
  allocate(zRight (lo(KAXIS):hi(KAXIS)) )

  xLeft = 0.0 ; xCenter(:) = 0.0 ; xRight(:) = 0.0
  yLeft = 0.0 ; yCenter(:) = 0.0 ; yRight(:) = 0.0
  zLeft = 0.0 ; zCenter(:) = 0.0 ; zRight(:) = 0.0

  call Grid_getDeltas(level, delta)
  dx_coarse = delta(IAXIS) * THORNADO_NNODESX
  dy_coarse = delta(JAXIS) * THORNADO_NNODESX
  dz_coarse = delta(KAXIS) * THORNADO_NNODESX

  call Grid_getGeometry(meshGeom)

  call tileDesc%boundBox(boundBox)

  call Grid_getCellCoords(IAXIS,LEFT_EDGE, level,lo,hi,xLeft  )
  call Grid_getCellCoords(IAXIS,CENTER,    level,lo,hi,xCenter)
  call Grid_getCellCoords(IAXIS,RIGHT_EDGE,level,lo,hi,xRight )

  call Grid_getCellCoords(JAXIS,LEFT_EDGE, level,lo,hi,yLeft  )
  call Grid_getCellCoords(JAXIS,CENTER,    level,lo,hi,yCenter)
  call Grid_getCellCoords(JAXIS,RIGHT_EDGE,level,lo,hi,yRight )

  call Grid_getCellCoords(KAXIS,LEFT_EDGE, level,lo,hi,zLeft  )
  call Grid_getCellCoords(KAXIS,CENTER,    level,lo,hi,zCenter)
  call Grid_getCellCoords(KAXIS,RIGHT_EDGE,level,lo,hi,zRight )

  nX = hi - lo + 1
  nX(1:NDIM) = nX(1:NDIM) / 2

  u_lo = 1
  u_hi = nX

  do kc = u_lo(KAXIS), u_hi(KAXIS)
     do jc = u_lo(JAXIS), u_hi(JAXIS)
        do ic = u_lo(IAXIS), u_hi(IAXIS)
           i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
           j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
           k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

           xCenter_coarse = 0.5*(xLeft(i) + xRight(i+THORNADO_NNODESX-1))
           yCenter_coarse = 0.5*(yLeft(j) + yRight(j+THORNADO_NNODESX-1))
           zCenter_coarse = 0.5*(zLeft(k) + zRight(k+THORNADO_NNODESX-1))

           ! Initialize hydro data
           do iNodeX = 1, THORNADO_FLUID_NDOF
              ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
              jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
              kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k
              solnData(VELX_VAR,ii,jj,kk) = sim_velx_i
              solnData(VELY_VAR,ii,jj,kk) = sim_vely_i
              solnData(VELZ_VAR,ii,jj,kk) = sim_velz_i
              solnData(DENS_VAR,ii,jj,kk) = sim_dens_i
              solnData(TEMP_VAR,ii,jj,kk) = sim_temp_i
              solnData(PRES_VAR,ii,jj,kk) = sim_pres_i
              solnData(EINT_VAR,ii,jj,kk) = sim_eint_i
              solnData(ENER_VAR,ii,jj,kk) = sim_etot_i
              solnData(GAMC_VAR,ii,jj,kk) = sim_gamc_i
              solnData(GAME_VAR,ii,jj,kk) = sim_game_i
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,ii,jj,kk) = sim_xn_i(n)
              enddo
              solnData(YE_MSCALAR,ii,jj,kk) = sim_ye_i
           enddo

           ! Initialize neutrino data
           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS ; do iE = 1, THORNADO_NE

              ioff = THORNADO_BEGIN &
                 + (iS -1)*(THORNADO_NNODESE*THORNADO_NE*THORNADO_NMOMENTS) &
                 + (iCR-1)*(THORNADO_NNODESE*THORNADO_NE) &
                 + (iE -1)*(THORNADO_NNODESE)

              do iNode = 1, THORNADO_RAD_NDOF

                 iNodeE = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                 jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                 kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                 ! calculate the indices
                 ivar = ioff + iNodeE - 1

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = xCenter_coarse + NodesX_q(1,iNodeX)
                 ynode = yCenter_coarse + NodesX_q(2,iNodeX)
                 znode = zCenter_coarse + NodesX_q(3,iNodeX)

                 ss = 1.0e0 + sin(2.0*PI*xnode * conv_x)

                 ! J moment, iCR = 1
                 if (iCR == iCR_N) solnData(ivar,ii,jj,kk) = ss! / conv_J

                 ! H_x moment, iCR = 2
                 if (iCR == iCR_G1) solnData(ivar,ii,jj,kk) = ss! / conv_H

                 ! H_y moment, iCR = 3
                 if (iCR == iCR_G2) solnData(ivar,ii,jj,kk) = 0.0e0

                 ! H_z moment, iCR = 4
                 if (iCR == iCR_G3) solnData(ivar,ii,jj,kk) = 0.0e0

              enddo
           enddo ; enddo ; enddo

        enddo
     enddo
  enddo

  ! cleanup
  deallocate(xLeft,xCenter,xRight)
  deallocate(yLeft,yCenter,yRight)
  deallocate(zLeft,zCenter,zRight)

  return
end subroutine Simulation_initBlock
