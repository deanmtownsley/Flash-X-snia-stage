!!****if* source/Simulation/SimulationMain/DiffusionSineWave/Simulation_initBlock
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
!!   Initialize solution data in one block for diffusion sine wave
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

#include "constants.h"
#include "Simulation.h"

subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getGeometry

  use KindModule, ONLY : TwoPi
  use GeometryFieldsModule, ONLY : uGF, iGF_Gm_dd_11, iGF_Gm_dd_22, iGF_Gm_dd_33
  use MeshModule, ONLY : NodeCoordinate, MeshX
  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use UnitsModule, ONLY : Centimeter, Second

#if defined(THORNADO_ORDER_V)
  use TwoMoment_UtilitiesModule, ONLY : ComputeConserved_TwoMoment
#endif

  implicit none
  
  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  real, dimension(LOW:HIGH,MDIM) :: boundBox

  integer :: meshGeom

  integer, dimension(1:MDIM) :: lo, hi, u_lo, u_hi
  integer :: i, j, k, n, ii, jj, kk
  integer :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff, ivar
  integer :: nX(3), swX(3)
  real :: xL(3), xR(3), L_x, Sigma
  real :: xnode, ynode, znode, ss_D, ss_I
  real :: Nnu, Gnu1, Gnu2, Gnu3
  real :: Dnu, Inu1, Inu2, Inu3

  real, parameter :: conv_x = Centimeter
  real, parameter :: UnitV  = Centimeter / Second
  real, parameter :: A = 0.5, B = 0.49

  ! get dimensions/limits and coordinates
  lo(1:MDIM) = tileDesc % limits(LOW ,1:MDIM)
  hi(1:MDIM) = tileDesc % limits(HIGH,1:MDIM)

  call Grid_getGeometry(meshGeom)
  nX  = 1
  swX = 0
  xL  = 0.0
  if ( meshGeom == CARTESIAN ) then
     xR = 1.0
  else
     call Driver_abort("Geometry not supported")
  end if

  nX (1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
  swX(1:NDIM) = NGUARD / THORNADO_NNODESX
  u_lo = 1  - swX
  u_hi = nX + swX

  call tileDesc%boundBox(boundBox)
  xL(1:NDIM) = boundBox(LOW, 1:NDIM)
  xR(1:NDIM) = boundBox(HIGH,1:NDIM)

  ! convert cm to m for Thornado (cartesian geometry assumed)
  xL = xL * conv_x
  xR = xR * conv_x

  call InitThornado_Patch &
        ( nX, swX, xL, xR, THORNADO_NSPECIES, 'cartesian' )

  L_x = ( sim_xmax - sim_xmin ) * conv_x

  Sigma = sim_Sigma / Centimeter

  do iX3 = u_lo(3), u_hi(3)
     do iX2 = u_lo(2), u_hi(2)
        do iX1 = u_lo(1), u_hi(1)

           i = lo(IAXIS) + THORNADO_NNODESX*(iX1-1)
           j = lo(JAXIS) + THORNADO_NNODESX*(iX2-1)
           k = lo(KAXIS) + THORNADO_NNODESX*(iX3-1)

           ! Initialize hydro data
           do iNodeX = 1, THORNADO_FLUID_NDOF

              ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
              jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
              kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

              solnData(DENS_VAR,ii,jj,kk) = sim_dens_i
              solnData(VELX_VAR,ii,jj,kk) = sim_velx_i
              solnData(VELY_VAR,ii,jj,kk) = sim_vely_i
              solnData(VELZ_VAR,ii,jj,kk) = sim_velz_i
              solnData(TEMP_VAR,ii,jj,kk) = sim_temp_i
              solnData(PRES_VAR,ii,jj,kk) = sim_pres_i

#if NSPECIES > 0
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,ii,jj,kk) = sim_xn_i(n)
              enddo
#endif
              solnData(YE_MSCALAR,ii,jj,kk) = sim_ye_i
           enddo

           ! Initialize neutrino data
           do iS = 1, THORNADO_NSPECIES ; do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE
              do iNode = 1, THORNADO_RAD_NDOF

                 ! calculate the indices
                 iNodeE  = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX  = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 iNodeX1 = mod((iNodeX-1)                    ,THORNADO_NNODESX) + 1
                 iNodeX2 = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + 1
                 iNodeX3 = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + 1

                 ii      = iNodeX1 + i - 1
                 jj      = iNodeX2 + j - 1
                 kk      = iNodeX3 + k - 1

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                 ynode = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
                 znode = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

                 ss_D = A + B * sin( TwoPi * xnode / L_x )
                 ss_I = - B * ( TwoPi / L_x ) / ( 3.0 * Sigma ) &
                            * cos( TwoPi * xnode / L_x )

#if   defined(THORNADO_ORDER_1)
                 Nnu  = ss_D
                 Gnu1 = ss_I
                 Gnu2 = 0.0
                 Gnu3 = 0.0
#elif defined(THORNADO_ORDER_V)
                 Dnu  = ss_D
                 Inu1 = ss_I
                 Inu2 = 0.0
                 Inu3 = 0.0
                 CALL ComputeConserved_TwoMoment &
                    ( Dnu, Inu1, Inu2, Inu3, &
                      Nnu, Gnu1, Gnu2, Gnu3, &
                      sim_velx_i * UnitV, &
                      sim_vely_i * UnitV, &
                      sim_velz_i * UnitV, &
                      uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_11), &
                      uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_22), &
                      uGF(iNodeX,iX1,iX2,iX3,iGF_Gm_dd_33) )
#endif

                 ! J moment, iCR = 1
                 iCR  = iCR_N
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Nnu

                 ! H_x moment, iCR = 2
                 iCR  = iCR_G1
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu1

                 ! H_y moment, iCR = 3
                 iCR  = iCR_G2
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu2

                 ! H_z moment, iCR = 4
                 iCR  = iCR_G3
                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 ivar = ioff + iNodeE - 1
                 solnData(ivar,ii,jj,kk) = Gnu3

              enddo
           enddo ; enddo

        enddo
     enddo
  enddo

  ! cleanup
  call FreeThornado_Patch()

  return
end subroutine Simulation_initBlock
