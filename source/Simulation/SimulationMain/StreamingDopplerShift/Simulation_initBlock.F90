!!****if* source/Simulation/SimulationMain/StreamingDopplerShift/Simulation_initBlock
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
!!   Initialize solution data in one block for streaming doppler wave test
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
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_coordTransfm
  use Eos_interface, ONLY : Eos_multiDim

  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use KindModule, ONLY : DP, Zero, One, Three, TwoPi
  use UnitsModule, ONLY : Centimeter, Gram, Second, MeV, Kelvin
  use MeshModule, ONLY : NodeCoordinate, MeshE, MeshX
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch

  use rt_data, ONLY : rt_UpperBry1

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
  
  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  real, allocatable, dimension(:) :: xCenter
  real, allocatable, dimension(:) :: yCenter
  real, allocatable, dimension(:) :: zCenter
  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: rcc, tcc, pcc, velr, velt, velp

  integer, dimension(1:MDIM)   :: lo, hi
  integer, dimension(1:MDIM+1) :: u_lo, u_hi
  integer :: i, j, k, n, ii, jj, kk, ic, jc, kc
  integer :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff, ivar
  integer :: nX(3), swX(3) ! Element # in Thornado DG grid
  real :: xL(3), xR(3)
  real :: xnode, ynode, znode, enode, Psi0, Psi1
  real :: dens_buffer, temp_buffer, ye_buffer, velr_buffer
  real :: vel_buffer
  integer :: faces(LOW:HIGH, 1:MDIM)

  real, parameter :: conv_x = Centimeter
  real, parameter :: UnitD  = Gram / Centimeter**3
  real, parameter :: UnitT  = Kelvin
  real, parameter :: UnitY  = 1.0
  real, parameter :: UnitV  = Centimeter/Second
  real, parameter :: conv_J = Gram/Second**2/Centimeter
  real, parameter :: conv_H = Gram/Second**3
  real, parameter :: conv_e = MeV

  ! get boundary information
  call tileDesc%faceBCs(faces)

  ! get dimensions/limits and coordinates
  lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
  hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)

  !! allocate all needed space
  allocate(xCenter(lo(IAXIS):hi(IAXIS)))
  allocate(yCenter(lo(JAXIS):hi(JAXIS)))
  allocate(zCenter(lo(KAXIS):hi(KAXIS)))

  xCenter = 0.0
  yCenter = 0.0
  zCenter = 0.0

  call Grid_getCellCoords(IAXIS,CENTER,tileDesc%level,lo,hi,xCenter)
  if ( NDIM > 1 ) call Grid_getCellCoords(JAXIS,CENTER,tileDesc%level,lo,hi,yCenter)
  if ( NDIM > 2 ) call Grid_getCellCoords(KAXIS,CENTER,tileDesc%level,lo,hi,zCenter)

  nX = 1
  swX = 0
  xL = 0.0
  if ( sim_geometry == CARTESIAN ) then
     xR = 1.0
  else if ( sim_geometry == CYLINDRICAL ) then
     xR = [ 1.0, 1.0, 2.0*PI ]
  else if ( sim_geometry == SPHERICAL ) then
     xR = [ 1.0, PI, 2.0*PI ]
  else
     call Driver_abort("Geometry not supported")
  end if

  nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
  swX(1:NDIM) = 2
  u_lo(2:4) = 1 - swX
  u_hi(2:4) = nX + swX
  u_lo(1)   = 1 - THORNADO_SWE
  u_hi(1)   = THORNADO_NE + THORNADO_SWE

  call tileDesc%boundBox(boundBox)
  xL(1:NDIM) = boundBox(LOW, 1:NDIM)
  xR(1:NDIM) = boundBox(HIGH,1:NDIM)

  ! convert cm to m for Thornado
  xL(1) = xL(1) * conv_x
  xR(1) = xR(1) * conv_x
  if ( sim_geometry /= SPHERICAL ) then
     xL(2) = xL(2) * conv_x
     xR(2) = xR(2) * conv_x
  end if
  if ( sim_geometry == CARTESIAN ) then
     xL(3) = xL(3) * conv_x
     xR(3) = xR(3) * conv_x
  end if

  do iX3 = 1, nX(3)
     do iX2 = 1, nX(2)
        do iX1 = 1, nX(1)

           i = lo(IAXIS) + THORNADO_NNODESX*(iX1-1)
           j = lo(JAXIS) + THORNADO_NNODESX*(iX2-1)
           k = lo(KAXIS) + THORNADO_NNODESX*(iX3-1)

           ! Initialize hydro data
           do iNodeX = 1, THORNADO_FLUID_NDOF

              ii = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
              jj = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
              kk = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

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

              select case( sim_rad_direction )
              case( 'X' )

                solnData(VELX_VAR,ii,jj,kk) &
                  = InitialDopplerVelocity( xCenter(ii) , sim_velx_i )

                solnData(VELY_VAR,ii,jj,kk) = sim_vely_i
                solnData(VELZ_VAR,ii,jj,kk) = sim_velz_i

              case( 'Y' )

                solnData(VELY_VAR,ii,jj,kk) &
                  = InitialDopplerVelocity( yCenter(jj) , sim_vely_i )

                solnData(VELX_VAR,ii,jj,kk) = sim_velx_i
                solnData(VELZ_VAR,ii,jj,kk) = sim_velz_i

              case( 'Z' )

                solnData(VELZ_VAR,ii,jj,kk) &
                  = InitialDopplerVelocity( zCenter(kk) , sim_velz_i )

                solnData(VELX_VAR,ii,jj,kk) = sim_velx_i
                solnData(VELY_VAR,ii,jj,kk) = sim_vely_i

              end select

           end do

           ! Initialize neutrino data
           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS ; do iE = u_lo(1), u_hi(1)

              ioff = THORNADO_BEGIN &
                 + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE) &
                           *THORNADO_NMOMENTS) &
                 + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                 + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

              do iNode = 1, THORNADO_RAD_NDOF

                 iNodeE  = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX  = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 iNodeX1 = mod((iNodeX-1)                    ,THORNADO_NNODESX) + 1
                 ii      = iNodeX1 + i - 1

                 iNodeX2 = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + 1
                 jj      = iNodeX2 + j - 1

                 iNodeX3 = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + 1
                 kk      = iNodeX3 + k - 1

                 ! calculate the indices
                 ivar = ioff + iNodeE - 1

                 ! J moment, iCR = 1
                 if (iCR == iCR_N) solnData(ivar,ii,jj,kk)  = 1.0e-40
            
                 ! H_x moment, iCR = 2
                 if (iCR == iCR_G1) solnData(ivar,ii,jj,kk) = 0.0e0

                 ! H_y moment, iCR = 3
                 if (iCR == iCR_G2) solnData(ivar,ii,jj,kk) = 0.0e0

                 ! H_z moment, iCR = 4
                 if (iCR == iCR_G3) solnData(ivar,ii,jj,kk) = 0.0e0

              end do
           end do ; end do ; end do

        end do
     end do
  end do

  ! cleanup
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)

  return

contains

  pure real function InitialDopplerVelocity( X1, V_0 )

    real, intent(in) :: X1, V_0

    if( X1 .lt. X_0 )then
      InitialDopplerVelocity &
        = 0.0d0
    elseif( X1 .ge. X_0 .and. X1 .lt. X_1 )then
      InitialDopplerVelocity &
        = V_0 * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
    elseif( X1 .ge. X_1 .and. X1 .lt. X_2 )then
      InitialDopplerVelocity &
        = V_0
    elseif( X1 .ge. X_2 .and. X1 .lt. X_3 )then
      InitialDopplerVelocity &
        = V_0 * SIN( TwoPi * ( X1 - X_0 ) / L_X )**2
    else
      InitialDopplerVelocity &
        = 0.0d0
    end if

    return

  end function InitialDopplerVelocity

end subroutine Simulation_initBlock
