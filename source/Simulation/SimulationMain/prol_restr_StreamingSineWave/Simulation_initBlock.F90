!!****if* source/Simulation/SimulationMain/StreamingSineWave/Simulation_initBlock
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
  use Eos_interface, ONLY : Eos
  use Grid_tile, ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getGeometry

  use KindModule, ONLY : TwoPi
  use MeshModule, ONLY : NodeCoordinate, MeshX
  use RadiationFieldsModule, ONLY : iCR_N, iCR_G1, iCR_G2, iCR_G3
  use ThornadoInitializationModule, ONLY : InitThornado_Patch, FreeThornado_Patch
  use UnitsModule, ONLY : Centimeter, Gram, Second

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"
  
  real, dimension(:,:,:,:), pointer :: solnData
  type(Grid_tile_t), intent(in)     :: tileDesc

  real, dimension(LOW:HIGH,MDIM) :: boundBox

  real, dimension(SPECIES_BEGIN:SPECIES_END) :: massFraction
  real, dimension(EOS_NUM) :: eosData

  integer :: meshGeom

  integer, dimension(1:MDIM) :: lo, hi, u_lo, u_hi
  integer :: i, j, k, n, ii, jj, kk
  integer :: iX1, iX2, iX3, iNodeX1, iNodeX2, iNodeX3
  integer :: iS, iCR, iE, iNode, iNodeX, iNodeE, ioff, ivar
  integer :: nX(3), swX(3)
  real :: xL(3), xR(3)
  real :: xnode, ynode, znode, ss, ye

  real, parameter :: conv_x = Centimeter
  real, parameter :: conv_J = Gram/Second**2/Centimeter
  real, parameter :: conv_H = Gram/Second**3

  ! get dimensions/limits and coordinates
  lo(1:MDIM) = tileDesc%limits(LOW,1:MDIM)
  hi(1:MDIM) = tileDesc%limits(HIGH,1:MDIM)

  call Grid_getGeometry(meshGeom)
  nX = 1
  swX = 0
  xL = 0.0
  if ( meshGeom == CARTESIAN ) then
     xR = 1.0
  else
     call Driver_abort("Geometry not supported")
  end if

  nX(1:NDIM) = (hi(1:NDIM) - lo(1:NDIM) + 1) / THORNADO_NNODESX
  swX(1:NDIM) = 2
  u_lo = 1 - swX
  u_hi = nX + swX

  call tileDesc%boundBox(boundBox)
  xL(1:NDIM) = boundBox(LOW, 1:NDIM)
  xR(1:NDIM) = boundBox(HIGH,1:NDIM)

  ! convert cm to m for Thornado (cartesian geometry assumed)
  xL = xL * conv_x
  xR = xR * conv_x

  call InitThornado_Patch &
       (nX, swX, xL, xR, THORNADO_NSPECIES, 'cartesian' )

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

              if ( MeshX(1) % Center(iX1) / conv_x < 0.5 * ( sim_xmin + sim_xmax ) ) then
                solnData(DENS_VAR,ii,jj,kk) = sim_dens_lo_i
              else
                solnData(DENS_VAR,ii,jj,kk) = sim_dens_hi_i
              end if
              solnData(VELX_VAR,ii,jj,kk) = sim_velx_i
              solnData(VELY_VAR,ii,jj,kk) = sim_vely_i
              solnData(VELZ_VAR,ii,jj,kk) = sim_velz_i
              solnData(TEMP_VAR,ii,jj,kk) = sim_temp_i
              solnData(PRES_VAR,ii,jj,kk) = sim_pres_i
              do n = SPECIES_BEGIN,SPECIES_END
                 solnData(n,ii,jj,kk) = sim_xn_i(n)
              enddo
              solnData(YE_MSCALAR,ii,jj,kk) = sim_ye_i

              ! Let the inital Eos_wrapped that follows Simulation_initBlock do this

              !eosData(EOS_DENS) = solnData(DENS_VAR,ii,jj,kk)
              !eosData(EOS_TEMP) = solnData(TEMP_VAR,ii,jj,kk)
              !massFraction = sim_xn_i

              !call Eos(MODE_DENS_TEMP,1,eosData,massFraction)

              !solnData(PRES_VAR,ii,jj,kk) = eosData(EOS_PRES)
              !solnData(EINT_VAR,ii,jj,kk) = eosData(EOS_EINT)
              !solnData(ENER_VAR,ii,jj,kk) = eosData(EOS_EINT) + 0.5*(sim_velx_i**2 + sim_vely_i**2 + sim_velz_i**2)
              !solnData(GAMC_VAR,ii,jj,kk) = eosData(EOS_GAMC)
              !solnData(GAME_VAR,ii,jj,kk) = eosData(EOS_PRES)/(eosData(EOS_EINT)*solnData(DENS_VAR,ii,jj,kk)) + 1.0e0
           enddo

           ! Initialize neutrino data
           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS ; do iE = 1, THORNADO_NE

              ioff = THORNADO_BEGIN &
                 + (iS -1)*(THORNADO_NNODESE*THORNADO_NE*THORNADO_NMOMENTS) &
                 + (iCR-1)*(THORNADO_NNODESE*THORNADO_NE) &
                 + (iE -1)*(THORNADO_NNODESE)

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

                 ivar    = ioff + iNodeE - 1

                 ! calculate actual positions of the nodes used for the gaussian quadrature
                 xnode = NodeCoordinate( MeshX(1), iX1, iNodeX1 )
                 ynode = NodeCoordinate( MeshX(2), iX2, iNodeX2 )
                 znode = NodeCoordinate( MeshX(3), iX3, iNodeX3 )

                 ss = 0.5 + 0.49 * sin( TwoPi * xnode )

                 ! J moment, iCR = 1
                 if (iCR == iCR_N) solnData(ivar,ii,jj,kk) = ss

                 ! H_x moment, iCR = 2
                 if (iCR == iCR_G1) solnData(ivar,ii,jj,kk) = ss * ( 1.0e0 - 1.0e-12 )

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
  call FreeThornado_Patch()

  return
end subroutine Simulation_initBlock
