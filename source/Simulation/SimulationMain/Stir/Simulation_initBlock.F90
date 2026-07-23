!!****if* source/Simulation/SimulationMain/Stir/Simulation_initBlock
!!
!! NAME
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            Grid_tile_t(IN) :: tileDesc)
!!
!! DESCRIPTION
!!  Initializes data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!! ARGUMENTS
!!  solnData - pointer to solution data
!!  tileDesc - describes the tile or block to initialize
!!
!! AUTHOR
!!   Christoph Federrath, 2008-2023
!!
!!***

subroutine Simulation_initBlock(solnData, tileDesc)

  use Simulation_data, ONLY: sim_dens, sim_cs, sim_magz
  use Grid_tile, ONLY : Grid_tile_t
  use Driver_data, ONLY : dr_globalMe
  use Eos_interface, ONLY : Eos_multiDim

  implicit none 

#include "constants.h"
#include "Simulation.h"

  real,              pointer    :: solnData(:,:,:,:)
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: i, j, k

  logical, parameter :: Debug = .false.

  if (Debug) print *, '[', dr_globalMe, '] Simulation_initBlock entering.'

   ! loop over cells in block
   do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
      do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
        do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)

          solnData(DENS_VAR,i,j,k) = sim_dens
          solnData(VELX_VAR,i,j,k) = 0.0
          solnData(VELY_VAR,i,j,k) = 0.0
          solnData(VELZ_VAR,i,j,k) = 0.0
#ifdef PRES_VAR
          solnData(PRES_VAR,i,j,k) = solnData(DENS_VAR,i,j,k)*sim_cs**2
#endif

#ifdef MAGX_VAR
          solnData(MAGX_VAR,i,j,k) = 0.0
#endif
#ifdef MAGY_VAR
          solnData(MAGY_VAR,i,j,k) = 0.0
#endif
#ifdef MAGZ_VAR
          solnData(MAGZ_VAR,i,j,k) = sim_magz
#endif
#ifdef MAGP_VAR
          solnData(MAGP_VAR,i,j,k) = 0.5 * (solnData(MAGX_VAR,i,j,k)**2 + &
                                            solnData(MAGY_VAR,i,j,k)**2 + &
                                            solnData(MAGZ_VAR,i,j,k)**2   )
#endif

        enddo ! i
      enddo ! j
   enddo ! k

  ! update temp, eint, ener, etc.. by calling EOS
  call Eos_multiDim(MODE_DENS_PRES, tileDesc%limits, tileDesc%blkLimitsGC(LOW,:), solnData)

  if (Debug) print *, '[', dr_globalMe, '] Simulation_initBlock exiting.'

end subroutine Simulation_initBlock
