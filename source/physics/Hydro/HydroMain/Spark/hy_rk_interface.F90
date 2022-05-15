!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_interface
!!
!! NAME
!!   hy_rk_interface
!!
!! SYNOPSIS
!!   use hy_rk_interface : ONLY
!!
!!  DESCRIPTION
!!    Interface for internal Spark Hydro subroutines
!!
!!***
module hy_rk_interface

#include "Simulation.h"
#include "constants.h"

  interface
     subroutine hy_rk_getFaceFlux (blockDesc,limits)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t) :: blockDesc
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
     end subroutine hy_rk_getFaceFlux
  end interface

  interface
     subroutine hy_rk_updateSoln (blockDesc, dt, dtOld, limits, coeffs)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t)    :: blockDesc
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
       real, intent(IN) :: dt, dtOld
       real, dimension(3), intent(IN) :: coeffs
     end subroutine hy_rk_updateSoln
  end interface

  interface
     subroutine hy_rk_eos(limits)
       implicit none
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
     end subroutine hy_rk_eos
  end interface

  interface
     subroutine hy_rk_getGravAccel(blockDesc,limits)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t)   :: blockDesc
       integer, intent(IN) :: limits(LOW:HIGH,MDIM)
     end subroutine hy_rk_getGravAccel
  end interface

end module hy_rk_interface
