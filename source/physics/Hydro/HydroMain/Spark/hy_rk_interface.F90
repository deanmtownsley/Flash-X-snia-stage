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
     subroutine hy_rk_getFaceFlux (limits,blkLimits,blkLimitsGC)
       implicit none
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits,blkLimits,blkLimitsGC
     end subroutine hy_rk_getFaceFlux
  end interface

  interface
     subroutine hy_rk_updateSoln (Uin,blkLimits,blklimitsGC,level,hy_del, dt, dtOld, limits, coeffs)
       implicit none
       real, pointer :: Uin(:,:,:,:)
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits, blkLimits, blkLimitsGC
       real, intent(IN) :: dt, dtOld
       real, dimension(3), intent(IN) :: coeffs,hy_del
       integer, intent(IN) :: level
     end subroutine hy_rk_updateSoln
  end interface

  interface
     subroutine hy_rk_eos(limits)
       implicit none
       integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
     end subroutine hy_rk_eos
  end interface

  interface
     subroutine hy_rk_getGravAccel(limits, blkLimitsGC, deltas)
       implicit none
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: limits, blkLimitsGC
       real,dimension(MDIM), intent(in):: deltas
     end subroutine hy_rk_getGravAccel
  end interface

end module hy_rk_interface
